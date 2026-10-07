"""Bounded E5 vacuole/dipeptide diagnostics; all experimental additions are temporary."""
from __future__ import annotations
import argparse
import copy
import hashlib
import math
from contextlib import contextmanager
from datetime import datetime, timezone
from importlib import metadata
import json
from pathlib import Path
import subprocess
import sys
import time

from cobra import Reaction
from cobra.io import read_sbml_model
from tools.diagnose_closed_energy import close_model, configure_solver, workspace
from tools.diagnose_dipeptide_supply import balance, sha, signature, table
from scripts.gem_annotate.energy_candidates import model_definition, protected_definitions, solver_definition
from scripts.gem_annotate.execution import execution_limits
from tools.validate_energy_candidates import SolverBudget, load_model, write, add_ntp_dissipation, energy_verdict

from scripts.gem_annotate.model_layout import PLATFORM_ROOT
ROOT = Path(__file__).resolve().parents[2]
from scripts.gem_annotate.config import load_project_paths
TASK = load_project_paths().task_outputs / 'vacuole_open_supply_20260924'
HYDRO = ('R2021', 'R2029', 'R2034', 'R2039')
CONNECTIONS = ('R1363', 'R795', 'R871', 'R876')
NAMES = ('Gly-Asp', 'Gly-Glu', 'Ala-Gly', 'Gly-Pro')
NOMINAL_C = (6, 7, 5, 7)
REPORT_REACTIONS = (*HYDRO, *CONNECTIONS, 'R2030', 'R2035', 'R2040', 'xMAINTENANCE', 'R171')


def reaction_diff(left, right):
    a, b = model_definition(left), model_definition(right)
    assert a.keys() == b.keys()
    changes = []
    for rid, fields in a['reactions'].items():
        for field, value in fields.items():
            if value != b['reactions'][rid][field]:
                changes.append({'reaction': rid, 'field': field, 'before': value,
                                'after': b['reactions'][rid][field]})
    for key in a.keys() - {'reactions'}:
        if a[key] != b[key]:
            changes.append({'reaction': '', 'field': key, 'before': a[key], 'after': b[key]})
    return changes


def dipeptides(model):
    result = {}
    for rid, name, carbon in zip(HYDRO, NAMES, NOMINAL_C):
        r = model.reactions.get_by_id(rid)
        assert r.bounds == (0, 1000) and not r.gene_reaction_rule
        assert r.get_coefficient('m1384[C_va]') == -1
        candidates = [m for m, c in r.metabolites.items() if c < 0 and m.id != 'm1384[C_va]']
        assert len(candidates) == 1 and candidates[0].compartment == 'C_va'
        met = candidates[0]
        assert r.metabolites[met] == -1
        result[rid] = {'metabolite': met.id, 'stored_name': met.name, 'sequence_name': name,
                       'formula': met.formula, 'charge': met.charge,
                       'nominal_C': carbon, 'nominal_N': 2,
                       'carbon_nitrogen_evidence': 'name-based nominal amino-acid composition, not validated model formula'}
    return result


@contextmanager
def supply_case(model, pools, supplied=(), epsilon=0.01, close=None):
    before = signature(model)
    with model:
        for rid in supplied:
            source = Reaction('DIAG_SUP_' + rid, name='Artificial direct vacuole supply of ' + pools[rid]['sequence_name'],
                              lower_bound=0, upper_bound=epsilon)
            source.add_metabolites({model.metabolites.get_by_id(pools[rid]['metabolite']): 1})
            source.notes['evidence_status'] = 'artificial material input; no endogenous production claim'
            model.add_reactions([source])
        if close:
            model.reactions.get_by_id(close).bounds = (0, 0)
        yield
    assert signature(model) == before, 'Temporary supply/bounds/constraints leaked'


class Run:
    def __init__(self, config, out, budget_path, config_path=None, input_identity=None):
        self.config, self.out = config, out
        out.mkdir(parents=True, exist_ok=False)
        self.budget = SolverBudget(budget_path, config)
        self.initial_calls = len(self.budget.record['calls'])
        self.initial_solve_seconds = self.budget.record['solve_call_seconds']
        self.started = time.monotonic()
        self.tables = {name: [] for name in ('closed_energy_results', 'growth_and_hydrolysis',
            'joint_feasibility', 'targeted_fva', 'proton_water_atp_ledger', 'ledger_totals',
            'related_fluxes', 'dipeptide_fates', 'source_carbon_nitrogen', 'monotonicity')}
        self.manifest = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
            'mode': 'fresh_validation', 'config': config,
            'config_path': str(config_path) if config_path else None,
            'config_sha256': sha(config_path) if config_path else None,
            'explicit_inputs': input_identity, 'budget_path': str(budget_path),
            'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            'software': {'python': sys.version, **{p: metadata.version(p) for p in ('cobra', 'optlang', 'gurobipy', 'memote')}},
            'source_sha256': {str(Path(m.__file__).relative_to(ROOT)): sha(m.__file__)
                for m in list(sys.modules.values()) if getattr(m, '__file__', None)
                and Path(m.__file__).is_relative_to(PLATFORM_ROOT)},
            'loaders': {}, 'outcomes': [], 'stages': [], 'primary_optima': {}}
        self.save()

    def save(self):
        self.manifest['wall_seconds'] = time.monotonic() - self.started
        self.manifest['optimization_calls'] = len(self.budget.record['calls']) - self.initial_calls
        self.manifest['solve_call_seconds'] = self.budget.record['solve_call_seconds'] - self.initial_solve_seconds
        self.manifest['aggregate_optimization_calls'] = len(self.budget.record['calls'])
        self.manifest['aggregate_solve_call_seconds'] = self.budget.record['solve_call_seconds']
        write(self.out/'manifest.json', self.manifest)
        for name, rows in self.tables.items():
            table(self.out/(name+'.tsv'), rows)

    def solve(self, model, label, purpose):
        configure_solver(model, self.config)
        model.solver.update()
        path = self.out/(label+'.json')
        row, flux = self.budget.solve(model, self.out.name+'/'+label, path, purpose)
        # Retain actual split-variable bounds/values and objective for independent audit.
        data = json.loads(path.read_text())
        native = model.solver.problem
        data['actual_solver'] = {'variables': native.VarName, 'lower_bounds': native.LB,
            'upper_bounds': native.UB, 'objective_coefficients': native.Obj,
            'objective_constant': native.ObjCon, 'objective_sense': native.ModelSense,
            'constraint_sense': native.Sense, 'constraint_rhs': native.RHS,
            'primal': native.X if flux is not None else None}
        write(path, data)
        self.manifest['outcomes'].append({'label': label, **row})
        self.save()
        return row, flux

    def require(self, result):
        if result['status'] != 'optimal':
            raise RuntimeError('Expected feasible control unresolved: '+result['label'])

    def loaded(self, path, key):
        raw = read_sbml_model(path)
        locks = protected_definitions(raw)
        sim = load_model(path, self.config)
        assert protected_definitions(sim.model) == locks
        assert model_definition(sim.model)['objective'] == {'biomass_C': 1.0}
        assert sim.model.reactions.xMAINTENANCE.bounds == (7.8625, 1000)
        assert 'R1219' not in sim.active_medium
        for rid in CONNECTIONS:
            assert sim.model.reactions.get_by_id(rid).bounds == raw.reactions.get_by_id(rid).bounds
        self.manifest['loaders'][key] = {'path': str(path.relative_to(ROOT)), 'sha256': sha(path),
            'provenance': sim.provenance(), 'active_medium': sim.active_medium,
            'strain_overlay': sim.strain_overlay_audit, 'file_to_loaded_diff': reaction_diff(raw, sim.model),
            'energy_locks_preserved': True}
        write(self.out/(key+'_effective.json'), signature(sim.model))
        self.save()
        return sim.model

    def ledger(self, model, flux, label, pools):
        for rid in REPORT_REACTIONS:
            r = model.reactions.get_by_id(rid)
            self.tables['related_fluxes'].append({'witness': label, 'reaction': rid,
                'flux': flux[rid], 'equation': r.reaction, 'kind': 'requested'})
        for r in model.boundary:
            if r.id.startswith('DIAG_SUP_') or abs(flux[r.id]) > self.config['tolerance']:
                self.tables['related_fluxes'].append({'witness': label, 'reaction': r.id,
                    'flux': flux[r.id], 'equation': r.reaction, 'kind': 'actual_material_boundary'})
        va_amino = {m.id for rid in HYDRO for m, c in model.reactions.get_by_id(rid).metabolites.items() if c > 0}
        cy_amino = set()
        for rid in ('R871', 'R876', 'R2030', 'R2035', 'R2040'):
            cy_amino.update(m.id for m in model.reactions.get_by_id(rid).metabolites
                            if m.compartment == 'C_cy' and m.id != 'm10[C_cy]')
        mids = va_amino | cy_amino | {'m1384[C_va]', 'm1007[C_va]', 'm141[C_cy]'}
        for mid in sorted(mids):
            met = model.metabolites.get_by_id(mid)
            rows = []
            for r in sorted(met.reactions, key=lambda r: r.id):
                contribution = r.metabolites[met]*flux[r.id]
                rows.append({'witness': label, 'metabolite': mid, 'name': met.name,
                    'compartment': met.compartment, 'formula': met.formula, 'charge': met.charge,
                    'reaction': r.id, 'coefficient': r.metabolites[met], 'flux': flux[r.id],
                    'contribution': contribution, 'equation': r.reaction})
            self.tables['proton_water_atp_ledger'].extend(rows)
            self.tables['ledger_totals'].append({'witness': label, 'metabolite': mid, 'name': met.name,
                'production': sum(max(x['contribution'], 0) for x in rows),
                'consumption': -sum(min(x['contribution'], 0) for x in rows),
                'residual': sum(x['contribution'] for x in rows)})
        for rid, pool in pools.items():
            met = model.metabolites.get_by_id(pool['metabolite'])
            for r in sorted(met.reactions, key=lambda r:r.id):
                self.tables['dipeptide_fates'].append({'witness': label, 'dipeptide': pool['sequence_name'],
                    'metabolite': met.id, 'reaction': r.id, 'flux': flux[r.id],
                    'pool_contribution': r.metabolites[met]*flux[r.id], 'equation': r.reaction})
            source = 'DIAG_SUP_'+rid
            value = flux.get(source, 0)
            self.tables['source_carbon_nitrogen'].append({'witness': label, 'source': source,
                'sequence_name': pool['sequence_name'], 'source_flux': value, 'formula': pool['formula'],
                'chemistry_status': 'unverifiable_missing_formula' if not pool['formula'] else 'stored_formula_present',
                'nominal_C_input': value*pool['nominal_C'], 'nominal_N_input': value*pool['nominal_N'],
                'evidence': pool['carbon_nitrogen_evidence']})
        self.save()

    def growth_row(self, model, label, mode, epsilon, supplied, closed, result, flux, pools):
        self.require(result)
        row = {'scenario': label, 'objective_kind': mode, 'epsilon': epsilon, 'supplied': list(supplied),
            'closed_connection': closed or '', 'status': result['status'],
            'objective_value': result['objective'], 'biomass_C': flux['biomass_C'],
            'xMAINTENANCE': flux['xMAINTENANCE'], **{rid: flux[rid] for rid in REPORT_REACTIONS if rid!='xMAINTENANCE'},
            **{'source_'+rid: flux.get('DIAG_SUP_'+rid, 0) for rid in HYDRO},
            'actual_boundary_fluxes': {r.id: flux[r.id] for r in model.boundary if abs(flux[r.id])>self.config['tolerance']}}
        self.tables['growth_and_hydrolysis'].append(row)
        self.save()
        if mode in ('growth', 'pfba'):
            self.ledger(model, flux, label+'_'+mode, pools)

    def energies(self, base, candidate):
        for name, model in [('C0', base), ('C1', candidate), ('C2', candidate)]:
            before = signature(model)
            with model:
                if name == 'C2':
                    for rid in CONNECTIONS:model.reactions.get_by_id(rid).bounds=(0,1000)
                closed, changes = close_model(model, self.config)
                assert balance(closed.reactions.xMAINTENANCE)['element_status'] == 'balanced'
                assert balance(closed.reactions.xMAINTENANCE)['charge_status'] == 'balanced_as_stored'
                self.manifest.setdefault('closures', {})[name] = changes
                write(self.out/(name+'_closed_definition.json'), signature(closed))
                closed_before = signature(closed)
                with closed:
                    for r in closed.reactions:r.bounds=(0,0)
                    closed.objective=closed.problem.Objective(0)
                    row,flux=self.solve(closed,name+'_zero',{'kind':'closed_exact_zero','scenario':name})
                    self.require(row)
                    assert all(abs(v)<=self.config['tolerance'] for v in flux.values())
                    self.tables['closed_energy_results'].append({'scenario':name,'carrier':'zero_feasibility',
                        'status':row['status'],'maximum':row['objective'],'verdict':'feasible'})
                for ntp in ('ATP','GTP','UTP','CTP'):
                    with closed:
                        if ntp != 'ATP':add_ntp_dissipation(closed,ntp)
                        else:closed.objective=closed.reactions.xMAINTENANCE
                        row,flux=self.solve(closed,name+'_'+ntp,{'kind':'closed_energy_max','scenario':name,'carrier':ntp})
                        verdict=energy_verdict(row['status'],row['objective'],self.config['tolerance'])
                        self.tables['closed_energy_results'].append({'scenario':name,'carrier':ntp,
                            'status':row['status'],'maximum':row['objective'],'verdict':verdict})
                        self.save()
                        if verdict != 'within_tolerance':
                            if flux and row['objective']>self.config['tolerance']:
                                drain=closed.reactions.get_by_id('xMAINTENANCE' if ntp=='ATP' else 'DIAG_D_'+ntp)
                                with closed:
                                    demand=min(1.,row['objective']/2)
                                    drain.bounds=(demand,demand)
                                    closed.objective=closed.problem.Objective(sum(r.forward_variable+r.reverse_variable for r in closed.reactions),direction='min')
                                    self.solve(closed,name+'_'+ntp+'_witness',{'kind':'positive_energy_L1_witness','demand':demand})
                                for rid in CONNECTIONS:
                                    with closed:
                                        closed.reactions.get_by_id(rid).bounds=(0,0)
                                        self.solve(closed,name+'_'+ntp+'_close_'+rid,{'kind':'energy_connection_isolation','closed':rid})
                            raise RuntimeError('Stop candidate acceptance: '+name+' '+ntp+' '+verdict)
                assert signature(closed)==closed_before
            assert signature(model)==before
        self.manifest['stages'].append('closed_energy_passed');self.save()

    def primary(self, base, candidate, pools):
        cases=[('G0',base,(),None),('G1',candidate,(),None),('G2',base,HYDRO,None),('G3',candidate,HYDRO,None)]
        cases += [('G'+str(i+4),candidate,(rid,),None) for i,rid in enumerate(HYDRO)]
        cases += [('G'+str(i+8),candidate,HYDRO,rid) for i,rid in enumerate(CONNECTIONS)]
        for label,model,supplied,closed in cases:
            with supply_case(model,pools,supplied,self.config['main_epsilon'],closed):
                model.objective=model.reactions.biomass_C
                row,flux=self.solve(model,label+'_growth',{'kind':'culture_growth','scenario':label,'sources':list(supplied),'closed':closed})
                self.growth_row(model,label,'growth',self.config['main_epsilon'],supplied,closed,row,flux,pools)
                self.manifest['primary_optima'][label]=flux['biomass_C']
        mu=self.manifest['primary_optima'];tol=self.config['tolerance']
        for high,low in [('G1','G0'),('G2','G0'),('G3','G1')]+[('G'+str(i),'G1') for i in range(4,8)]+[('G3','G'+str(i)) for i in range(4,12)]:
            ok=mu[high]>=mu[low]-tol
            self.tables['monotonicity'].append({'relaxed_scenario':high,'restricted_scenario':low,'difference':mu[high]-mu[low],'passes':ok})
            if not ok:raise RuntimeError('Relaxation decreased optimum: '+high+' vs '+low)
        self.manifest['stages'].append('primary_growth_complete');self.save()
        return cases

    def joint(self, candidate, pools):
        mu=self.manifest['primary_optima']['G3'];eps=self.config['main_epsilon']
        floors={'none':None,'99pct':self.config['fva_fraction']*mu,'near_strict':mu-self.config['strict_growth_delta']}
        for requirement,floor in floors.items():
            for fraction in (1.,.5):
                with supply_case(candidate,pools,HYDRO,eps):
                    for rid in HYDRO:candidate.reactions.get_by_id(rid).lower_bound=eps*fraction
                    if floor is not None:candidate.add_cons_vars(candidate.problem.Constraint(candidate.reactions.biomass_C.flux_expression,lb=floor,name='DIAG_joint_growth'))
                    candidate.objective=candidate.reactions.biomass_C
                    label='G3_joint_'+requirement+('_full' if fraction==1 else '_half')
                    row,flux=self.solve(candidate,label,{'kind':'joint_hydrolysis','growth_requirement':requirement,'growth_floor':floor,'hydrolysis_minimum':eps*fraction})
                    self.tables['joint_feasibility'].append({'case':label,'growth_requirement':requirement,'growth_floor':floor,'hydrolysis_minimum':eps*fraction,
                        'status':row['status'],'biomass_C':flux['biomass_C'] if flux else None,
                        'hydrolysis_fluxes':{rid:flux[rid] for rid in HYDRO} if flux else None,
                        'R795':flux['R795'] if flux else None,'R1363':flux['R1363'] if flux else None})
                    if flux:self.ledger(candidate,flux,label,pools)
                    self.save()
                if row['status']=='optimal':break
        self.manifest['stages'].append('joint_feasibility_complete');self.save()

    def usage(self, cases, pools):
        for label,model,supplied,closed in cases:
            with supply_case(model,pools,supplied,self.config['main_epsilon'],closed):
                with model:
                    mu=self.manifest['primary_optima'][label]
                    model.add_cons_vars(model.problem.Constraint(model.reactions.biomass_C.flux_expression,lb=mu-self.config['pfba_growth_slack'],name='DIAG_pfba_growth'))
                    model.objective=model.problem.Objective(sum(r.forward_variable+r.reverse_variable for r in model.reactions),direction='min')
                    row,flux=self.solve(model,label+'_pfba',{'kind':'pfba_secondary','growth_floor':mu,'primary_optimum':mu})
                    self.growth_row(model,label,'pfba',self.config['main_epsilon'],supplied,closed,row,flux,pools)
                for rid in HYDRO:
                    with model:
                        model.objective=model.reactions.get_by_id(rid);model.objective.direction='max'
                        row,flux=self.solve(model,label+'_max_'+rid,{'kind':'individual_hydrolysis_max','hydrolase':rid,'growth_requirement':None})
                        self.growth_row(model,label,'max_'+rid,self.config['main_epsilon'],supplied,closed,row,flux,pools)
        self.manifest['stages'].append('pfba_and_individual_hydrolysis_complete');self.save()

    def fva(self, candidate, pools):
        mu=self.manifest['primary_optima']['G3']
        for requirement,floor in [('99pct',.99*mu),('near_strict',mu-self.config['strict_growth_delta'])]:
            with supply_case(candidate,pools,HYDRO,self.config['main_epsilon']):
                candidate.add_cons_vars(candidate.problem.Constraint(candidate.reactions.biomass_C.flux_expression,lb=floor,name='DIAG_fva_growth'))
                for rid in (*HYDRO,'R795'):
                    for direction in ('min','max'):
                        with candidate:
                            candidate.objective=candidate.reactions.get_by_id(rid);candidate.objective.direction=direction
                            row,flux=self.solve(candidate,'G3_FVA_'+requirement+'_'+rid+'_'+direction,{'kind':'targeted_fva','reaction':rid,'direction':direction,'growth_floor':floor})
                            self.require(row)
                            self.tables['targeted_fva'].append({'scenario':'G3','requirement':requirement,'growth_floor':floor,'delta':self.config['strict_growth_delta'] if requirement=='near_strict' else None,
                                'reaction':rid,'direction':direction,'status':row['status'],'flux':flux[rid],'biomass_C':flux['biomass_C']})
                            self.save()
        self.manifest['stages'].append('targeted_fva_complete');self.save()

    def sensitivity(self, candidate, pools):
        previous=self.manifest['primary_optima']['G1']
        for eps in self.config['sensitivity_epsilon']:
            label='G3_eps_'+str(eps)
            with supply_case(candidate,pools,HYDRO,eps):
                candidate.objective=candidate.reactions.biomass_C
                row,flux=self.solve(candidate,label+'_growth',{'kind':'supply_sensitivity_growth','epsilon':eps})
                self.growth_row(candidate,label,'growth',eps,HYDRO,None,row,flux,pools)
                mu=flux['biomass_C']
                assert mu>=previous-self.config['tolerance'] and mu<=self.manifest['primary_optima']['G3']+self.config['tolerance']
                previous=mu
                with candidate:
                    candidate.add_cons_vars(candidate.problem.Constraint(candidate.reactions.biomass_C.flux_expression,lb=mu,name='DIAG_pfba_growth'))
                    candidate.objective=candidate.problem.Objective(sum(r.forward_variable+r.reverse_variable for r in candidate.reactions),direction='min')
                    row,flux=self.solve(candidate,label+'_pfba',{'kind':'supply_sensitivity_pfba','epsilon':eps,'growth_floor':mu})
                    self.growth_row(candidate,label,'pfba',eps,HYDRO,None,row,flux,pools)
        self.manifest['stages'].append('sensitivity_complete');self.save()



def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def input_identity(path, expected_sha):
    """Bind the requested path, never another candidate with a matching name."""
    path = workspace(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    if not isinstance(expected_sha, str) or len(expected_sha) != 64:
        raise ValueError('An explicit complete expected SHA256 is required: ' + str(path))
    actual = sha(path)
    if actual != expected_sha:
        raise ValueError('Input SHA256 differs: ' + str(path))
    return {'path': str(path), 'sha256': actual}


def bind_inputs(baseline_model, baseline_sha256, candidate_model, candidate_sha256=None, build_manifest=None):
    baseline = input_identity(baseline_model, baseline_sha256)
    candidate_path = workspace(candidate_model)
    manifest_record = None
    if build_manifest is not None:
        manifest_path = workspace(build_manifest)
        record = json.loads(manifest_path.read_text())
        if workspace(record['output']) != candidate_path:
            raise ValueError('Build manifest belongs to a different candidate path')
        if (workspace(record['source']) != Path(baseline['path'])
                or record['source_sha256'] != baseline['sha256']):
            raise ValueError('Build manifest belongs to a different baseline')
        if candidate_sha256 is not None and candidate_sha256 != record['output_sha256']:
            raise ValueError('Candidate SHA and build manifest disagree')
        candidate_sha256 = record['output_sha256']
        manifest_record = {'path': str(manifest_path), 'sha256': sha(manifest_path),
                           'source_and_output_bind_requested_paths': True}
    candidate = input_identity(candidate_path, candidate_sha256)
    return {'baseline': baseline, 'candidate': candidate, 'build_manifest': manifest_record}


def math_signature(model):
    """Actual LP definition, without chemical names, formulae, or charges."""
    return solver_definition(model)


def preflight_models(base, candidate, config, expected_diff):
    if any(r.id.startswith('DIAG_') for model in (base, candidate) for r in model.reactions):
        raise ValueError('Input XML contains diagnostic reactions')
    if protected_definitions(base) != protected_definitions(candidate) or len(protected_definitions(base)) != 8:
        raise ValueError('Eight E5 energy protection definitions must match')
    diff = reaction_diff(base, candidate)
    if expected_diff == 'four_connections':
        expected = [{'reaction': rid, 'field': 'bounds', 'before': [0., 0.],
                     'after': config['connection_bounds'][rid]} for rid in CONNECTIONS]
        if sorted(diff, key=lambda x: x['reaction']) != sorted(expected, key=lambda x: x['reaction']):
            raise ValueError('Candidate changes are not exactly the four reviewed connection bounds')
    elif expected_diff == 'chemistry_only':
        if math_signature(base) != math_signature(candidate):
            raise ValueError('Chemical metadata patch changes the actual optimization problem')
        a, b = model_definition(base), model_definition(candidate)
        targets = {'m1869[C_ex]', 'm1870[C_cy]', 'm1871[C_va]',
                   'm1860[C_ex]', 'm1861[C_cy]', 'm1862[C_va]',
                   'm1876[C_ex]', 'm1877[C_cy]', 'm1878[C_va]',
                   'm1864[C_ex]', 'm1865[C_cy]', 'm1866[C_va]', 'm1863[C_va]'}
        if a['reactions'] != b['reactions'] or a['metabolites'].keys() != b['metabolites'].keys():
            raise ValueError('Chemical candidate changes reactions/GPR or metabolite identifiers')
        for mid, value in a['metabolites'].items():
            other = b['metabolites'][mid]
            if (mid not in targets and value != other) or value['compartment'] != other['compartment']:
                raise ValueError('Chemical candidate changes an unauthorized metabolite field: ' + mid)
    else:
        raise ValueError('An explicit expected-diff kind is required')
    return {'expected_diff': expected_diff, 'definition_diff': diff,
            'mathematical_signature_sha256': {name: digest(math_signature(model))
                for name, model in [('baseline', base), ('candidate', candidate)]},
            'chemical_metadata_sha256': {name: digest(model_definition(model)['metabolites'])
                for name, model in [('baseline', base), ('candidate', candidate)]},
            'counts': {name: {'reactions': len(model.reactions), 'metabolites': len(model.metabolites)}
                for name, model in [('baseline', base), ('candidate', candidate)]},
            'energy_locks_match': True}


def incremental_validation(runner, base, candidate, pools):
    """Eleven targeted solves; no repeat of the historical 114-solve matrix."""
    starting = {'baseline': signature(base), 'candidate': signature(candidate)}
    closed, changes = close_model(candidate, runner.config)
    runner.manifest['closures'] = {'candidate': changes}
    closed_before = signature(closed)
    with closed:
        for reaction in closed.reactions:
            reaction.bounds = (0, 0)
        closed.objective = closed.problem.Objective(0)
        row, flux = runner.solve(closed, 'candidate_zero', {'kind': 'closed_exact_zero'})
        runner.require(row)
        if any(abs(value) > runner.config['tolerance'] for value in flux.values()):
            raise ValueError('Closed zero control has nonzero flux')
        runner.tables['closed_energy_results'].append({'scenario': 'candidate', 'carrier': 'zero_feasibility',
            'status': row['status'], 'maximum': row['objective'], 'verdict': 'feasible'})
    for ntp in ('ATP', 'GTP', 'UTP', 'CTP'):
        with closed:
            if ntp == 'ATP':
                closed.objective = closed.reactions.xMAINTENANCE
            else:
                add_ntp_dissipation(closed, ntp)
            row, flux = runner.solve(closed, 'candidate_' + ntp,
                                     {'kind': 'closed_energy_max', 'carrier': ntp})
            verdict = energy_verdict(row['status'], row['objective'], runner.config['tolerance'])
            runner.tables['closed_energy_results'].append({'scenario': 'candidate', 'carrier': ntp,
                'status': row['status'], 'maximum': row['objective'], 'verdict': verdict})
            runner.save()
            if verdict != 'within_tolerance':
                raise RuntimeError('Stop acceptance: ' + ntp + ' ' + verdict)
    if signature(closed) != closed_before:
        raise ValueError('Closed diagnostic constraints leaked')
    for label, model, supplied in [('baseline_no_source', base, ()), ('candidate_no_source', candidate, ()),
                                    ('G3', candidate, HYDRO)]:
        with supply_case(model, pools, supplied, runner.config['main_epsilon']):
            model.objective = model.reactions.biomass_C
            row, flux = runner.solve(model, label + '_growth', {'kind': 'culture_growth', 'sources': list(supplied)})
            runner.growth_row(model, label, 'growth', runner.config['main_epsilon'], supplied, None, row, flux, pools)
            runner.manifest['primary_optima'][label] = flux['biomass_C']
    mu = runner.manifest['primary_optima']['G3']
    with supply_case(candidate, pools, HYDRO, runner.config['main_epsilon']):
        candidate.add_cons_vars(candidate.problem.Constraint(candidate.reactions.biomass_C.flux_expression,
            lb=mu-runner.config['pfba_growth_slack'], name='DIAG_pfba_growth'))
        candidate.objective = candidate.problem.Objective(sum(r.forward_variable+r.reverse_variable
            for r in candidate.reactions), direction='min')
        row, flux = runner.solve(candidate, 'G3_pfba', {'kind': 'pfba_secondary',
            'primary_optimum': mu, 'growth_floor': mu-runner.config['pfba_growth_slack']})
        runner.growth_row(candidate, 'G3', 'pfba', runner.config['main_epsilon'], HYDRO, None, row, flux, pools)
    with supply_case(candidate, pools, HYDRO, runner.config['main_epsilon']):
        for rid in HYDRO:
            candidate.reactions.get_by_id(rid).lower_bound = runner.config['main_epsilon']
        row, flux = runner.solve(candidate, 'G3_joint', {'kind': 'joint_hydrolysis',
            'hydrolysis_minimum': runner.config['main_epsilon'], 'growth_requirement': None})
        runner.require(row)
        runner.tables['joint_feasibility'].append({'case': 'G3_joint', 'status': row['status'],
            'biomass_C': flux['biomass_C'], 'hydrolysis_fluxes': {rid: flux[rid] for rid in HYDRO},
            'R795': flux['R795'], 'R1363': flux['R1363']})
        runner.ledger(candidate, flux, 'G3_joint', pools)
    if starting != {'baseline': signature(base), 'candidate': signature(candidate)}:
        raise ValueError('Temporary sources or constraints leaked')
    row, flux = runner.solve(candidate, 'candidate_restored_growth', {'kind': 'restoration_growth'})
    runner.growth_row(candidate, 'candidate_restored', 'growth', 0., (), None, row, flux, pools)
    optima = runner.manifest['primary_optima']
    tol = runner.config['tolerance']
    if abs(flux['biomass_C'] - optima['candidate_no_source']) > tol:
        raise ValueError('Restored growth differs from the same candidate before supply')
    if optima['G3'] < optima['candidate_no_source'] - tol:
        raise ValueError('Adding optional supply decreased growth')
    if (runner.manifest['preflight']['expected_diff'] == 'chemistry_only'
            and abs(optima['baseline_no_source']-optima['candidate_no_source']) > tol):
        raise ValueError('Mathematically identical baseline and chemical candidate differ in growth')
    runner.manifest.update(all_temporary_changes_removed=True, energy_locks_match_E5=True)
    runner.manifest['stages'].append('incremental_fresh_validation_complete')


def audit_saved_run(saved_run, source_archive, out, budget_path, config_path, identities):
    """Audit historical evidence against its archived code, never call optimize."""
    saved_run, source_archive, out = map(workspace, (saved_run, source_archive, out))
    if out.exists():
        raise FileExistsError(out)
    manifest = json.loads((saved_run/'manifest.json').read_text())
    budget = json.loads(workspace(budget_path).read_text())
    config_path = workspace(config_path)
    if manifest['status'] != 'complete' or manifest['config_sha256'] != sha(config_path):
        raise ValueError('Saved run incomplete or configuration mismatch')
    for role, saved_key in [('baseline', 'E5'), ('candidate', 'candidate')]:
        loader = manifest['loaders'][saved_key]
        if (workspace(loader['path']) != Path(identities[role]['path'])
                or loader['sha256'] != identities[role]['sha256']):
            raise ValueError('Saved run used a different ' + role)
    for path, expected in manifest['config']['input_configuration_sha256'].items():
        input_identity(path, expected)
    source_rows = []
    for path, expected in manifest['source_sha256'].items():
        archived = source_archive / path
        if not archived.is_file() or sha(archived) != expected:
            raise ValueError('Archived executed source missing or changed: ' + path)
        current = workspace(path)
        source_rows.append({'path': path, 'executed_sha256': expected, 'archived_sha256': sha(archived),
            'current_sha256': sha(current) if current.is_file() else None,
            'current_matches_executed': current.is_file() and sha(current) == expected})
    if (len(budget['calls']) != manifest['optimization_calls']
            or len(budget['calls']) != len(manifest['outcomes'])):
        raise ValueError('Saved call count differs')
    if abs(sum(call['solve_call_seconds'] for call in budget['calls'])-budget['solve_call_seconds']) > 1e-8:
        raise ValueError('Saved timing ledger differs')
    for outcome, call in zip(manifest['outcomes'], budget['calls']):
        if call['status'] not in ('optimal', 'infeasible') or call.get('postprocessing') != 'complete':
            raise ValueError('Unresolved historical solve')
        record = json.loads((saved_run/(outcome['label'].rsplit('/', 1)[-1] + '.json')).read_text())
        if record['result'] != call or record['result']['status'] != outcome['status']:
            raise ValueError('Saved raw solve differs from budget/outcome: ' + outcome['label'])
        if call['status'] == 'optimal' and (record['fluxes'] is None
                or not all(math.isfinite(v) for v in record['fluxes'].values())):
            raise ValueError('Missing or nonfinite saved optimal flux')
    result = {'mode': 'audit_saved_run', 'status': 'complete', 'optimization_calls_this_audit': 0,
              'saved_optimization_calls': len(budget['calls']), 'saved_run': str(saved_run),
              'saved_manifest_sha256': sha(saved_run/'manifest.json'), 'explicit_inputs': identities,
              'source_identity': source_rows,
              'scope': 'input/code archive, call/status/timing and raw-result identity audit; no new optimization'}
    out.mkdir(parents=True)
    write(out/'manifest.json', result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', required=True, choices=('audit_saved_run', 'fresh_validation'))
    parser.add_argument('--config', required=True)
    parser.add_argument('--baseline-model', required=True)
    parser.add_argument('--baseline-sha256', required=True)
    parser.add_argument('--candidate-model', required=True)
    parser.add_argument('--candidate-sha256')
    parser.add_argument('--build-manifest')
    parser.add_argument('--expected-diff', choices=('four_connections', 'chemistry_only'))
    parser.add_argument('--output', required=True)
    parser.add_argument('--budget', required=True)
    parser.add_argument('--saved-run')
    parser.add_argument('--source-archive')
    args = parser.parse_args(argv)
    if not args.candidate_sha256 and not args.build_manifest:
        parser.error('--candidate-sha256 or --build-manifest is required')
    if args.mode == 'fresh_validation' and not args.expected_diff:
        parser.error('--expected-diff is required for fresh_validation')
    if args.mode == 'audit_saved_run' and (not args.saved_run or not args.source_archive):
        parser.error('--saved-run and --source-archive are required for audit_saved_run')
    with execution_limits(no_solve=args.mode == 'audit_saved_run', allow_network=False):
        identities = bind_inputs(args.baseline_model, args.baseline_sha256, args.candidate_model,
                                 args.candidate_sha256, args.build_manifest)
        config_path = workspace(args.config)
        if args.mode == 'audit_saved_run':
            audit_saved_run(args.saved_run, args.source_archive, args.output, args.budget, config_path, identities)
            return
        config = json.loads(config_path.read_text())
        for path, expected in config['input_configuration_sha256'].items():
            input_identity(path, expected)
        config = copy.deepcopy(config)
        config.update(model=identities['baseline']['path'], model_sha256=identities['baseline']['sha256'],
                      candidate=identities['candidate']['path'], candidate_sha256=identities['candidate']['sha256'])
        raw_a = read_sbml_model(config['model'])
        raw_b = read_sbml_model(config['candidate'])
        preflight = preflight_models(raw_a, raw_b, config, args.expected_diff)
        runner = Run(config, workspace(args.output), workspace(args.budget), config_path, identities)
        try:
            runner.manifest['preflight'] = preflight
            base = runner.loaded(workspace(config['model']), 'baseline')
            candidate = runner.loaded(workspace(config['candidate']), 'candidate')
            runner.manifest['loaded_preflight'] = preflight_models(base, candidate, config, args.expected_diff)
            pools = dipeptides(candidate)
            write(runner.out/'dipeptide_identities.json', pools)
            runner.manifest['candidate_sha256'] = identities['candidate']['sha256']
            incremental_validation(runner, base, candidate, pools)
            if bind_inputs(args.baseline_model, args.baseline_sha256, args.candidate_model,
                           args.candidate_sha256, args.build_manifest) != identities:
                raise ValueError('Input changed during validation')
            runner.manifest['status'] = 'complete'
        except Exception as exc:
            runner.manifest.update(status='incomplete', error=repr(exc))
            raise
        finally:
            runner.save()


if __name__ == '__main__':
    main()
