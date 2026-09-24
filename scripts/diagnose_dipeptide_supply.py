"""Bounded, isolated dipeptide supply diagnostics; never exports a changed model."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from importlib import metadata
import inspect
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from unittest.mock import patch

import cobra
import pandas as pd
from cobra import Reaction
from cobra.util.solver import linear_reaction_coefficients
from memote.support import consistency, helpers

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

ROOT = Path(__file__).resolve().parents[1]
MODEL = ROOT / 'model_metadata_trna_r1159_leak.xml'
CONDITIONS = ROOT / 'artifacts/reference_pipeline_restore_20260909/research/state'
TARGETS = {
    'R2021': ('m1871[C_va]', 'Gly-Asp', 6),
    'R2029': ('m1862[C_va]', 'Gly-Glu', 7),
    'R2034': ('m1878[C_va]', 'Ala-Gly', 5),
    'R2039': ('m1866[C_va]', 'Gly-Pro', 7),
}
WATER = 'm1384[C_va]'
LEVELS = (1e-4, 1e-3, 1e-2)
TOL = 1e-7
GROWTH_SLACK = 1e-8
ENDO_REQUIRED = (
    'precursor_identity_and_sequence', 'precursor_formula_charge',
    'residue_composition', 'source_synthesis_reactions', 'synthesis_energy_cost',
    'degradation_rate_mmol_precursor_gdw_h', 'vacuolar_fraction',
    'nonoverlapping_dipeptide_yields', 'remaining_products',
    'balanced_degradation_equation', 'localization_and_transport_evidence',
    'net_growth_vs_replacement_accounting', 'condition_matched_parameter_sources',
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def table(path, rows):
    rows = list(rows)
    if not rows:
        return
    with Path(path).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in row.items()})


def formula_elements(formula):
    """Reject unknown/generic formulae instead of treating missing atoms as zero."""
    if not formula or not re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', formula):
        return None
    pairs = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    if any(e not in {'C', 'H', 'O', 'N', 'P', 'S', 'Fe', 'Mg', 'Zn', 'Na', 'K', 'Cl', 'Ca', 'Cu', 'Mn', 'Co', 'Mo', 'Se'} for e, _ in pairs):
        return None
    result = {}
    for element, n in pairs:
        result[element] = result.get(element, 0) + int(n or 1)
    return result


def balance(reaction):
    missing = [m.id for m in reaction.metabolites if formula_elements(m.formula) is None]
    elements = {}
    if not missing:
        for m, c in reaction.metabolites.items():
            for e, n in formula_elements(m.formula).items():
                elements[e] = elements.get(e, 0) + c * n
    charge = None if any(m.charge is None for m in reaction.metabolites) else sum(c * m.charge for m, c in reaction.metabolites.items())
    return dict(element_status='unverifiable' if missing else ('balanced' if all(abs(v) < 1e-9 for v in elements.values()) else 'imbalanced'),
                missing_formula=missing, element_residual=elements if not missing else None,
                stored_charge_residual=charge, charge_status='unknown' if charge is None else ('balanced_as_stored' if abs(charge) < 1e-9 else 'imbalanced'))


def signature(model):
    # Native sparse coefficients avoid repeatedly rendering 1877 symbolic rows.
    # Names and ordering are included, so equal matrices cannot hide renamed rows.
    model.solver.update()
    native = model.solver.problem
    matrix = native.getA().tocsr()
    return {
        'reactions': {r.id: [r.lower_bound, r.upper_bound, r.gene_reaction_rule,
                           {m.id: c for m, c in r.metabolites.items()}] for r in model.reactions},
        'metabolites': {m.id: [m.name, m.formula, m.charge, m.compartment] for m in model.metabolites},
        'objective': str(model.objective.expression), 'direction': model.objective.direction,
        'constraints': {'names': native.ConstrName, 'bounds': [[c.lb,c.ub] for c in model.constraints],
                        'variables': native.VarName, 'matrix_indptr':matrix.indptr.tolist(),
                        'matrix_indices':matrix.indices.tolist(), 'matrix_data':matrix.data.tolist()},
    }


def validate_endogenous_spec(spec):
    if not isinstance(spec.get('enabled'), bool):
        raise ValueError('enabled must be explicitly true or false')
    if spec.get('enabled') is not True:
        return {'enabled': False, 'added_reactions': [], 'missing_parameters': [k for k in ENDO_REQUIRED if spec.get(k) is None]}
    missing = [k for k in ENDO_REQUIRED if spec.get(k) is None]
    if missing:
        raise ValueError('Endogenous module has no accepted parameterization: ' + ', '.join(missing))
    raise ValueError('Specification only: complete fields still require source, chemistry, cost and flux validation before implementation')


def cached_solution(row, witness, extra):
    """Reuse a completed label only under exactly the same declared constraints."""
    if any(row.get(k) != v for k,v in (extra or {}).items()):
        raise ValueError('Saved solve does not match requested scenario/constraints')
    if row['status']=='optimal' and (not row['numerically_valid'] or witness is None):
        raise ValueError('Missing valid saved flux witness')
    if witness is not None and witness.get('extra') != extra:
        raise ValueError('Saved witness context mismatch')
    return row, cobra.Solution(row['objective'],row['status'],pd.Series(witness['fluxes'] if witness else {},dtype=float))


def add_boundary(model, rid, mid, coefficient, cap, purpose):
    if not math.isfinite(cap) or not 0 < cap <= 0.1:
        raise ValueError('Diagnostic boundary cap must be finite, positive and <=0.1')
    if rid in model.reactions:
        raise ValueError('Diagnostic ID already exists: ' + rid)
    r = Reaction(rid, name=purpose, lower_bound=0, upper_bound=cap)
    r.add_metabolites({model.metabolites.get_by_id(mid): coefficient})
    r.notes = {'evidence': 'Artificial diagnostic boundary; no biological supply or native transporter claim'}
    model.add_reactions([r])


@contextmanager
def scenario(model, supplied=(), epsilon=0.0, water=False, drains=False, pump=False, outputs=False):
    """All mutations, including constraints/objectives added by callers, roll back."""
    before = signature(model)
    with model:
        changes = {'sources': {}, 'drains': {}, 'bounds': {}}
        for rid in supplied:
            mid = TARGETS[rid][0]
            r = model.reactions.get_by_id(rid)
            if r.metabolites.get(model.metabolites.get_by_id(mid)) != -1 or r.lower_bound != 0:
                raise ValueError('Input no longer encodes positive-direction hydrolysis: ' + rid)
            source = 'DIAG_SUP_' + rid
            add_boundary(model, source, mid, 1, epsilon, 'Diagnostic boundary source')
            changes['sources'][source] = {'metabolite': mid, 'cap': epsilon}
        capacity = max(1, len(supplied)) * epsilon
        def bound(rid, cap):
            r = model.reactions.get_by_id(rid)
            if r.bounds != (0, 0):
                raise ValueError('Expected closed original boundary: ' + rid)
            changes['bounds'][rid] = {'before': list(r.bounds), 'after': [0, cap]}
            r.bounds = (0, cap)
        if water:
            bound('R1363', capacity)
        if pump:
            bound('R795', capacity)
        if outputs:
            for rid in ('R871', 'R876'):
                bound(rid, epsilon)
        if drains:
            products = {}
            for rid in supplied:
                for m, c in model.reactions.get_by_id(rid).metabolites.items():
                    if c > 0:
                        products[m.id] = products.get(m.id, 0) + c * epsilon
            for i, (mid, cap) in enumerate(sorted(products.items())):
                rid = 'DIAG_DRAIN_' + str(i)
                add_boundary(model, rid, mid, -1, cap, 'Artificial diagnostic product drain')
                changes['drains'][rid] = {'metabolite': mid, 'cap': cap}
        yield changes
    assert signature(model) == before, 'Diagnostic context leaked into baseline'


def audit(model, folder):
    metabolites = set()
    target_rows, neighbor_rows = [], []
    for rid, (mid, name, _) in TARGETS.items():
        r = model.reactions.get_by_id(rid)
        target_rows.append(dict(reaction_id=rid, name=r.name, equation=r.reaction,
            named_equation=r.build_reaction_string(use_metabolite_names=True), lower_bound=r.lower_bound,
            upper_bound=r.upper_bound, reversible=r.reversibility, gpr=r.gene_reaction_rule,
            annotation=r.annotation, notes=r.notes, **balance(r)))
        # ID stem binds identical compartment copies; name also checked for other duplicates/isomers.
        normalize = lambda x: re.sub(r'[^a-z]', '', x.lower())
        pair = name.split('-')
        related = [m for m in model.metabolites if normalize(m.name).startswith(normalize(name)) or normalize(m.name).startswith(normalize('-'.join(reversed(pair))))]
        related += list(r.metabolites)
        for m in set(related):
            metabolites.add(m)
            for n in sorted(m.reactions, key=lambda r: r.id):
                neighbor_rows.append(dict(target=rid, metabolite=m.id, metabolite_name=m.name,
                    role='hydrolysis_product' if r.metabolites.get(m, 0) > 0 else ('hydrolysis_substrate' if m in r.metabolites else 'name_or_reverse_order_match'),
                    neighbor_reaction=n.id, coefficient=n.metabolites[m], equation=n.reaction,
                    bounds=list(n.bounds), gpr=n.gene_reaction_rule, **balance(n)))
                metabolites.update(n.metabolites)
    table(folder/'reaction_audit.tsv', target_rows)
    table(folder/'upstream_downstream.tsv', neighbor_rows)
    table(folder/'metabolite_audit.tsv', [dict(id=m.id, name=m.name, compartment=m.compartment,
        formula=m.formula, charge=m.charge, annotation=m.annotation, notes=m.notes) for m in sorted(metabolites, key=lambda m:m.id)])
    # Stoichiometric duplicates, irrespective of names; exact reverse columns are also recorded.
    duplicates = []
    for rid in TARGETS:
        a = model.reactions.get_by_id(rid)
        col = {m.id:c for m,c in a.metabolites.items()}
        for r in model.reactions:
            if r.id != rid and ({m.id:c for m,c in r.metabolites.items()} == col or {m.id:-c for m,c in r.metabolites.items()} == col):
                duplicates.append({'target':rid, 'other':r.id, 'gpr':r.gene_reaction_rule})
    write_json(folder/'duplicates.json', duplicates)
    trna = [r for r in model.reactions if r.id.startswith('TRNA_BIOMASS_')]
    write_json(folder/'protein_module_inventory.json', {
        'trna_net_growth_incorporation': [dict(id=r.id, equation=r.reaction, notes=r.notes) for r in trna],
        'biomass_equation':model.reactions.biomass_C.reaction, 'biomass_notes':model.reactions.biomass_C.notes,
        'keyword_matches': [dict(id=r.id,name=r.name,equation=r.reaction) for r in model.reactions if any(s in r.name.lower() for s in ('protein','peptide','autophag','turnover'))],
        'conclusion':'20 tRNA-coupled net-growth residue requirements; no connected proteolysis source to the four target pools. Keyword matches alone are not evidence of turnover.',
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reuse-run', type=Path, help='Continue only uncomputed labels from a budget-stopped batch')
    args = parser.parse_args()
    folder = args.output.resolve()
    if not folder.is_relative_to(ROOT):
        raise ValueError('Output must remain in project workspace')
    folder.mkdir(parents=True, exist_ok=True)
    if (folder/'run_manifest.json').exists():
        raise ValueError('Use a new output directory; preserve prior run')
    os.environ['IYALI26_RESEARCH_ROOT'] = str(CONDITIONS.parent)
    paths = [MODEL, MODEL.with_suffix('.build.json'), CONDITIONS/'media/sd_leu.csv', CONDITIONS/'strain_profiles/po1f_sd_leu.json']
    expected = ('aad701126d12d113816fda4b872333b614b4469ee1b6d8ab8c419231c89e965f',
                None, 'ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85',
                '35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383')
    for p, h in zip(paths, expected):
        if not p.resolve().is_relative_to(ROOT) or (h and sha(p) != h):
            raise ValueError('Input identity changed: ' + str(p))
    build = json.loads(paths[1].read_text())
    assert build['requested_build_complete'] and not build['options']['vatpase_gpr_hypothesis']
    sim = load_effective_simulation_context(model_path=paths[0], media_path=paths[2], strain_profile_path=paths[3])
    m = sim.model
    m.solver = 'gurobi'
    m.solver.configuration.presolve = False
    m.solver.problem.Params.Threads = 1
    m.solver.problem.Params.TimeLimit = 60
    assert {r.id:c for r,c in linear_reaction_coefficients(m).items()} == {'biomass_C':1.0}
    assert all(m.reactions.get_by_id(r).lower_bound == 0 for r in ('R2018','R2026','R2031','R2036'))
    base = signature(m)
    audit(m, folder)
    write_json(folder/'effective_model.json', base)
    spec = {'enabled':False, 'status':'specification_only', **dict.fromkeys(ENDO_REQUIRED)}
    write_json(folder/'endogenous_module_spec.json', spec)
    write_json(folder/'endogenous_validation.json', validate_endogenous_spec(spec))
    started = time.monotonic()
    record = dict(status='running', started_utc=datetime.now(timezone.utc).isoformat(),
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths}, script_sha256=sha(__file__),
        source_sha256={str(Path(x.__file__).relative_to(ROOT)):sha(x.__file__) for x in list(sys.modules.values()) if getattr(x,'__file__',None) and Path(x.__file__).is_relative_to(ROOT/'scripts') and Path(x.__file__).is_file()},
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        software={'python':sys.version,**{p:metadata.version(p) for p in ('cobra','optlang','gurobipy','memote')}},
        solver={k:getattr(m.solver.problem.Params,k) for k in ('Threads','TimeLimit','Presolve','FeasibilityTol','OptimalityTol','IntFeasTol')},
        flux_tolerance=TOL, growth_slack=GROWTH_SLACK, supply_levels=LEVELS,
        simulation_context=sim.provenance(), active_medium=sim.active_medium, strain_overlay=sim.strain_overlay_audit,
        baseline_objective={'biomass_C':1}, primary_solves=0, scenarios=[], outcomes=[],
        limits={'solves':900,'wall_seconds':1200,'per_solve_seconds':60},
        native_parameters_measured=False)
    fluxfile = gzip.open(folder/'flux_witnesses.jsonl.gz', 'wt')
    rows, sensitivity, pathways = [], [], []
    prior, witnesses, completed_cases = {}, {}, set()
    if args.reuse_run:
        previous = args.reuse_run.resolve()
        if not previous.is_relative_to(ROOT):
            raise ValueError('Prior run must remain inside workspace')
        old=json.loads((previous/'run_manifest.json').read_text())
        if old['status']!='incomplete' or old.get('error')!="RuntimeError('Declared compute budget reached')":
            raise ValueError('Only a budget stop can be continued; do not retry scientific/solver failures')
        driver=str(Path(__file__).relative_to(ROOT))
        loader_sources=lambda x:{k:v for k,v in x.items() if k!=driver}
        archived=ROOT/'artifacts/dipeptide_supply_repair_20260923/diagnose_v1.py'
        if old['script_sha256'] not in (sha(__file__),sha(archived)):
            raise ValueError('Unreviewed prior diagnostic driver')
        if old['input_sha256']!=record['input_sha256'] or old['solver']!=record['solver'] or loader_sources(old['source_sha256'])!=loader_sources(record['source_sha256']):
            raise ValueError('Inputs, loader or solver changed; cannot reuse previous outcomes')
        if old['simulation_context']!=record['simulation_context'] or old['supply_levels']!=list(LEVELS):
            raise ValueError('Conditions or supply levels changed')
        completed_cases={c['name'] for c in old['scenarios'] if c.get('context_completed')}
        record['scenarios']=[c for c in old['scenarios'] if c['name'] in completed_cases]
        record['outcomes']=old['outcomes']
        prior={r['label']:r for r in old['outcomes']}
        assert len(prior)==len(old['outcomes'])
        for k in ('baseline_WT','baseline_exchange_fluxes','baseline_glucose_C_uptake'):
            record[k]=old[k]
        record['reused_outcomes']=len(prior)
        record['resume_provenance']={'run':str(previous.relative_to(ROOT)),
            'manifest_sha256':sha(previous/'run_manifest.json'),
            'fluxes_sha256':sha(previous/'flux_witnesses.jsonl.gz'),
            'previous_script_sha256':old['script_sha256'],
            'completed_label_sha256':hashlib.sha256(json.dumps(sorted(prior)).encode()).hexdigest(),
            'reason':'Continue uncomputed labels after declared wall-time stop; no completed LP repeated'}
        record['limits']['wall_seconds']=600
        with gzip.open(previous/'flux_witnesses.jsonl.gz','rt') as f:
            for line in f:
                w=json.loads(line);witnesses[w['label']]=w;fluxfile.write(line)
        with (previous/'scenario_summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        with (previous/'supply_sensitivity.tsv').open() as f:
            sensitivity=[r for r in csv.DictReader(f,delimiter='\t') if r['scenario'] in completed_cases]
        with (previous/'product_flux_routes.tsv').open() as f:
            pathways=[r for r in csv.DictReader(f,delimiter='\t') if any(r['label'].startswith(c+'/') for c in completed_cases)]

    def save():
        record['elapsed_seconds'] = time.monotonic()-started
        write_json(folder/'run_manifest.json', record)

    def solve(label, model=m, extra=None):
        if label in prior:
            return cached_solution(prior[label],witnesses.get(label),extra)
        if record['primary_solves']+len(prior) >= 900 or time.monotonic()-started > record['limits']['wall_seconds']:
            raise RuntimeError('Declared compute budget reached')
        record['primary_solves'] += 1
        sol = model.optimize()
        status = str(sol.status)
        row = {'label':label,'status':status,'objective':float(sol.objective_value) if sol.objective_value is not None and math.isfinite(sol.objective_value) else None,
               'objective_expression':str(model.objective.expression) if extra and extra.get('kind')=='energy' else None,
               **(extra or {})}
        if status == 'optimal':
            flux = {r.id:float(sol.fluxes[r.id]) for r in model.reactions}
            finite = all(math.isfinite(v) for v in flux.values())
            if not finite:
                raise RuntimeError('Nonfinite flux result')
            mass = max(abs(sum(r.metabolites[x]*flux[r.id] for r in x.reactions)) for x in model.metabolites)
            bound = max(max(r.lower_bound-flux[r.id],flux[r.id]-r.upper_bound,0) for r in model.reactions)
            row.update(growth=flux.get('biomass_C'),hydrolysis={r:flux.get(r) for r in TARGETS},
                       max_mass_residual=mass,max_bound_violation=bound)
            row['max_constraint_violation'] = max([0.0]+[max((c.lb-c.primal) if c.lb is not None else 0,(c.primal-c.ub) if c.ub is not None else 0,0) for c in model.constraints])
            row['numerically_valid'] = max(mass,bound,row['max_constraint_violation']) <= TOL
            fluxfile.write(json.dumps({'label':label,'fluxes':flux,'bounds_changes':{r.id:list(r.bounds) for r in model.reactions if r.id not in base['reactions'] or list(r.bounds)!=base['reactions'][r.id][:2]},'extra':extra},allow_nan=False)+'\n')
            fluxfile.flush()
        else:
            row.update(growth=None,hydrolysis=None,numerically_valid=False)
        record['outcomes'].append(row)
        rows.append({k:row.get(k) for k in ('label','status','objective','growth','hydrolysis','numerically_valid','max_mass_residual','max_bound_violation','max_constraint_violation')})
        save()
        if status not in ('optimal','infeasible') or (status=='optimal' and not row['numerically_valid']):
            raise RuntimeError('Unresolved solver outcome retained; no automatic retry: '+label)
        return row, sol

    def floor(value):
        m.add_cons_vars(m.problem.Constraint(m.reactions.biomass_C.flux_expression, lb=value, name='DIAG_GROWTH_FLOOR'))

    def target_solve(label, rid, direction, growth_floor=0.0):
        with m:
            floor(growth_floor)
            m.objective = m.reactions.get_by_id(rid)
            m.objective.direction = direction
            row, sol = solve(label, extra={'target':rid,'direction':direction,'growth_floor':growth_floor})
            if row['status']=='optimal' and direction=='max' and row['hydrolysis'][rid]>TOL:
                for metabolite,c in m.reactions.get_by_id(rid).metabolites.items():
                    if c>0:
                        for r in sorted(metabolite.reactions,key=lambda r:r.id):
                            v=float(sol.fluxes[r.id]);rate=v*r.metabolites[metabolite]
                            if abs(rate)>TOL:
                                pathways.append(dict(label=label,target=rid,product=metabolite.id,reaction=r.id,
                                    flux=v,product_rate=rate,equation=r.reaction,role='consumption' if rate<0 else 'production'))
            return row

    def run_case(name, supplied=(), epsilon=0, **options):
        if name in completed_cases:
            return
        with scenario(m,supplied,epsilon,**options) as changes:
            case = dict(name=name,supplied=list(supplied),epsilon=epsilon,options=options,changes=changes)
            record['scenarios'].append(case)
            fba, sol = solve(name+'/FBA',extra={'scenario':name,'kind':'growth'})
            if fba['status']!='optimal':
                case['unresolved']='Growth optimization infeasible; no essentiality classification'
                return
            wt=fba['growth']
            if name=='baseline':
                record['baseline_WT']=wt
                record['baseline_exchange_fluxes']={r.id:float(sol.fluxes[r.id]) for r in m.exchanges}
                record['baseline_glucose_C_uptake']=max(0,-float(sol.fluxes['R1070']))*6
            with m:
                floor(max(0,wt-GROWTH_SLACK))
                m.objective=m.problem.Objective(sum(r.forward_variable+r.reverse_variable for r in m.reactions),direction='min')
                pfba,_=solve(name+'/pFBA',extra={'scenario':name,'kind':'pFBA','growth_floor':max(0,wt-GROWTH_SLACK)})
            targets=list(supplied) or list(TARGETS)
            for rid in targets:
                maximum=target_solve(name+'/'+rid+'/free_max',rid,'max')
                near=[]
                for direction in ('min','max'):
                    near.append(target_solve(name+'/'+rid+'/near_'+direction,rid,direction,max(0,0.99*wt-GROWTH_SLACK)))
                sensitivity.append(dict(scenario=name,target=rid,epsilon=epsilon,
                    growth=wt,delta_growth=wt-record['baseline_WT'],free_status=maximum['status'],
                    free_max=maximum['hydrolysis'][rid] if maximum['hydrolysis'] else None,
                    near_min=near[0]['hydrolysis'][rid] if near[0]['hydrolysis'] else None,
                    near_max=near[1]['hydrolysis'][rid] if near[1]['hydrolysis'] else None,
                    fba_flux=fba['hydrolysis'][rid],pfba_flux=pfba['hydrolysis'][rid] if pfba['hydrolysis'] else None,
                    supplied_C_model=None if supplied else 0,supplied_N_model=None if supplied else 0,
                    nominal_total_C_cap=sum(TARGETS[r][2]*epsilon for r in supplied),
                    nominal_total_N_cap=2*len(supplied)*epsilon,
                    nominal_C_fraction_of_WT_glucose_C=sum(TARGETS[r][2]*epsilon for r in supplied)/record['baseline_glucose_C_uptake'],
                    input_accounting='Model dipeptide formula absent: complete balance unverified; nominal C/N assumes named standard dipeptide identity'))
                if name=='baseline':
                    for fraction in (1.0,0.0):
                        for direction in ('min','max'):
                            target_solve(name+'/'+rid+'/FVA_'+str(fraction)+'_'+direction,rid,direction,max(0,fraction*wt-GROWTH_SLACK))
            if len(supplied)==4:
                with m:
                    for rid in supplied:
                        m.add_cons_vars(m.problem.Constraint(m.reactions.get_by_id(rid).flux_expression,lb=epsilon/2,name='DIAG_JOINT_'+rid))
                    floor(max(0,0.99*wt-GROWTH_SLACK))
                    case['joint_test']=solve(name+'/joint_half_supply',extra={'scenario':name,'kind':'joint','hydrolysis_min':epsilon/2,'growth_floor':max(0,0.99*wt-GROWTH_SLACK)})[0]
            case['context_completed']=True
        print(name,'growth',wt,'done',flush=True)

    try:
        run_case('baseline')
        variants={
            'supply':{},
            'supply_drains':{'drains':True},
            'supply_water':{'water':True},
            'supply_water_drains':{'water':True,'drains':True},
            'supply_water_pump':{'water':True,'pump':True},
            'supply_water_pump_outputs':{'water':True,'pump':True,'outputs':True},
        }
        for epsilon in LEVELS:
            for label,options in variants.items():
                for supplied in [(r,) for r in TARGETS]+[tuple(TARGETS)]:
                    run_case(f'{label}/{"all" if len(supplied)==4 else supplied[0]}/{epsilon:g}',supplied,epsilon,**options)
        # No dipeptide input controls distinguish repair connections from an input source.
        run_case('water_pump_outputs_no_source',(),0.01,water=True,pump=True,outputs=True)
        # Reuse the project's memote closed-boundary ATP protocol, with explicit ATPM control.
        record['energy_protocol_sha256']=sha(ROOT/'artifacts/trna_biomass_restore_20260910/static_check.py')
        record['memote_energy_function_sha256']=hashlib.sha256(inspect.getsource(consistency.detect_energy_generating_cycles).encode()).hexdigest()
        record['energy_checks']=[]
        for variant in ('baseline','hypothesis_connections'):
            for relax in (False,True):
                with scenario(m,(),0.01,water=variant!='baseline',pump=variant!='baseline',outputs=variant!='baseline'):
                    with m:
                        if relax:
                            m.reactions.xMAINTENANCE.lower_bound=0
                        label=f'energy/{variant}/ATPM_{"zero" if relax else "preserved"}'
                        item={'label':label,'ATPM_bounds':list(m.reactions.xMAINTENANCE.bounds),'closed_boundaries':[r.id for r in m.boundary]}
                        def close_boundaries(model):
                            for r in model.boundary:r.bounds=(0,0)
                            return model
                        raw=m.optimize
                        def capture(*a,**kw):
                            m.optimize=raw
                            try:
                                result,sol=solve(label,extra={'kind':'energy','ATPM_relaxed':relax})
                                item.update(result)
                                item['dissipation_equation']=m.reactions.Dissipation.reaction
                                item['dissipation_balance']=balance(m.reactions.Dissipation)
                                return sol
                            finally:m.optimize=capture
                        with patch.object(helpers,'close_boundaries_sensibly',close_boundaries),patch.object(m,'optimize',capture):
                            item['cycle_reactions']=consistency.detect_energy_generating_cycles(m,'MNXM3')
                        record['energy_checks'].append(item)
        # Stored neutral species may differ from memote's proton convention.
        # Also maximize the existing maintenance chemistry, without adding a column.
        record['closed_local_cycles']=[]
        for variant in ('baseline','hypothesis_connections'):
            with scenario(m,(),0.01,water=variant!='baseline',pump=variant!='baseline',outputs=variant!='baseline'):
                with m:
                    for r in m.boundary:r.bounds=(0,0)
                    m.reactions.xMAINTENANCE.lower_bound=0
                    for rid in ['xMAINTENANCE',*TARGETS,'R795','R1363','R871','R876','R2030','R2035','R2040']:
                        for direction in (('max',) if rid=='xMAINTENANCE' else ('min','max')):
                            m.objective=m.reactions.get_by_id(rid)
                            m.objective.direction=direction
                            label=f'closed_local/{variant}/{rid}/{direction}'
                            row,_=solve(label,extra={'kind':'closed_local_cycle','target':rid,'direction':direction,
                                'ATPM_bounds':list(m.reactions.xMAINTENANCE.bounds),'all_boundaries_closed':True})
                            record['closed_local_cycles'].append(row)
        record['baseline_restored']=signature(m)==base
        assert record['baseline_restored']
        record['inputs_unchanged']=all(sha(ROOT/p)==h for p,h in record['input_sha256'].items())
        assert record['inputs_unchanged']
        record.update(status='complete',completed_utc=datetime.now(timezone.utc).isoformat())
    except BaseException as error:
        record.update(status='incomplete',error=repr(error))
        raise
    finally:
        fluxfile.close()
        table(folder/'scenario_summary.tsv',rows)
        table(folder/'supply_sensitivity.tsv',sensitivity)
        table(folder/'product_flux_routes.tsv',pathways)
        save()


if __name__=='__main__':
    main()
