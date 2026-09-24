"""Bounded validation of exported energy candidates; no scientific patches here."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib import metadata
import json
import math
import os
from pathlib import Path
import sys
import time

from cobra import Reaction

from scripts.diagnose_closed_energy import close_model, configure_solver, block_direction, workspace
from scripts.diagnose_dipeptide_supply import balance, sha, signature, table
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

ROOT = Path(__file__).resolve().parents[1]
TASK = ROOT / 'artifacts/atp_candidate_repair_20260924'


def json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return 'Infinity' if value > 0 else '-Infinity' if value < 0 else 'NaN'
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def write(path, value):
    Path(path).write_text(json.dumps(json_safe(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def energy_verdict(status, value, tolerance):
    if status != 'optimal' or value is None or not math.isfinite(value) or value < -tolerance:
        return 'unresolved'
    return 'within_tolerance' if abs(value) <= tolerance else 'positive_energy_regeneration'


class SolverBudget:
    """One persistent serial ledger for diagnostics and behavior tests."""
    def __init__(self, path, config):
        self.path, self.config = Path(path), config
        if self.path.exists():
            self.record = json.loads(self.path.read_text())
            assert self.record['limits'] == config['limits']
            assert self.record['solver'] == config['solver']
            if any(r['status'] == 'running' for r in self.record['calls']):
                raise RuntimeError('Unresolved interrupted solver call; do not erase its budget')
        else:
            self.record = {'limits': config['limits'], 'solver': config['solver'],
                           'calls': [], 'solve_call_seconds': 0.0}
            write(self.path, self.record)

    def solve(self, model, label, out, purpose):
        if len(self.record['calls']) >= self.config['limits']['solves'] or (
            self.config['limits']['solve_wall_seconds']-self.record['solve_call_seconds'] <
            self.config['limits']['per_solve_seconds']):
            raise RuntimeError('Declared aggregate solver budget reached')
        if any(r['label'] == label for r in self.record['calls']):
            raise ValueError('Solver label already used: '+label)
        params = configure_solver(model, self.config)
        model.solver.update()
        problem = signature(model)
        row = {'label': label, 'status': 'running', 'objective': None,
               'purpose': purpose, 'solver': params,
               'started_utc': datetime.now(timezone.utc).isoformat()}
        self.record['calls'].append(row)
        write(self.path, self.record)
        solution = None
        start = time.monotonic()
        try:
            solution = model.optimize()
            row['status'] = str(solution.status)
        except Exception as exc:
            row.update(status='software_error', error=repr(exc))
            raise
        finally:
            row['solve_call_seconds'] = time.monotonic()-start
            self.record['solve_call_seconds'] += row['solve_call_seconds']
            write(self.path, self.record)
            write(out, {'result': row, 'problem': problem, 'fluxes': None})
        flux = None
        try:
            value = solution.objective_value
            row['objective'] = float(value) if value is not None and math.isfinite(value) else None
            if row['status'] == 'optimal':
                flux = {r.id: float(solution.fluxes[r.id]) for r in model.reactions}
                if row['objective'] is None or not all(math.isfinite(v) for v in flux.values()):
                    raise RuntimeError('Nonfinite optimal result')
                mass = max(abs(sum(r.metabolites[m]*flux[r.id] for r in m.reactions))
                           for m in model.metabolites)
                bounds = max(max(r.lower_bound-flux[r.id], flux[r.id]-r.upper_bound, 0)
                             for r in model.reactions)
                constraints = max([0]+[max(c.lb-c.primal if c.lb is not None else 0,
                                           c.primal-c.ub if c.ub is not None else 0, 0)
                                         for c in model.constraints])
                row.update(max_mass_residual=mass, max_bound_violation=bounds,
                           max_constraint_violation=constraints)
                if max(mass, bounds, constraints) > self.config['tolerance']:
                    raise RuntimeError('Numerical validation failed: '+label)
            row['postprocessing'] = 'complete'
        except Exception as exc:
            row.update(postprocessing='failed', postprocessing_error=repr(exc))
            raise
        finally:
            write(out, {'result': row, 'problem': problem, 'fluxes': flux})
            write(self.path, self.record)
        print(label, row['status'], row['objective'], flush=True)
        if row['status'] not in ('optimal', 'infeasible'):
            raise RuntimeError('Unresolved solver status: '+row['status'])
        return row, flux


def load_model(path, config):
    sim = load_effective_simulation_context(model_path=workspace(path),
        media_path=workspace(config['media']), strain_profile_path=workspace(config['strain_profile']))
    sim.model.solver = 'gurobi'
    configure_solver(sim.model, config)
    assert sim.model.reactions.get_by_id(config['maintenance']).lower_bound == 7.8625
    assert sim.model.objective.direction == 'max'
    return sim


def reaction_record(r):
    return {'reaction': r.id, 'name': r.name, 'equation': r.reaction,
            'named_equation': r.build_reaction_string(use_metabolite_names=True),
            'stoichiometry': {m.id: c for m, c in r.metabolites.items()},
            'bounds': list(r.bounds), 'compartments': sorted(r.compartments),
            'gpr': r.gene_reaction_rule, 'gpr_evidence': 'stored model assignment; not native validation',
            'annotation': r.annotation, 'notes': r.notes, **balance(r)}


def export_witness(model, flux, label, config, output):
    dissipation = config['maintenance']
    rows = [{**reaction_record(r), 'witness': label, 'flux': flux[r.id],
             'actual_direction': 'forward' if flux[r.id] > 0 else 'reverse'}
            for r in model.reactions if abs(flux[r.id]) > config['display_threshold']]
    ledger = []
    # Match species by stable ID; no model-copy Python object identity assumptions.
    columns = {r.id: {m.id: c for m, c in r.metabolites.items()} for r in model.reactions}
    for met in model.metabolites:
        contributions = {rid: stoich.get(met.id, 0)*flux[rid] for rid, stoich in columns.items()
                         if rid != dissipation and stoich.get(met.id, 0)*flux[rid] != 0}
        d = columns[dissipation].get(met.id, 0)*flux[dissipation]
        if contributions or d:
            ledger.append({'witness': label, 'metabolite': met.id, 'name': met.name,
                'compartment': met.compartment, 'formula': met.formula, 'charge': met.charge,
                'contributions_excluding_D': contributions, 'net_excluding_D': sum(contributions.values()),
                'D_contribution': d, 'total_residual': sum(contributions.values())+d})
    assert max(abs(r['total_residual']) for r in ledger) <= config['tolerance']
    table(output/'atp_witness_fluxes.tsv', rows)
    table(output/'atp_net_reaction_ledger.tsv', ledger)
    write(output/'witness.json', {'label': label, 'reactions': rows, 'ledger': ledger,
         'minimization': 'L1 flux, not minimum reaction count', 'D': flux[dissipation]})


def add_ntp_dissipation(model, ntp):
    # Explicit IDs and formulas pinned to the audited current neutral nucleotide convention.
    species = {
        'GTP': ('m266[C_cy]', 'm268[C_cy]', 'C10H16N5O14P3', 'C10H15N5O11P2'),
        'UTP': ('m439[C_cy]', 'm11[C_cy]', 'C9H15N2O15P3', 'C9H14N2O12P2'),
        'CTP': ('m406[C_cy]', 'm500[C_cy]', 'C9H16N3O14P3', 'C9H15N3O11P2'),
    }
    trip, dip, tf, df = species[ntp]
    for mid, formula in [(trip, tf), (dip, df), ('m32[C_cy]', 'H2O'), ('m35[C_cy]', 'H3O4P')]:
        m = model.metabolites.get_by_id(mid)
        if m.formula != formula or m.charge != 0:
            raise ValueError('Unverified nucleotide test species: '+repr((mid,m.formula,m.charge)))
    r = Reaction('DIAG_D_'+ntp, name='Diagnostic balanced '+ntp+' dissipation', lower_bound=0, upper_bound=1000)
    r.add_metabolites({model.metabolites.get_by_id(trip):-1, model.metabolites.get_by_id('m32[C_cy]'):-1,
                      model.metabolites.get_by_id(dip):1, model.metabolites.get_by_id('m35[C_cy]'):1})
    if balance(r)['element_status'] != 'balanced' or balance(r)['charge_status'] != 'balanced_as_stored':
        raise ValueError('Dissipation is not balanced')
    model.add_reactions([r])
    model.objective = r
    model.objective.direction = 'max'
    return r


def run(args):
    config = json.loads(workspace(args.config).read_text())
    assert sha(workspace(config['model'])) == config['model_sha256']
    out = workspace(args.output); out.mkdir(parents=True, exist_ok=False)
    budget = SolverBudget(workspace(args.budget), config)
    started = time.monotonic()
    manifest = {'status':'running', 'config':config, 'config_sha256':sha(workspace(args.config)),
        'script_sha256':sha(__file__), 'started_utc':datetime.now(timezone.utc).isoformat(),
        'software': {'python':sys.version, **{p:metadata.version(p) for p in ('cobra','optlang','gurobipy','memote')}},
        'source_sha256':{str(Path(m.__file__).relative_to(ROOT)):sha(m.__file__) for m in list(sys.modules.values())
                        if getattr(m,'__file__',None) and Path(m.__file__).is_relative_to(ROOT/'scripts')},
        'variants':{}, 'outcomes':[]}
    energies, growth = [], []

    def save():
        manifest['wall_seconds'] = time.monotonic()-started
        write(out/'manifest.json', manifest)
        table(out/'energy_iterations.tsv', energies)
        table(out/'growth_comparison.tsv', growth)

    def solve(model, label, purpose):
        row, flux = budget.solve(model, out.name+'/'+label, out/(label+'.json'), purpose)
        manifest['outcomes'].append(row); save()
        return row, flux

    def maximum(closed, label, variant, carrier='ATP', blocks=None):
        row, flux = solve(closed, label, {'kind':'closed_max', 'variant':variant,'carrier':carrier,'blocks':blocks or []})
        energies.append({'variant':variant,'label':label,'carrier':carrier,'blocks':blocks or [],
                         'status':row['status'],'maximum':row['objective'],
                         'verdict':energy_verdict(row['status'],row['objective'],config['tolerance'])})
        save(); return row, flux

    try:
        for entry in args.model:
            name, path = entry.split('=',1)
            sim = load_model(path, config); original = sim.model; before = signature(original)
            closed, changes = close_model(original, config)
            maintenance = closed.reactions.get_by_id(config['maintenance'])
            assert balance(maintenance)['element_status']=='balanced'
            assert balance(maintenance)['charge_status']=='balanced_as_stored'
            manifest['variants'][name] = {'path':path,'sha256':sha(workspace(path)),
                'simulation_context':sim.provenance(),'active_medium':sim.active_medium,
                'strain_overlay':sim.strain_overlay_audit,'closure':changes,
                'maintenance':reaction_record(maintenance)}
            write(out/(name+'_effective.json'), before)
            write(out/(name+'_closed.json'), signature(closed))
            with closed:
                for r in closed.reactions: r.bounds=(0,0)
                closed.objective=closed.problem.Objective(0)
                zero, zflux=solve(closed,name+'_zero',{'kind':'all_zero_feasibility','variant':name})
                if zero['status']!='optimal': raise RuntimeError('Closed zero flux infeasible')
            atp,_=maximum(closed,name+'_ATP',name)
            if name in args.witness and atp['status']=='optimal' and atp['objective']>config['tolerance']:
                d = min(1.0, atp['objective']/2) if atp['objective']<1 else 1.0
                if d<=config['tolerance']: raise RuntimeError('No normalization safely above tolerance')
                with closed:
                    maintenance.bounds=(d,d)
                    closed.objective=closed.problem.Objective(sum(r.forward_variable+r.reverse_variable
                                                for r in closed.reactions if r.id!=maintenance.id),direction='min')
                    row,flux=solve(closed,name+'_witness',{'kind':'fixed_D_L1','variant':name,'D':d})
                    if row['status']!='optimal': raise RuntimeError('Positive normalization failed')
                    folder=out/name;folder.mkdir()
                    export_witness(closed,flux,name,config,folder)
            if name in args.other_energy:
                for ntp in ('GTP','UTP','CTP'):
                    with closed:
                        test=add_ntp_dissipation(closed,ntp)
                        manifest['variants'][name].setdefault('energy_test_chemistry',{})[ntp]=reaction_record(test)
                        maximum(closed,name+'_'+ntp,name,ntp)
            if args.blocks:
                for case in json.loads(workspace(args.blocks).read_text()):
                    if case['variant']!=name: continue
                    with closed:
                        for rid,direction in case['blocks']: block_direction(closed,rid,direction)
                        maximum(closed,name+'_'+case['label'],name,blocks=case['blocks'])
            if not args.no_growth:
                row,flux=solve(original,name+'_growth',{'kind':'real_culture_growth','variant':name})
                growth.append({'variant':name,'status':row['status'],'growth':row['objective'],
                    'maintenance':flux[config['maintenance']] if flux else None,
                    'maintenance_bounds':list(original.reactions.get_by_id(config['maintenance']).bounds),
                    'objective':str(original.objective.expression),
                    'actual_boundary_fluxes':{r.id:flux[r.id] for r in original.boundary if abs(flux[r.id])>config['tolerance']} if flux else None,
                    'ATP_reaction_fluxes':{r.id:flux[r.id] for r in original.metabolites.get_by_id('m141[C_cy]').reactions if abs(flux[r.id])>config['tolerance']} if flux else None,
                    'R1372_flux':flux.get('R1372') if flux else None})
            assert signature(original)==before, 'Temporary diagnostic modifications leaked'
            manifest['variants'][name]['original_restored']=True
            save()
        manifest['status']='complete'
    except Exception as exc:
        manifest.update(status='incomplete',error=repr(exc))
        raise
    finally:
        save()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default=str(TASK/'config.json'))
    parser.add_argument('--budget',default=os.environ.get('IYALI26_ENERGY_BUDGET',str(TASK/'budget.json')))
    parser.add_argument('--output',required=True)
    parser.add_argument('--model',action='append',required=True,help='Unique version name=workspace XML path')
    parser.add_argument('--witness',action='append',default=[])
    parser.add_argument('--other-energy',action='append',default=[])
    parser.add_argument('--blocks')
    parser.add_argument('--no-growth',action='store_true')
    run(parser.parse_args())


if __name__=='__main__':
    main()
