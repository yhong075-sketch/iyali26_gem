"""Bounded closed-energy witnesses; all scientific edits stay in model copies."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from importlib import metadata
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from scripts.diagnose_dipeptide_supply import balance, sha, signature, table, write_json
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'artifacts/dipeptide_energy_audit_20260924'


def configure_solver(model, config):
    # COBRA copies reset native Gurobi parameters: set and verify on every solve.
    for key, value in config['solver'].items():
        setattr(model.solver.problem.Params, key, value)
    actual = {k:getattr(model.solver.problem.Params,k) for k in config['solver']}
    assert actual == config['solver'], actual
    return actual


def workspace(path):
    path = (ROOT / path).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError('Path escapes authorized workspace: ' + str(path))
    return path


def load_baseline(config):
    paths = {k: workspace(config[k]) for k in ('model', 'media', 'strain_profile')}
    assert sha(paths['model']) == config['model_sha256'], 'Model identity changed'
    sim = load_effective_simulation_context(model_path=paths['model'], media_path=paths['media'],
                                             strain_profile_path=paths['strain_profile'])
    sim.model.solver = 'gurobi'
    configure_solver(sim.model, config)
    return sim


def close_model(original, config):
    model = original.copy()
    configure_solver(model, config)
    changes = []
    biomass = set(config['biomass_reactions'])
    assert biomass <= set(r.id for r in model.reactions)
    for r in model.reactions:
        reasons = []
        if not r.reactants or not r.products:
            reasons.append('single_sided_source_sink_exchange')
        if r.id in biomass:
            reasons.append('biomass_synthesis')
        if reasons:
            changes.append({'reaction': r.id, 'reasons': reasons, 'before': list(r.bounds), 'after': [0, 0]})
            r.bounds = (0, 0)
        elif r.id == config['maintenance']:
            changes.append({'reaction': r.id, 'reasons': ['release_forced_maintenance'],
                            'before': list(r.bounds), 'after': [0, r.upper_bound]})
            r.lower_bound = 0
        elif r.lower_bound > 0 or r.upper_bound < 0:
            raise ValueError('Unreviewed forced internal flux: ' + r.id)
    model.objective = model.reactions.get_by_id(config['maintenance'])
    model.objective.direction = 'max'
    model.solver.update()
    extra = [c.name for c in model.constraints if c.name not in model.metabolites]
    if extra:
        raise ValueError('Unreviewed custom constraints: ' + repr(extra))
    assert len(model.constraints) == len(model.metabolites)
    zero_violations = [v.name for v in model.variables if
                       (v.lb is not None and v.lb > 0) or (v.ub is not None and v.ub < 0)]
    zero_violations += [c.name for c in model.constraints if
                       (c.lb is not None and c.lb > 0) or (c.ub is not None and c.ub < 0)]
    assert not zero_violations, zero_violations
    return model, changes


def block_direction(model, rid, direction):
    r = model.reactions.get_by_id(rid)
    if direction == 'forward':
        r.upper_bound = min(r.upper_bound, 0)
    elif direction == 'reverse':
        r.lower_bound = max(r.lower_bound, 0)
    else:
        raise ValueError(direction)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=DEFAULT / 'energy')
    parser.add_argument('--config', type=Path, default=DEFAULT / 'config.json')
    parser.add_argument('--reuse-run', type=Path)
    parser.add_argument('--extra-cases', type=Path)
    args = parser.parse_args()
    out, config_path = workspace(args.output), workspace(args.config)
    out.mkdir(parents=True, exist_ok=False)
    config = json.loads(config_path.read_text())
    sim = load_baseline(config)
    original = sim.model
    protected = signature(original)
    previous = json.loads((ROOT/'artifacts/dipeptide_supply_repair_20260923/run_complete/run_manifest.json').read_text())
    assert sim.provenance() == previous['simulation_context'], 'Historical/current context differs'
    closed, changes = close_model(original, config)
    dissipation = closed.reactions.get_by_id(config['maintenance'])
    chemistry = balance(dissipation)
    assert chemistry['element_status'] == 'balanced' and chemistry['charge_status'] == 'balanced_as_stored'
    manifest = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256': {config[k]: sha(workspace(config[k])) for k in ('model','media','strain_profile')},
        'config': config, 'config_sha256': sha(config_path), 'script_sha256': sha(__file__),
        'source_sha256': {str(Path(m.__file__).relative_to(ROOT)): sha(m.__file__) for m in list(sys.modules.values())
                          if getattr(m,'__file__',None) and Path(m.__file__).is_relative_to(ROOT/'scripts')},
        'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'software': {'python':sys.version, **{p:metadata.version(p) for p in ('cobra','optlang','gurobipy','memote')}},
        'simulation_context': sim.provenance(), 'active_medium': sim.active_medium,
        'strain_overlay': sim.strain_overlay_audit, 'original_objective': str(original.objective.expression),
        'solver': {k:getattr(closed.solver.problem.Params,k) for k in config['solver']},
        'closure': changes, 'zero_flux_feasible_by_all_variable_and_constraint_bounds': True,
        'custom_constraints': [], 'maintenance': {'id':dissipation.id,'equation':dissipation.reaction,
                         'named_equation':dissipation.build_reaction_string(use_metabolite_names=True), **chemistry},
        'solves': [], 'solve_wall_seconds': 0.0, 'witnesses': [],
        'prior_attempt_budget': config.get('prior_attempt_budget', {'solves':0,'reserved_seconds':0})}
    write_json(out/'closed_model.json', signature(closed))
    interventions, flux_rows, ledger = [], [], []
    reused = {}
    if args.reuse_run:
        prior_path = workspace(args.reuse_run)
        prior = json.loads((prior_path/'closed_energy_manifest.json').read_text())
        assert prior['config'] == config and prior['input_sha256'] == manifest['input_sha256']
        assert prior['simulation_context'] == sim.provenance()
        assert prior['script_sha256'] in [sha(DEFAULT/p) for p in
            ('diagnose_closed_energy_attempt2.py','diagnose_closed_energy_initial_complete.py','diagnose_closed_energy_extended.py')]
        assert json.loads((prior_path/'closed_model.json').read_text()) == signature(closed)
        reused = {r['label']:json.loads((prior_path/(r['label']+'.json')).read_text())
                  for r in prior['solves'] if r['status']=='optimal'}
        manifest['solve_wall_seconds'] = prior['solve_wall_seconds']
        manifest['reuse'] = {'path':str(prior_path),'manifest_sha256':sha(prior_path/'closed_energy_manifest.json'),
                            'reason':'Ledger-only object-identity fix; saved optimization problems unchanged',
                            'labels':list(reused)}

    def save():
        write_json(out/'closed_energy_manifest.json', manifest)
        table(out/'direction_interventions.tsv', interventions)
        table(out/'atp_witness_fluxes.tsv', flux_rows)
        table(out/'atp_net_reaction_ledger.tsv', ledger)

    def solve(model, label, purpose):
        if label in reused:
            saved = reused[label]
            row, flux = saved['result'], saved['fluxes']
            assert row['purpose']==purpose and row['solver']==config['solver']
            assert saved['bounds']=={r.id:list(r.bounds) for r in model.reactions}
            value = sum(abs(v) for r,v in flux.items() if r!=config['maintenance']) if purpose['kind']=='fixed_dissipation_L1' else flux[config['maintenance']]
            assert abs(value-row['objective']) <= config['tolerance']
            assert max(abs(sum(r.metabolites[m]*flux[r.id] for r in m.reactions)) for m in model.metabolites) <= config['tolerance']
            manifest['solves'].append(row)
            write_json(out/(label+'.json'),saved)
            save()
            print(label,'reused',row['objective'],flush=True)
            return row,flux
        prior = manifest['prior_attempt_budget']
        if len(manifest['solves'])+prior['solves'] >= config['limits']['solves'] or manifest['solve_wall_seconds']+prior['reserved_seconds'] >= config['limits']['solve_wall_seconds']:
            raise RuntimeError('Declared compute budget reached')
        solver_parameters = configure_solver(model, config)
        start = time.monotonic()
        solution = model.optimize()
        manifest['solve_wall_seconds'] += time.monotonic()-start
        row = {'label':label, 'purpose':purpose, 'status':str(solution.status),
               'objective':float(solution.objective_value) if solution.objective_value is not None and math.isfinite(solution.objective_value) else None,
               'solver':solver_parameters}
        manifest['solves'].append(row)
        if solution.status == 'optimal':
            flux = {r.id:float(solution.fluxes[r.id]) for r in model.reactions}
            assert all(math.isfinite(v) for v in flux.values())
            mass = max(abs(sum(r.metabolites[m]*flux[r.id] for r in m.reactions)) for m in model.metabolites)
            bound = max(max(r.lower_bound-flux[r.id], flux[r.id]-r.upper_bound, 0) for r in model.reactions)
            cons = max([0]+[max(c.lb-c.primal if c.lb is not None else 0,c.primal-c.ub if c.ub is not None else 0,0) for c in model.constraints])
            row.update(maintenance_flux=flux[config['maintenance']], max_mass_residual=mass,
                       max_bound_violation=bound, max_constraint_violation=cons)
            write_json(out/(label+'.json'), {'result':row,'fluxes':flux,
                       'bounds':{r.id:list(r.bounds) for r in model.reactions}})
            save()
            if max(mass,bound,cons) > config['tolerance']:
                raise RuntimeError('Numerical validation failed: '+label)
        else:
            flux = None
        save()
        if row['status'] not in ('optimal','infeasible'):
            raise RuntimeError('Unresolved solver outcome: '+label)
        print(label, row['status'], row['objective'], flush=True)
        return row, flux

    def witness(blocks, label):
        model = closed.copy()
        for rid,direction in blocks:
            block_direction(model,rid,direction)
        model.reactions.get_by_id(config['maintenance']).bounds = (config['fixed_dissipation'],)*2
        model.objective = model.problem.Objective(sum(r.forward_variable+r.reverse_variable for r in model.reactions
                                          if r.id != config['maintenance']), direction='min')
        result, flux = solve(model, label, {'kind':'fixed_dissipation_L1', 'blocks':blocks,
                                'fixed_dissipation':config['fixed_dissipation'], 'minimizes':'sum |net flux| excluding maintenance'})
        if flux is None:
            return None
        active = [r for r in model.reactions if abs(flux[r.id]) > config['display_threshold']]
        for r in active:
            flux_rows.append({'witness':label,'reaction':r.id,'name':r.name,'equation':r.reaction,
                'named_equation':r.build_reaction_string(use_metabolite_names=True),'flux':flux[r.id],
                'actual_direction':'forward' if flux[r.id]>0 else 'reverse','bounds':list(r.bounds),
                'compartments':sorted(r.compartments),'gpr':r.gene_reaction_rule,
                'stoichiometry':{m.id:c for m,c in r.metabolites.items()},'annotation':r.annotation,'notes':r.notes,**balance(r)})
        maximum = 0.0
        for m in model.metabolites:
            contributions = {r.id:r.metabolites[m]*flux[r.id] for r in sorted(m.reactions,key=lambda r:r.id)
                             if r.id!=config['maintenance'] and flux[r.id]!=0}
            d = model.reactions.get_by_id(config['maintenance']).metabolites.get(m,0)*flux[config['maintenance']]
            net = sum(contributions.values())
            maximum=max(maximum,abs(net+d))
            if contributions or d:
                ledger.append({'witness':label,'metabolite':m.id,'name':m.name,'compartment':m.compartment,
                    'formula':m.formula,'charge':m.charge,'nonmaintenance_contributions':contributions,
                    'nonmaintenance_net':net,'maintenance_contribution':d,'total_residual':net+d})
        assert maximum <= config['tolerance']
        manifest['witnesses'].append({'label':label,'blocks':blocks,'nonzero_reactions':len(active),
                 'L1':result['objective'],'net_identity_max_residual':maximum,'strict_minimum_support':False})
        save()
        return [(r.id,'forward' if flux[r.id]>0 else 'reverse') for r in active if r.id!=config['maintenance']]

    try:
        maximum,_ = solve(closed, 'closed_max', {'kind':'max_maintenance','blocks':[]})
        assert maximum['status']=='optimal' and maximum['maintenance_flux']>config['tolerance']
        if maximum['maintenance_flux'] < config['fixed_dissipation']:
            raise RuntimeError('Positive maximum below fixed normalization; explicit smaller diagnostic value needed')
        directions = witness([], 'witness_1')
        assert directions
        alternative_blocks = []
        for rid,direction in directions:
            model=closed.copy()
            block_direction(model,rid,direction)
            result,_=solve(model,'block_'+rid+'_'+direction,{'kind':'direction_intervention','blocks':[[rid,direction]]})
            outcome = ('closed_model_infeasible' if result['status']=='infeasible' else
                       'ATP_within_tolerance' if result['maintenance_flux'] <= config['tolerance'] else 'alternative_ATP_path')
            interventions.append({'reaction':rid,'blocked_direction':direction,'blocks':[[rid,direction]],
                'status':result['status'],'max_ATP':result.get('maintenance_flux'),'outcome':outcome,
                'full_closed_model':True,'fixed_dissipation_removed':True})
            if outcome=='alternative_ATP_path':
                alternative_blocks.append([[rid,direction]])
            save()
        for i,blocks in enumerate(alternative_blocks[:config['limits']['max_witnesses']-1],2):
            witness(blocks, 'witness_'+str(i))
        if args.extra_cases:
            cases_path = workspace(args.extra_cases)
            manifest['extra_cases'] = {'path':str(cases_path),'sha256':sha(cases_path)}
            for case in json.loads(cases_path.read_text()):
                model = closed.copy()
                for rid,direction in case['blocks']:
                    block_direction(model,rid,direction)
                result,_ = solve(model,'extra_'+case['name'],{'kind':'direction_intervention','blocks':case['blocks']})
                outcome = ('closed_model_infeasible' if result['status']=='infeasible' else
                           'ATP_within_tolerance' if result['maintenance_flux']<=config['tolerance'] else 'alternative_ATP_path')
                interventions.append({'reaction':';'.join(r for r,d in case['blocks']),
                    'blocked_direction':';'.join(d for r,d in case['blocks']), 'blocks':case['blocks'],
                    'status':result['status'],'max_ATP':result.get('maintenance_flux'),'outcome':outcome,
                    'full_closed_model':True,'fixed_dissipation_removed':True})
                save()
                if case.get('extract_witness') and outcome=='alternative_ATP_path':
                    witness(case['blocks'],'witness_'+case['name'])
        manifest['status']='complete_initial_witness_and_interventions'
    except Exception as exc:
        manifest.update(status='incomplete',error=repr(exc))
        raise
    finally:
        manifest['baseline_restored']=signature(original)==protected
        manifest['inputs_unchanged']=all(sha(workspace(p))==h for p,h in manifest['input_sha256'].items())
        manifest['completed_utc']=datetime.now(timezone.utc).isoformat()
        save()
        assert manifest['baseline_restored'] and manifest['inputs_unchanged']


if __name__=='__main__':
    main()
