"""Temporary R795/R1363 opening comparison; no model is exported."""
import csv
import hashlib
import json
import math
import os
import sys
import time
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

base = Path(__file__).resolve().parent
root = base.parents[1]
code = root / 'artifacts/r608_engineering_20260907/code'
research = root / 'artifacts/reference_pipeline_restore_20260909/research'
os.environ['IYALI26_RESEARCH_ROOT'] = str(research)
sys.path.insert(0, str(code))
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
old = json.loads((root / 'artifacts/vatpase_gpr_hypothesis_20260911/screen/run_manifest.json').read_text())
for field in ('model', 'medium', 'strain_profile'):
    assert sha(Path(old[field]['path'])) == old[field]['sha256']
simulation = load_effective_simulation_context(model_path=Path(old['model']['path']),
    media_path=Path(old['medium']['path']), strain_profile_path=Path(old['strain_profile']['path']))
m = simulation.model
m.solver = 'gurobi'
m.solver.problem.Params.Threads = 1
m.solver.problem.Params.TimeLimit = 60
assert m.solver.configuration.presolve is False
assert simulation.active_medium == old['medium']['active_medium']
assert simulation.strain_overlay_audit == old['strain_profile']['audit']
assert m.reactions.R795.bounds == (0, 0)
assert m.reactions.R1363.bounds == (0, 0)
previous = json.loads((root / 'artifacts/r795_open_diagnostic_20260911/results.json').read_text())
for name, digest in previous['source_sha256'].items():
    assert sha(root / name) == digest, name
required = ['YALI1D00581g', 'YALI0E16192g', 'YALI1F38820g']
start = time.monotonic()
record = {'started_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'running',
          'input': old['model'], 'medium': old['medium'], 'strain_profile': old['strain_profile'],
          'simulation_context': simulation.provenance(),
          'solver': {'interface': m.solver.interface.__name__, 'presolve': m.solver.configuration.presolve,
                     'Threads': m.solver.problem.Params.Threads, 'TimeLimit': m.solver.problem.Params.TimeLimit,
                     'tolerances': {k: getattr(m.solver.configuration.tolerances, k) for k in ('feasibility', 'optimality', 'integrality')}},
          'software': {'python': sys.version, **{k: metadata.version(k) for k in ('cobra', 'optlang', 'gurobipy')}},
          'source_sha256': {str(Path(mod.__file__).relative_to(root)): sha(Path(mod.__file__))
                            for mod in list(sys.modules.values()) if getattr(mod, '__file__', None)
                            and Path(mod.__file__).is_relative_to(code) and Path(mod.__file__).is_file()},
          'script_sha256': sha(Path(__file__)), 'primary_solves': 0, 'runs': {},
          'fva_biomass_absolute_relaxation': 1e-8, 'temporary_open_bounds': {'R795': [0, 1000], 'R1363': [0, 1000]}, 'gene_ids': required}
def save():
    record['elapsed_seconds'] = time.monotonic() - start
    (base / 'results.json').write_text(json.dumps(record, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
def solve(label, opened=(), gene=None, flux_target=None, flux_direction=None, biomass_floor=None):
    assert record['primary_solves'] < 12 and time.monotonic() - start < 540
    with m:
        for rid in opened:
            assert rid in ('R795', 'R1363')
            m.reactions.get_by_id(rid).bounds = (0, 1000)
        if gene:
            m.genes.get_by_id(gene).knock_out()
        if biomass_floor is not None:
            m.reactions.biomass_C.lower_bound = biomass_floor
        if flux_direction:
            m.objective = m.reactions.get_by_id(flux_target)
            m.objective.direction = flux_direction
        record['primary_solves'] += 1
        solution = m.optimize()
        row = {'opened': opened, 'gene_knockout': gene, 'objective': flux_target if flux_direction else 'biomass_C',
               'direction': m.objective.direction, 'biomass_floor': biomass_floor,
               'raw_status': str(solution.status),
               'raw_objective': float(solution.objective_value) if solution.objective_value is not None and math.isfinite(float(solution.objective_value)) else None,
               'R794_bounds': list(m.reactions.R794.bounds), 'R795_bounds': list(m.reactions.R795.bounds), 'R1363_bounds': list(m.reactions.R1363.bounds)}
        record['runs'][label] = row
        save()
        assert solution.status == 'optimal' and row['raw_objective'] is not None, row
        row.update(biomass=float(solution.fluxes['biomass_C']), R794=float(solution.fluxes['R794']), R795=float(solution.fluxes['R795']), R1363=float(solution.fluxes['R1363']))
        row['max_mass_balance_residual'] = max(abs(sum(r.metabolites[met] * float(solution.fluxes[r.id]) for r in met.reactions)) for met in m.metabolites)
        assert all(math.isfinite(float(v)) for v in solution.fluxes)
        row['max_bound_violation'] = max(max(r.lower_bound - float(solution.fluxes[r.id]), float(solution.fluxes[r.id]) - r.upper_bound, 0) for r in m.reactions)
        assert row['max_mass_balance_residual'] <= 1e-7 and row['max_bound_violation'] <= 1e-7, row
        row['vacuolar_proton_rates'] = {r.id: {'flux': float(solution.fluxes[r.id]), 'coefficient': r.metabolites[m.metabolites.get_by_id('m1007[C_va]')]}
                                       for r in sorted(m.metabolites.get_by_id('m1007[C_va]').reactions, key=lambda r: r.id)}
        with (base / (label + '_fluxes.tsv')).open('w') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['reaction_id', 'name', 'lower_bound', 'upper_bound', 'flux'])
            for r in m.reactions:
                writer.writerow([r.id, r.name, r.lower_bound, r.upper_bound, repr(float(solution.fluxes[r.id]))])
        save()
        print(label, row['raw_status'], 'biomass', row['biomass'], 'R795', row['R795'], 'R1363', row['R1363'], flush=True)
        return row

try:
    closed = solve('both_closed_WT')
    pump_only = solve('only_R795_open_WT', ('R795',))
    water_only = solve('only_R1363_open_WT', ('R1363',))
    both = ('R795', 'R1363')
    opened = solve('both_open_WT', both)
    assert opened['biomass'] > 0
    assert math.isclose(closed['biomass'], previous['runs']['closed_WT']['biomass'], rel_tol=1e-6, abs_tol=1e-8)
    for gene in required:
        row = solve('KO_' + gene, both, gene=gene)
        row['KO_WT_ratio'] = row['biomass'] / opened['biomass']
        row['essential_at_cutoff'] = {str(c): row['KO_WT_ratio'] < c for c in (.01, .05, .10, .15)}
        save()
    floor = opened['biomass'] - record['fva_biomass_absolute_relaxation']
    solve('R795_min_at_optimal_growth', both, flux_target='R795', flux_direction='min', biomass_floor=floor)
    solve('R795_max_at_optimal_growth', both, flux_target='R795', flux_direction='max', biomass_floor=floor)
    solve('R795_max_without_growth_floor', both, flux_target='R795', flux_direction='max')
    solve('R1363_max_at_optimal_growth', both, flux_target='R1363', flux_direction='max', biomass_floor=floor)
    solve('R1363_max_without_growth_floor', both, flux_target='R1363', flux_direction='max')
    assert record['primary_solves'] == 12
    record['input_files_unchanged'] = all(sha(Path(old[f]['path'])) == old[f]['sha256'] for f in ('model', 'medium', 'strain_profile'))
    before = json.loads((base / 'before.json').read_text())
    record['protected_files_unchanged'] = all(sha(root / name) == digest for name, digest in before['protected_sha256'].items())
    assert record['input_files_unchanged'] and record['protected_files_unchanged']
    record['status'] = 'complete'
except Exception as error:
    record.update(status='incomplete', error=repr(error))
    raise
finally:
    save()
