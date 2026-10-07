"""Build explicit, no-source Gly-Pro compartment hypotheses; never optimize."""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.diagnose_dipeptide_supply import balance, sha, table, write_json
from scripts.gem_annotate.config import REPO_ROOT
from scripts.gem_annotate.dipeptide_chemistry import protected_chemistry
from scripts.gem_annotate.energy_candidates import (
    export_candidate, model_definition, protected_definitions, signature,
    solver_definition,
)

SOURCE_SHA256 = 'a0307b9b00c1ed6e981d605f65ba15a0e313fe95e0e664dc05cfb95ef4291848'
HYDROLYSIS_ID = 'HYP_GLYPRO_HYD_CY'
SCENARIOS = ('TEMPLATE', 'V_ONLY', 'C_ONLY', 'NONE')
CAPACITY = 0.01
STOICHIOMETRY = {'m1865[C_cy]': -1, 'm32[C_cy]': -1,
                 'm272[C_cy]': 1, 'm765[C_cy]': 1}
NOTES = {
    'candidate_gene': 'YALI1E16433g',
    'candidate_gene_name': 'native_formal_name_unverified',
    'candidate_protein_function': 'M24B/X-Pro peptidase; prolidase-like candidate',
    'functional_status': 'hypothesis; target_free_Gly-L-Pro_activity_not_experimentally_confirmed',
    'localization_status': 'prediction_supported_not_experimental',
    'nuclear_localization': 'prediction_candidate_retained; not_modelled_in_this_comparison',
    'hypothesis_scope': 'diagnostic_compartment_comparison; not_formal_GPR_or_native_localization',
    'default_status': 'disabled; activate_only_in_explicit_scenario',
}
CONNECTIONS = {'R1363': (0., .04), 'R795': (0., .04),
               'R871': (0., .01), 'R876': (0., .01)}


def exact_copy(model):
    """Preserve coefficients and isolate metadata; reapply solver settings later."""
    model.solver.update()
    result = model.copy()
    for before, after in [(model, result), *(
            (obj, getattr(result, collection).get_by_id(obj.id))
            for collection in ('metabolites', 'reactions', 'genes', 'groups')
            for obj in getattr(model, collection))]:
        after.notes = copy.deepcopy(before.notes)
        after.annotation = copy.deepcopy(before.annotation)
    result._solver = model.problem.Model.from_json(model.solver.to_json())
    if solver_definition(result) != solver_definition(model):
        raise ValueError('JSON solver copy changed the exact mathematical definition')
    return result


def _workspace(path):
    path = Path(path).resolve()
    if not path.is_relative_to(REPO_ROOT.resolve()):
        raise ValueError('Candidate paths must remain in this project workspace')
    return path


def _metadata(model):
    """Include stored locks verbatim, not just their parsed meaning."""
    def fields(obj, names):
        return {key: copy.deepcopy(getattr(obj, key)) for key in names}
    return {
        'model': fields(model, ('id', 'name', 'notes', 'annotation', 'compartments')),
        'metabolites': {m.id: fields(m, ('name', 'formula', 'charge', 'compartment', 'notes', 'annotation'))
                        for m in model.metabolites},
        'reactions': {r.id: fields(r, ('name', 'subsystem', 'notes', 'annotation')) for r in model.reactions},
        'genes': {g.id: fields(g, ('name', 'notes', 'annotation')) for g in model.genes},
        'groups': {g.id: {**fields(g, ('name', 'kind', 'notes', 'annotation')),
                         'members': sorted((type(m).__name__, m.id) for m in g.members)} for g in model.groups},
    }


def equivalent_hydrolyses(model):
    """Detect proportional columns in either stored direction, not name matches."""
    found = []
    for reaction in model.reactions:
        values = {m.id: c for m, c in reaction.metabolites.items()}
        if set(values) != set(STOICHIOMETRY):
            continue
        ratios = {values[mid] / coefficient for mid, coefficient in STOICHIOMETRY.items()}
        if len(ratios) == 1:
            found.append({'reaction': reaction.id, 'stoichiometric_factor': ratios.pop()})
    return found


def scenario_bounds(baseline, scenario):
    if scenario not in SCENARIOS:
        raise ValueError(f'Unknown scenario: {scenario}')
    return {
        'R2039': baseline.reactions.R2039.bounds if scenario == 'TEMPLATE'
                  else ((0., CAPACITY) if scenario == 'V_ONLY' else (0., 0.)),
        HYDROLYSIS_ID: (0., CAPACITY) if scenario == 'C_ONLY' else (0., 0.),
    }


def verify_candidate(baseline, candidate, scenario):
    """Validate only declared overlays; restore them in a copy for strict locks.

    Historical chemistry locks include adjacency. They remain byte-identical in
    every model; this bounded validator does not weaken the production guard.
    """
    bounds = scenario_bounds(baseline, scenario)
    energy_locks = protected_definitions(baseline)
    chemistry_locks = protected_chemistry(baseline)
    if len(energy_locks) != 8 or len(chemistry_locks) != 13:
        raise ValueError('Expected eight E5 energy locks and thirteen chemical locks')
    if protected_definitions(candidate) != energy_locks:
        raise ValueError('E5 protected definitions differ')
    expected = model_definition(baseline)
    expected['reactions']['R2039']['bounds'] = list(bounds['R2039'])
    expected['reactions'][HYDROLYSIS_ID] = {
        'name': 'Gly-L-Pro hydrolysis in cytosol (compartment hypothesis)',
        'stoichiometry': STOICHIOMETRY, 'bounds': list(bounds[HYDROLYSIS_ID]), 'gpr': '',
    }
    if model_definition(candidate) != expected:
        raise ValueError('Candidate has an undeclared model-definition change')
    metadata = _metadata(candidate)
    new_metadata = metadata['reactions'].pop(HYDROLYSIS_ID)
    if metadata != _metadata(baseline):
        raise ValueError('Candidate changed original metadata, evidence, groups or gene records')
    if new_metadata != {'name': expected['reactions'][HYDROLYSIS_ID]['name'],
                        'subsystem': '', 'notes': NOTES, 'annotation': {}}:
        raise ValueError('New hypothesis metadata differ from the declared definition')
    for rid, connection_bounds in CONNECTIONS.items():
        if candidate.reactions.get_by_id(rid).bounds != connection_bounds:
            raise ValueError(f'Shared diagnostic connection changed: {rid}')
    result = balance(candidate.reactions.get_by_id(HYDROLYSIS_ID))
    if result['element_status'] != 'balanced' or result['charge_status'] != 'balanced_as_stored':
        raise ValueError('Cytosolic hydrolysis is not fully balanced')
    restored = exact_copy(candidate)
    restored.remove_reactions([HYDROLYSIS_ID])
    restored.reactions.R2039.bounds = baseline.reactions.R2039.bounds
    if protected_chemistry(restored) != chemistry_locks:
        raise ValueError('Restored historical chemistry/adjacency locks differ')
    if solver_definition(restored) != solver_definition(baseline):
        raise ValueError('Removing the declared overlay does not restore the actual solver definition')
    return {'scenario': scenario, 'energy_locks_preserved': len(energy_locks),
            'chemical_locks_preserved_verbatim': len(chemistry_locks),
            'historical_lock_validation_after_overlay_removal': True,
            'overlay_removal_restores_exact_solver_definition': True,
            'original_metadata_unchanged': True, 'no_artificial_sources_or_test_constraints': True,
            'hydrolysis_balance': result, 'declared_bounds': bounds}


def _negative_checks(baseline, template):
    tests = {}
    baseline_metadata, template_metadata = _metadata(baseline), _metadata(template)
    mutations = {
        'reject_changed_original_stoichiometry': lambda m: m.reactions.R2039.add_metabolites({m.metabolites.get_by_id('m1384[C_va]'): -1}),
        'reject_changed_chemical_field': lambda m: setattr(m.metabolites.get_by_id('m1865[C_cy]'), 'charge', 1),
        'reject_changed_evidence_note': lambda m: m.metabolites.get_by_id('m1865[C_cy]').notes.update({'unauthorized': 'yes'}),
        'reject_changed_energy_reaction': lambda m: setattr(m.reactions.get_by_id(next(iter(protected_definitions(m)))), 'upper_bound', 999.),
        'reject_undeclared_custom_constraint': lambda m: m.add_cons_vars(m.problem.Constraint(m.reactions.R2039.flux_expression, ub=0., name='UNDECLARED_TEST')),
    }
    for name, mutate in mutations.items():
        altered = exact_copy(template)
        mutate(altered)
        try:
            verify_candidate(baseline, altered, 'TEMPLATE')
        except ValueError as error:
            tests[name] = {'passed': True, 'reason': str(error)}
        else:
            raise AssertionError(f'Negative check failed: {name}')
    if _metadata(baseline) != baseline_metadata or _metadata(template) != template_metadata:
        raise AssertionError('Diagnostic copy mutations escaped to original metadata')
    tests['diagnostic_copy_metadata_isolation'] = {'passed': True}
    return tests


def build(source, source_sha256, output_dir):
    source, output_dir = _workspace(source), _workspace(output_dir)
    if source_sha256 != SOURCE_SHA256 or sha(source) != source_sha256:
        raise ValueError('Explicit source SHA differs from the audited chemical candidate')
    source_manifest_path = source.with_suffix('.build.json')
    source_manifest = json.loads(source_manifest_path.read_text())
    if (Path(source_manifest['output']).resolve() != source
            or source_manifest['output_sha256'] != source_sha256
            or source_manifest.get('candidate') != 'E5_vacuole_open_chemistry'):
        raise ValueError('Source build manifest does not identify the actual chemical candidate')
    parent_outputs = [output_dir.parent / name for name in
                      ('hypothesis_definitions.json', 'model_differences.tsv')]
    if (output_dir.is_symlink() or (output_dir.exists() and any(output_dir.iterdir()))
            or any(p.exists() or p.is_symlink() for p in parent_outputs)):
        raise FileExistsError('Refusing to overwrite candidates or hypothesis/difference records')
    baseline = read_sbml_model(source)
    existing = equivalent_hydrolyses(baseline)
    if existing or HYDROLYSIS_ID in baseline.reactions:
        # This builder is pinned to a reviewed source without an equivalent column.
        # A changed input must be reviewed for reuse, never duplicated silently.
        raise ValueError(f'Equivalent hydrolysis or ID collision; do not add a duplicate: {existing}')
    protected_chemistry(baseline)
    protected_definitions(baseline)
    ordinary_copy = baseline.copy()
    original_solver, ordinary_solver = solver_definition(baseline), solver_definition(ordinary_copy)
    copy_differences = [{'row': key, 'before': value, 'ordinary_copy': ordinary_solver['constraints'][key]}
                        for key, value in original_solver['constraints'].items()
                        if value != ordinary_solver['constraints'][key]]
    template = exact_copy(baseline)
    reaction = Reaction(HYDROLYSIS_ID, name='Gly-L-Pro hydrolysis in cytosol (compartment hypothesis)')
    reaction.add_metabolites({template.metabolites.get_by_id(mid): c for mid, c in STOICHIOMETRY.items()})
    reaction.bounds = (0., 0.)
    reaction.notes = copy.deepcopy(NOTES)
    template.add_reactions([reaction])
    verify_candidate(baseline, template, 'TEMPLATE')
    negative_checks = _negative_checks(baseline, template)
    output_dir.mkdir(parents=True, exist_ok=True)
    implementation = {str(p.relative_to(REPO_ROOT)): sha(p) for p in (
        Path(__file__).resolve(), REPO_ROOT/'scripts/gem_annotate/energy_candidates.py',
        REPO_ROOT/'scripts/gem_annotate/dipeptide_chemistry.py', REPO_ROOT/'scripts/gem_annotate/sbml.py')}
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(),
                'source': str(source), 'source_sha256': source_sha256,
                'source_build_manifest': str(source_manifest_path), 'source_build_manifest_sha256': sha(source_manifest_path),
                'implementation_sha256': implementation, 'no_optimization_performed': True,
                'copy_preflight': {'ordinary_gurobi_copy_exact': not bool(copy_differences),
                                   'changed_rows': copy_differences,
                                   'resolution': 'optlang JSON solver clone; exact actual solver definition checked'},
                'duplicate_check': {'equivalent_columns': existing, 'id_collision': False},
                'REF': {'path': str(source), 'sha256': source_sha256, 'immutable_source_reference': True},
                'candidates': {}, 'static_negative_checks': negative_checks}
    differences = []
    for scenario in SCENARIOS:
        model = exact_copy(template)
        for rid, limits in scenario_bounds(baseline, scenario).items():
            model.reactions.get_by_id(rid).bounds = limits
        verify_candidate(baseline, model, scenario)
        path = output_dir / f'{scenario}.xml'
        loaded = export_candidate(model, path)
        check = verify_candidate(baseline, loaded, scenario)
        record = {'scenario': scenario, 'source': str(source), 'source_sha256': source_sha256,
                  'source_build_manifest_sha256': sha(source_manifest_path),
                  'output': str(path), 'output_sha256': sha(path), 'implementation_sha256': implementation,
                  'scope': 'no-source diagnostic compartment hypothesis; not formal model or GPR acceptance',
                  'export_reload_definition_match': True, 'verification': check}
        write_json(path.with_suffix('.build.json'), record)
        manifest['candidates'][scenario] = {'path': str(path), 'sha256': sha(path),
            'build_manifest': str(path.with_suffix('.build.json')),
            'build_manifest_sha256': sha(path.with_suffix('.build.json'))}
        for rid in ('R2039', HYDROLYSIS_ID):
            old = signature(baseline.reactions.get_by_id(rid)) if rid in baseline.reactions else None
            after = signature(loaded.reactions.get_by_id(rid))
            if old != after:
                differences.append({'scenario': scenario, 'reaction_id': rid,
                    'change': 'added_hypothesis_column' if old is None else 'diagnostic_route_isolation_or_capacity',
                    'before': old, 'after': after, 'formal_biological_assignment': False})
    if sha(source) != source_sha256 or any(sha(REPO_ROOT / p) != value for p, value in implementation.items()):
        raise ValueError('Input or implementation changed during candidate construction')
    manifest['source_unchanged'] = True
    write_json(output_dir / 'manifest.json', manifest)
    write_json(parent_outputs[0], {
        'REF': 'Immutable original chemical candidate; no cytosolic hydrolysis added',
        'scenarios': {s: {'bounds': scenario_bounds(baseline, s),
                        'route_isolation_only_not_localization_evidence': True} for s in SCENARIOS},
        'hydrolysis_id': HYDROLYSIS_ID, 'matched_capacity_mmol_gdw_h': CAPACITY,
        'sources': {'O_CY': 'm1865[C_cy]', 'O_VA': 'm1866[C_va]'},
        'diagnostic_sources_persisted_in_xml': False, 'shared_transport_bounds': CONNECTIONS,
        'biological_evidence': NOTES,
        'protection_scope': 'Historical chemical and energy locks remain exact; the explicit overlay is validated separately and removed in a copy for strict historical adjacency validation',
        'candidate_manifest': str(output_dir / 'manifest.json')})
    table(parent_outputs[1], differences)
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--source-sha256', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    result = build(args.source, args.source_sha256, args.output_dir)
    print(json.dumps({'source_sha256': result['source_sha256'],
                      'candidates': result['candidates'], 'static_negative_checks': len(result['static_negative_checks'])}, indent=2))
