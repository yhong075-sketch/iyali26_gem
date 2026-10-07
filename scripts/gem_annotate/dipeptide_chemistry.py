"""Opt-in neutral-parent metadata candidate; chemistry never changes the LP."""
import copy
import hashlib
import html
import json
from pathlib import Path

from cobra.io import read_sbml_model

from .config import REPO_ROOT
from .energy_candidates import export_candidate, model_definition, protected_definitions, signature, solver_definition

SPEC_PATH = REPO_ROOT / 'data/dipeptide_chemistry_patch.json'
LOCK_NOTE = 'dipeptide_chemistry_protected_definition'
FIELDS = ('name', 'formula', 'charge', 'compartment', 'annotation', 'notes')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def species_fields(met):
    result = {key: copy.deepcopy(getattr(met, key)) for key in FIELDS}
    result['notes'].pop(LOCK_NOTE, None)
    result['notes'] = {key: html.unescape(str(value)) for key, value in result['notes'].items()}
    return result


def adjacency(met):
    return {r.id: signature(r) for r in sorted(met.reactions, key=lambda r: r.id)}


def protected_chemistry(model):
    """Preflight persisted fields and their adjacent reactions before any metadata edit."""
    locks = {}
    for met in model.metabolites:
        if LOCK_NOTE not in met.notes:
            continue
        lock = json.loads(html.unescape(str(met.notes[LOCK_NOTE])))
        if (lock.get('schema_version') != 1 or lock.get('id') != met.id
                or lock.get('fields') != species_fields(met)
                or lock.get('adjacency') != adjacency(met)):
            raise ValueError(f'{met.id}: protected dipeptide chemistry or adjacency changed')
        locks[met.id] = lock
    return locks


def protected_reactions(locks):
    return {rid for lock in locks.values() for rid in lock['adjacency']}


def apply_chemistry_candidate(model, enabled=False, spec=None):
    if not enabled:
        return {'status': 'disabled', 'items': []}
    spec = json.loads(SPEC_PATH.read_text()) if spec is None else spec
    if spec.get('schema_version') != 1 or spec.get('enabled_by_default') is not False:
        raise ValueError('Invalid opt-in chemistry specification')
    if digest(solver_definition(model)) != spec['mathematical_sha256']:
        raise ValueError('Chemistry input mathematical definition differs')
    if protected_definitions(model) != spec['energy_protected_definitions']:
        raise ValueError('Eight E5 energy locks differ')
    protected_chemistry(model)
    planned = []
    for row in spec['entries']:
        met = model.metabolites.get_by_id(row['id'])
        current = species_fields(met)
        if current not in (row['before'], row['after']) or adjacency(met) != row['adjacency']:
            raise ValueError(f'{met.id}: chemistry precondition differs; no edits applied')
        if set(row['before']) != set(FIELDS) or set(row['after']) != set(FIELDS):
            raise ValueError('Only declared chemical metadata fields may change')
        if row['before']['compartment'] != row['after']['compartment']:
            raise ValueError('Compartment changes are outside this candidate')
        planned.append((met, row, current))
    # Every input was checked before the first mutation; absolute assignment is idempotent.
    result = []
    for met, row, current in planned:
        for key, value in row['after'].items():
            setattr(met, key, copy.deepcopy(value))
        lock = {'schema_version': 1, 'id': met.id, 'curation_id': spec['curation_id'],
                'fields': row['after'], 'adjacency': row['adjacency']}
        met.notes[LOCK_NOTE] = json.dumps(lock, sort_keys=True)
        result.append({'id': met.id, 'status': 'already_correct' if current == row['after'] else 'applied',
                       'evidence_status': row['evidence_status']})
    protected_chemistry(model)
    return {'status': 'complete', 'items': result}


def build_candidate_file(source, source_sha256, output, enabled=False, spec_path=SPEC_PATH):
    source, output, spec_path = map(Path, (source, output, spec_path))
    if not all(p.resolve().is_relative_to(REPO_ROOT) for p in (source, output, spec_path)):
        raise ValueError('Chemistry build paths must remain inside project workspace')
    manifest = output.with_suffix('.build.json')
    if any(p.exists() or p.is_symlink() for p in (output, manifest)):
        raise FileExistsError('Refusing to overwrite candidate or manifest')
    spec = json.loads(spec_path.read_text())
    actual_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    if actual_sha != source_sha256 or actual_sha != spec['source_sha256']:
        raise ValueError('Explicit chemistry input SHA differs')
    model = read_sbml_model(source)
    before, math_before = model_definition(model), solver_definition(model)
    met_before = {m.id: species_fields(m) for m in model.metabolites}
    reaction_meta = {r.id: [r.annotation, r.notes] for r in model.reactions}
    result = apply_chemistry_candidate(model, enabled, spec)
    expected = copy.deepcopy(before)
    expected_meta = copy.deepcopy(met_before)
    if enabled:
        for row in spec['entries']:
            expected['metabolites'][row['id']] = {k: row['after'][k] for k in before['metabolites'][row['id']]}
            expected_meta[row['id']] = row['after']
    if model_definition(model) != expected or solver_definition(model) != math_before:
        raise ValueError('Chemistry patch changed undeclared mathematical or species fields')
    loaded = export_candidate(model, output)
    if ({m.id: species_fields(m) for m in loaded.metabolites} != expected_meta
            or {r.id: [r.annotation, r.notes] for r in loaded.reactions} != reaction_meta
            or protected_chemistry(loaded) != protected_chemistry(model)):
        raise ValueError('Export/reload did not preserve complete chemistry metadata')
    if hashlib.sha256(source.read_bytes()).hexdigest() != actual_sha:
        raise ValueError('Source changed during candidate build')
    implementations = [Path(__file__), REPO_ROOT/'scripts/build_dipeptide_chemistry.py',
        *(Path(__file__).with_name(name+'.py') for name in
          ('energy_candidates', 'metabolites', 'microspecies', 'reaction_selection', 'sbml'))]
    record = {'candidate': 'E5_vacuole_open_chemistry', 'enabled': enabled, 'enabled_by_default': False,
        'source': str(source.resolve()), 'source_sha256': actual_sha,
        'output': str(output.resolve()), 'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
        'spec_path': str(spec_path.resolve()), 'spec_sha256': hashlib.sha256(spec_path.read_bytes()).hexdigest(),
        'implementation_sha256': {str(p.relative_to(REPO_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in implementations},
        'patch': result, 'mathematical_sha256': digest(math_before), 'mathematical_definition_unchanged': True,
        'all_reaction_fields_and_notes_unchanged': True, 'non_target_metabolites_unchanged': True,
        'energy_locks_unchanged': True, 'export_reload_definition_match': True, 'source_unchanged': True,
        'scope': 'chemical metadata candidate; twelve conditional dipeptide mappings, not native function validation'}
    manifest.write_text(json.dumps(record, indent=2)+'\n')
    return record
