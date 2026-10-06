"""Apply the authorized two-coefficient R305 candidate; exactly four test LPs."""
import copy
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from cobra import Reaction
from cobra.io import read_sbml_model
from artifacts.coq_candidate_checks_20261005.run import exact_residual
from scripts.build_coq_biomass_candidate import pool_balance
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.execution import execution_limits
from scripts.validate_energy_candidates import (
    SolverBudget, load_model, close_model, balance, sha, write, table, signature,
    energy_verdict,
)

SPEC = ROOT / 'data/coq_r305_candidate.json'
NS = {'s': 'http://www.sbml.org/sbml/level3/version1/core',
      'h': 'http://www.w3.org/1999/xhtml'}


def build(spec, manifest):
    source = ROOT / spec['source']
    assert sha(source) == spec['source_sha256'], 'Control identity changed'
    control = read_sbml_model(source)
    q = control.metabolites.get_by_id('m468[C_mi]')
    assert control.reactions.biomass_C.metabolites[q] == -spec['alpha_mmol_per_gDW']
    pool = {'R385': 1.0, 'biomass_C': -spec['alpha_mmol_per_gDW']}
    assert pool_balance(control) == pool
    reaction = control.reactions.get_by_id(spec['reaction'])
    before = {m.id: c for m, c in reaction.metabolites.items()}
    assert all(before[k] == v for k, v in spec['before'].items())
    assert reaction.notes[spec['note_key']] == spec['note_before']
    text = source.read_text()
    expected_xml = ET.fromstring(text)
    xml_r = next(r for r in expected_xml.findall('.//s:reaction', NS)
                 if r.get('id') == 'R_R305')
    species = {'m28[C_mi]': ('listOfReactants', 'M_m28__91__C_mi__93__'),
               'm10[C_cy]': ('listOfProducts', 'M_m10__91__C_cy__93__')}
    match = re.search(r'<reaction\b[^>]*\bid="R_R305"[\s\S]*?</reaction>', text)
    assert match is not None
    chunk = match.group()
    for mid, final in spec['final_coefficients'].items():
        side, sid = species[mid]
        node = next(n for n in xml_r.findall('s:' + side + '/s:speciesReference', NS)
                    if n.get('species') == sid)
        assert float(node.get('stoichiometry')) == abs(spec['before'][mid])
        value = format(abs(final), 'g')  # SBML reactant magnitudes are positive.
        node.set('stoichiometry', value)
        pattern = r'(<speciesReference\b[^>]*species="' + re.escape(sid) + r'"[^>]*stoichiometry=")[^"]+("[^>]*/>)'
        chunk, count = re.subn(pattern, lambda m: m.group(1) + value + m.group(2), chunk)
        assert count == 1
    old_note = spec['note_key'] + ': ' + spec['note_before']
    new_note = spec['note_key'] + ': ' + spec['note_after']
    note = next(p for p in xml_r.findall('.//h:p', NS) if p.text == old_note)
    note.text = new_note
    assert chunk.count(old_note) == 1
    chunk = chunk.replace(old_note, new_note)
    edited = text[:match.start()] + chunk + text[match.end():]
    assert ET.tostring(ET.fromstring(edited)) == ET.tostring(expected_xml)
    path = OUT / 'E5_coq9_alpha_1e-4_R305_qcycle.xml'
    with path.open('x') as f:
        f.write(edited)
    candidate = read_sbml_model(path)
    actual = model_definition(candidate)
    expected = model_definition(control)
    expected['reactions']['R305']['stoichiometry'].update(spec['final_coefficients'])
    assert actual == expected, 'Non-target model definition changed'
    notes = {r.id: copy.deepcopy(r.notes) for r in control.reactions}
    notes['R305'][spec['note_key']] = spec['note_after']
    assert notes == {r.id: r.notes for r in candidate.reactions}
    assert pool_balance(candidate) == pool
    # Compare by identifiers, never cross-model Metabolite object identity.
    original = model_definition(control)
    rows = [{'reaction': rid, 'metabolite': mid, 'before': old.get(mid, 0), 'after': stoich.get(mid, 0)}
            for rid, record in actual['reactions'].items()
            for stoich, old in [(record['stoichiometry'], original['reactions'][rid]['stoichiometry'])]
            for mid in sorted(set(stoich) | set(old)) if stoich.get(mid, 0) != old.get(mid, 0)]
    assert len(rows) == 2 and {r['metabolite']: r['after'] for r in rows} == spec['final_coefficients']
    table(OUT / 'stoichiometric_diff.tsv', rows)
    r = candidate.reactions.R305
    chemistry = exact_residual(candidate, {m.id: c for m, c in r.metabolites.items()})
    manifest['build'] = {'control_preserved': True, 'control_alpha_confirmed_before_edit': 0.0001,
        'output': str(path.relative_to(ROOT)), 'output_sha256': sha(path),
        'stoichiometric_diff': rows, 'full_XML_other_fields_unchanged': True,
        'reloaded_definition_and_notes_match': True, 'R305_equation': r.reaction,
        'R305_named_equation': r.build_reaction_string(use_metabolite_names=True),
        'R305_chemistry': chemistry, 'combined_CoQ_pool_row': pool}
    write(OUT / 'manifest.json', manifest)
    print('Reloaded R305:', r.reaction, flush=True)
    print('Stoichiometric diff:', rows, flush=True)
    print('R305 exact balance:', chemistry, flush=True)
    if not chemistry['balanced']:
        raise ValueError('R305 balance failed; no optimizations permitted')
    assert sha(source) == spec['source_sha256']
    return path


def cases(model, config):
    for disabled in (False, True):
        with model:
            if disabled:
                model.reactions.R385.bounds = (0, 0)
            yield ('R385_off' if disabled else 'WT'), model
    closed, changes = close_model(model, config)
    for compartment in ('C_cy', 'C_mi'):
        with closed:
            if compartment == 'C_cy':
                drain = closed.reactions.get_by_id(config['maintenance'])
            else:
                drain = Reaction('DIAG_ATP_DISSIPATION_mi', lower_bound=0, upper_bound=1000)
                drain.add_metabolites({closed.metabolites.get_by_id(mid): c for mid, c in
                    {'m46[C_mi]': -1, 'm26[C_mi]': -1, 'm197[C_mi]': 1, 'm58[C_mi]': 1}.items()})
                closed.add_reactions([drain])
            chem = balance(drain)
            assert chem['element_status'] == 'balanced' and chem['charge_status'] == 'balanced_as_stored'
            assert drain.bounds == (0, 1000)
            closed.objective = drain
            closed.objective.direction = 'max'
            yield 'ATP_closed_' + compartment, closed


def main():
    assert not (OUT / 'manifest.json').exists(), 'Do not overwrite or repeat a run'
    spec = json.loads(SPEC.read_text())
    previous = ROOT / spec['previous_tests']
    prior = json.loads((previous / 'manifest.json').read_text())
    config = copy.deepcopy(prior['config'])
    config['limits'] = {'solves': 4, 'solve_wall_seconds': 240, 'per_solve_seconds': 60}
    protected = {str(SPEC): sha(SPEC), str(Path(__file__)): sha(__file__)}
    for path, digest in prior['protected_sha256'].items():
        assert sha(path) == digest, 'Previous input or helper changed: ' + path
        protected[path] = digest
    for path in [previous / 'manifest.json'] + sorted(previous.glob('candidate_*.json')):
        protected[str(path)] = sha(path)
    manifest = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
        'scope': spec['scope'], 'spec': spec, 'protected_sha256': protected,
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'git_dirty': subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True).splitlines(),
        'software': {'python': sys.version, **{p: version(p) for p in ('cobra', 'optlang', 'gurobipy')}},
        'control_reuse': {}, 'results': [], 'actual_optimization_calls': 0}
    write(OUT / 'manifest.json', manifest)
    budget = SolverBudget(OUT / 'solver_budget.json', config)
    try:
        with execution_limits(no_solve=True, allow_network=False):
            candidate = build(spec, manifest)
            control = load_model(ROOT / spec['source'], config)
            assert control.provenance() == prior['results']['candidate']['simulation_context']
            assert control.active_medium == prior['results']['candidate']['medium']
            assert control.strain_overlay_audit == prior['results']['candidate']['strain_overlay']
            for name, model in cases(control.model, config):
                old = json.loads((previous / ('candidate_' + name + '.json')).read_text())
                assert signature(model) == old['problem'], 'Control problem changed: ' + name
                assert old['result']['solver'] == config['solver']
                manifest['control_reuse'][name] = {'file': str((previous / ('candidate_' + name + '.json')).relative_to(ROOT)),
                    'exact_problem_and_solver_match': True, 'result': old['result'],
                    'R305_flux': old['fluxes']['R305'] if old['fluxes'] is not None else None,
                    'R385_flux': old['fluxes']['R385'] if old['fluxes'] is not None else None}
            sim = load_model(candidate, config)
            assert sim.active_medium == control.active_medium
            assert sim.strain_overlay_audit == control.strain_overlay_audit
            manifest['simulation_context'] = sim.provenance()
            manifest['medium'] = sim.active_medium
            manifest['strain_overlay'] = sim.strain_overlay_audit
            _, closure = close_model(sim.model, config)
            assert closure == prior['results']['candidate']['energy_closure']
            manifest['energy_closure'] = closure
            config['model'], config['model_sha256'] = str(candidate), sha(candidate)
            manifest['config'] = config
            protected[str(candidate)] = sha(candidate)
        write(OUT / 'manifest.json', manifest)
        for name, model in cases(sim.model, config):
            result, flux = budget.solve(model, name, OUT / (name + '.json'), {'case': name, 'variant': 'R305_only_candidate'})
            is_atp = name.startswith('ATP_')
            mu = flux['biomass_C'] if flux is not None else None
            row = {'case': name, 'status': result['status'], 'objective': result['objective'],
                'growth_h_inverse': mu, 'R305_flux': flux['R305'] if flux is not None else None,
                'R385_flux': flux['R385'] if flux is not None else None,
                'CoQ_pool_residual': flux['R385'] - 1e-4 * mu if flux is not None else None,
                'energy_verdict': energy_verdict(result['status'], result['objective'], config['tolerance']) if is_atp else None}
            manifest['results'].append(row)
            write(OUT / 'manifest.json', manifest)
            table(OUT / 'results.tsv', manifest['results'])
        assert len(budget.record['calls']) == 4
        assert all(sha(p) == h for p, h in protected.items()), 'Protected file changed during run'
        manifest.update(status='complete', protected_files_unchanged=True)
    except BaseException as exc:
        manifest.update(status='failed', error=repr(exc))
        raise
    finally:
        manifest['actual_optimization_calls'] = len(budget.record['calls'])
        write(OUT / 'manifest.json', manifest)


if __name__ == '__main__':
    with execution_limits(allow_network=False):
        main()
