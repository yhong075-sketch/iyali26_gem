"""Metadata candidate and protection checks; optimization/network forbidden."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from cobra.io import read_sbml_model
from scripts.gem_annotate.dipeptide_chemistry import (
    SPEC_PATH, apply_chemistry_candidate, build_candidate_file, species_fields,
    protected_chemistry, protected_reactions,
)
from scripts.gem_annotate.energy_candidates import model_definition, solver_definition
from scripts.gem_annotate.execution import execution_limits
from scripts.gem_annotate.metabolites import annotate_metabolites
from scripts.gem_annotate.microspecies import apply_curated_microspecies, balance_protons_and_water
from scripts.gem_annotate.reaction_selection import apply_metadata_reaction_selection, reaction_fields

ROOT = Path(__file__).resolve().parents[1]
TASK = ROOT/'artifacts/dipeptide_chemistry_routes_20260924'
SOURCE = TASK/'E5_vacuole_open_rebuilt.xml'


class ChemistryCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.guard = execution_limits(no_solve=True, allow_network=False)
        cls.attempts = cls.guard.__enter__()
        cls.spec = json.loads(SPEC_PATH.read_text())

    @classmethod
    def tearDownClass(cls):
        cls.guard.__exit__(None, None, None)
        assert cls.attempts == {'optimization': 0, 'network': 0}

    def test_default_atomic_conflict_idempotency_and_unchanged_lp(self):
        model = read_sbml_model(SOURCE)
        before, lp = model_definition(model), solver_definition(model)
        self.assertEqual(apply_chemistry_candidate(model)['status'], 'disabled')
        self.assertEqual(model_definition(model), before)
        # A conflict in the LAST item must prevent changes to the first twelve.
        gly = model.metabolites.get_by_id('m1863[C_va]')
        gly.annotation['unexpected'] = 'conflict'
        snapshot = {m.id: species_fields(m) for m in model.metabolites}
        with self.assertRaisesRegex(ValueError, 'precondition'):
            apply_chemistry_candidate(model, True)
        self.assertEqual({m.id: species_fields(m) for m in model.metabolites}, snapshot)
        gly.annotation.pop('unexpected')
        apply_chemistry_candidate(model, True)
        self.assertEqual(solver_definition(model), lp)
        self.assertEqual(len(protected_chemistry(model)), 13)
        self.assertTrue(all(x['status'] == 'already_correct' for x in apply_chemistry_candidate(model, True)['items']))
        with model:
            model.reactions.R795.upper_bound = 1
            with self.assertRaisesRegex(ValueError, 'mathematical'):
                apply_chemistry_candidate(model, True)

    def test_metadata_and_auto_balance_cannot_overwrite_candidate(self):
        model = read_sbml_model(SOURCE)
        apply_chemistry_candidate(model, True)
        locked = protected_chemistry(model)
        lp = solver_definition(model)
        annotate_metabolites(model, {'by_source': {}, 'by_mnxid': {}},
            {'prop': {'fake': {'name': 'Gly-L-Asp', 'formula': 'C999', 'charge': 7}},
             'name_index': {'gly-l-asp': 'fake'}})
        self.assertEqual(protected_chemistry(model), locked)
        r = model.reactions.R2021
        before = reaction_fields(r); after = copy.deepcopy(before); after['bounds'] = [0, 0]
        spec = {'selection_id': 'conflicting_future_metadata', 'source_metadata_sha256': 'test',
                'field_counts': {}, 'reactions': {r.id: {'before': before, 'after': after,
                    'fields': ['bounds'], 'species': {}}}}
        result = apply_metadata_reaction_selection(model, spec)
        self.assertEqual(result['items'][0]['status'], 'preserved_candidate')
        ids = protected_reactions(locked)
        report = balance_protons_and_water(model, ids)
        self.assertEqual(set(report['skipped_curated_lock_reaction_ids']), ids)
        self.assertEqual(report['changes'], [])
        self.assertEqual(solver_definition(model), lp)
        gly = model.metabolites.get_by_id('m1863[C_va]')
        row = SimpleNamespace(status='active', family_id='test', target_formula='C2H4NO2',
                              target_charge=-1, allowed_current_pairs={(gly.formula, gly.charge)})
        with patch('scripts.gem_annotate.microspecies.load_curated_microspecies', return_value=[row]), \
             patch('scripts.gem_annotate.microspecies._resolve_pinned_targets', return_value=[gly]), \
             patch('scripts.gem_annotate.microspecies._is_allowed_current_pair', return_value=True), \
             self.assertRaisesRegex(ValueError, 'separate microspecies review'):
            apply_curated_microspecies(model)
        self.assertEqual(protected_chemistry(model), locked)

    def test_export_and_explicit_identity_failures(self):
        with tempfile.TemporaryDirectory(dir=TASK/'tmp') as directory:
            output = Path(directory)/'candidate.xml'
            with self.assertRaisesRegex(ValueError, 'SHA'):
                build_candidate_file(SOURCE, '0'*64, output, True)
            self.assertFalse(output.exists())
            record = build_candidate_file(SOURCE, self.spec['source_sha256'], output, True)
            self.assertTrue(record['mathematical_definition_unchanged'])
            loaded = read_sbml_model(output)
            self.assertEqual(len(protected_chemistry(loaded)), 13)
            for row in self.spec['entries']:
                self.assertEqual(species_fields(loaded.metabolites.get_by_id(row['id'])), row['after'])
            with self.assertRaises(FileExistsError):
                build_candidate_file(SOURCE, self.spec['source_sha256'], output, True)


if __name__ == '__main__':
    unittest.main()
