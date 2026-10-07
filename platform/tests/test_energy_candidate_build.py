"""No-solve candidate build guards; optimization behavior is tested separately."""

import copy
import json
from pathlib import Path
import tempfile
import unittest

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.energy_candidates import (
    LOCK_NOTE, apply_energy_candidate, export_candidate, load_spec,
    model_definition, protected_definitions, signature, target_definition,
)
from scripts.gem_annotate.execution import execution_limits
from scripts.gem_annotate.reaction_selection import SELECTION_PATH, apply_metadata_reaction_selection

ROOT = Path(__file__).resolve().parents[1]


class EnergyCandidateBuildTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.guard = execution_limits(no_solve=True, allow_network=False)
        cls.attempts = cls.guard.__enter__()
        cls.base = read_sbml_model(ROOT / 'model_metadata_trna_r1159_leak.xml')
        cls.spec = load_spec()

    @classmethod
    def tearDownClass(cls):
        cls.guard.__exit__(None, None, None)
        assert cls.attempts == {'optimization': 0, 'network': 0}, cls.attempts

    def test_default_noop_and_variant_exact_scope(self):
        self.assertEqual(parse_args([]).energy_candidate, 'E0')
        model = self.base.copy(); before = model_definition(model)
        self.assertEqual(apply_energy_candidate(model)['status'], 'disabled')
        self.assertEqual(model_definition(model), before)
        for variant in ['E1', 'E2', 'E3', 'E4', 'E5']:
            with self.subTest(variant=variant):
                model = self.base.copy()
                apply_energy_candidate(model, variant)
                expected = copy.deepcopy(before)
                for rid, rule in self.spec['reactions'].items():
                    expected['reactions'][rid] = target_definition(rule, variant)
                self.assertEqual(model_definition(model), expected)
                if variant in ('E3', 'E4', 'E5'):
                    self.assertTrue(all(not model.reactions.get_by_id(rid).check_mass_balance()
                                        for rid in protected_definitions(model)))

    def test_cumulative_equivalence_idempotency_and_no_downgrade(self):
        direct = self.base.copy(); cumulative = self.base.copy()
        apply_energy_candidate(direct, 'E5')
        apply_energy_candidate(cumulative, 'E3'); apply_energy_candidate(cumulative, 'E4'); apply_energy_candidate(cumulative, 'E5')
        self.assertEqual(model_definition(direct), model_definition(cumulative))
        self.assertEqual(protected_definitions(direct), protected_definitions(cumulative))
        before = model_definition(cumulative)
        self.assertTrue(all(row['status'] == 'already_correct' for row in apply_energy_candidate(cumulative, 'E5')['items']))
        self.assertEqual(model_definition(cumulative), before)
        with self.assertRaisesRegex(ValueError, 'cannot remove'):
            apply_energy_candidate(cumulative, 'E1')
        self.assertEqual(model_definition(cumulative), before)

    def test_metadata_priority_and_export_reload(self):
        # The actual builder starts from SBML. Gurobi's LP-based model.copy()
        # rounds some tiny coefficients; the strict export guard rejects that drift.
        model = read_sbml_model(ROOT / 'model_metadata_trna_r1159_leak.xml')
        apply_energy_candidate(model, 'E5')
        spec = json.loads(SELECTION_PATH.read_text())
        spec['reactions'] = {rid: row for rid, row in spec['reactions'].items() if rid in self.spec['reactions']}
        before = model_definition(model)
        result = apply_metadata_reaction_selection(model, spec)
        self.assertTrue(result['complete'])
        self.assertEqual(model_definition(model), before)
        self.assertTrue(any(f['status'] == 'preserved_candidate' for row in result['items'] for f in row['fields']))
        with tempfile.TemporaryDirectory(dir=ROOT / 'artifacts/atp_candidate_repair_20260924') as directory:
            path = Path(directory) / 'candidate.xml'
            loaded = export_candidate(model, path)
            self.assertEqual(model_definition(loaded), before)
            self.assertEqual(protected_definitions(loaded), protected_definitions(model))
            with self.assertRaises(FileExistsError):
                export_candidate(model, path)

    def test_signature_conflicts_are_atomic(self):
        for kind in ['stoichiometry', 'gpr', 'species']:
            with self.subTest(kind=kind):
                model = self.base.copy(); r = model.reactions.R72
                if kind == 'stoichiometry':
                    r.add_metabolites({next(iter(r.metabolites)): .1})
                elif kind == 'gpr':
                    r.gene_reaction_rule = ''
                else:
                    model.metabolites.get_by_id('m170[C_cy]').charge = 999
                before = model_definition(model)
                with self.assertRaisesRegex(ValueError, 'signature differs'):
                    apply_energy_candidate(model, 'E4')
                self.assertEqual(model_definition(model), before)
                self.assertFalse(any(LOCK_NOTE in r.notes for r in model.reactions))

    def test_corrupt_persisted_lock_fails_before_metadata_edits(self):
        model = self.base.copy(); apply_energy_candidate(model, 'E4')
        model.reactions.R_NDP1.upper_bound = 1
        before = model_definition(model)
        with self.assertRaisesRegex(ValueError, 'protected energy candidate definition changed'):
            apply_metadata_reaction_selection(model)
        self.assertEqual(model_definition(model), before)

    def test_temporary_constraints_and_mass_row_edits_rejected(self):
        model = self.base.copy(); apply_energy_candidate(model, 'E4')
        with tempfile.TemporaryDirectory(dir=ROOT / 'artifacts/atp_candidate_repair_20260924') as directory:
            path = Path(directory) / 'invalid.xml'
            model.add_cons_vars(model.problem.Constraint(model.reactions.R72.flux_expression, ub=0, name='diagnostic_only'))
            with self.assertRaisesRegex(ValueError, 'custom constraints'):
                export_candidate(model, path)
            self.assertFalse(path.exists())
            model.remove_cons_vars(model.constraints.diagnostic_only)
            row = model.constraints[0]; row.ub = 1
            with self.assertRaisesRegex(ValueError, 'steady-state row bounds'):
                export_candidate(model, path)
            self.assertFalse(path.exists())
            model = read_sbml_model(ROOT / 'model_metadata_trna_r1159_leak.xml')
            apply_energy_candidate(model, 'E5')
            model.constraints['m170[C_cy]'].set_linear_coefficients({model.reactions.R72.forward_variable: -1.1})
            with self.assertRaisesRegex(ValueError, 'actual solver definition'):
                export_candidate(model, path)


if __name__ == '__main__':
    unittest.main()
