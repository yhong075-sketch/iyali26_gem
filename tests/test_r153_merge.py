"""Authorized duplicate removal, capacity, preservation and conflict guard."""
import copy
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model
from scripts.gem_annotate.patches import merge_r153_r2176
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics, annotations

ROOT = Path(__file__).resolve().parents[1]


class R153MergeTests(unittest.TestCase):
    def test_merge_roundtrip_idempotence_and_conflicts(self):
        original = read_sbml_model(str(ROOT / 'model_metadata_trna_r153_single_gpr.xml'))
        for conflict in ('bounds', 'gpr', 'stoichiometry', 'notes', 'constraint'):
            bad = original.copy()
            r = bad.reactions.R2176
            if conflict == 'bounds':
                r.upper_bound = 999
            elif conflict == 'gpr':
                r.gene_reaction_rule = ''
            elif conflict == 'stoichiometry':
                r.add_metabolites({next(iter(r.metabolites)): 1})
            elif conflict == 'notes':
                r.notes['unreviewed'] = 'new source'
            else:
                bad.add_cons_vars(bad.problem.Constraint(r.flux_expression, ub=1, name='extra'))
            before = semantics(bad), annotations(bad)
            with self.subTest(conflict=conflict), self.assertRaises(ValueError):
                merge_r153_r2176(bad)
            self.assertEqual((semantics(bad), annotations(bad)), before)

        pre_export = original.copy()
        for rid in ('R153', 'R2176'):
            reaction = pre_export.reactions.get_by_id(rid)
            reaction.annotation = {k: [v] if isinstance(v, str) and k != 'sbo' else v
                                   for k, v in reaction.annotation.items()}
        self.assertEqual(merge_r153_r2176(pre_export)['status'], 'applied')

        # Gurobi copy changes bound string forms (0 vs 0.0); compare direct loads.
        merged = read_sbml_model(str(ROOT / 'model_metadata_trna_r153_single_gpr.xml'))
        self.assertEqual(merge_r153_r2176(merged)['status'], 'applied')
        self.assertNotIn('R2176', merged.reactions)
        self.assertEqual(merged.reactions.R153.bounds, (-1000, 1000))
        self.assertEqual(merged.reactions.R153.notes['gpr_evidence_label'],
                         original.reactions.R153.notes['gpr_evidence_label'])
        expected = read_sbml_model(str(ROOT / 'model_metadata_trna_r153_single_gpr.xml'))
        expected.remove_reactions(['R2176'], remove_orphans=False)
        expected.reactions.R153.notes = copy.deepcopy(merged.reactions.R153.notes)
        self.assertTrue(semantics(merged) == semantics(expected), 'Unexpected mathematical change')
        self.assertEqual(annotations(merged), annotations(expected))
        for group in merged.groups:
            old = {r.id for r in original.groups.get_by_id(group.id).members}
            self.assertEqual({r.id for r in group.members},
                             (old - {'R2176'}) | ({'R153'} if 'R2176' in old else set()))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'merged.xml'
            write_deterministic_sbml_model(merged, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(semantics(loaded) == semantics(merged), 'SBML mathematical roundtrip differs')
        self.assertTrue(annotations(loaded) == annotations(merged), 'SBML annotation roundtrip differs')
        self.assertEqual(merge_r153_r2176(loaded)['status'], 'already_correct')


if __name__ == '__main__':
    unittest.main()
