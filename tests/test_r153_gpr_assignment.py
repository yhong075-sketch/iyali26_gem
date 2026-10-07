"""Remove the placeholder while preserving contrary biological evidence."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.patches import R153_ASSIGNMENT_PATH, apply_r153_gpr_assignment
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

ROOT = Path(__file__).resolve().parents[1]


class R153AssignmentTests(unittest.TestCase):
    def test_single_gene_knockout_roundtrip_and_conflict_guards(self):
        model = read_sbml_model(str(ROOT / 'model_metadata_trna_ntp1_hydrolysis.xml'))
        before, original_notes = semantics(model), annotations(model)
        self.assertEqual(apply_r153_gpr_assignment(model)['status'], 'applied')
        rx = model.reactions.R153
        self.assertEqual(rx.gene_reaction_rule, 'YALI1D17462g')
        self.assertNotIn('YALIUNK2', model.genes)
        self.assertEqual(len(model.genes), len(before['genes']) - 1)
        self.assertEqual(rx.notes['gpr_evidence_label'], '用户指定赋值；功能证据冲突，未实验确认')
        self.assertEqual(rx.notes['gpr_experimental_confirmation'], 'required_not_confirmed')
        original_bounds = {r.id: r.bounds for r in model.reactions}
        for deleted, closed in [((), set()), (('YALI1D17462g',), {'R153', 'R2176'})]:
            with self.subTest(deleted=deleted), model:
                for gid in deleted:
                    model.genes.get_by_id(gid).knock_out()
                self.assertEqual({r.id: r.bounds for r in model.reactions
                                  if r.bounds != original_bounds[r.id]},
                                 {rid: (0, 0) for rid in closed})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'r153.xml'
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual(loaded.reactions.R153.notes, rx.notes)
        self.assertNotIn('YALIUNK2', loaded.genes)
        self.assertEqual(apply_r153_gpr_assignment(loaded)['status'], 'already_correct')
        loaded.reactions.R153.gene_reaction_rule = 'YALIUNK2'
        loaded.reactions.R153.notes = copy.deepcopy(original_notes['R153'][2])
        self.assertEqual(semantics(loaded), before)
        self.assertEqual(annotations(loaded), original_notes)

        for conflict in ('gpr', 'bounds', 'species', 'gene_identity', 'evidence',
                         'scope', 'removal_scope', 'other_placeholder_reaction'):
            bad, spec = model.copy(), json.loads(R153_ASSIGNMENT_PATH.read_text())
            if conflict == 'gpr':
                bad.reactions.R153.gene_reaction_rule = 'YALIUNK2 or YALI1D17462g'
            elif conflict == 'bounds':
                bad.reactions.R153.lower_bound = 0
            elif conflict == 'species':
                next(iter(bad.reactions.R153.metabolites)).compartment = 'C_nu'
            elif conflict == 'gene_identity':
                bad.genes.YALI1D17462g.annotation['uniprot'] = 'conflicting_accession'
            elif conflict == 'evidence':
                bad.reactions.R153.notes['gpr_experimental_confirmation'] = 'experimentally_verified'
            elif conflict == 'scope':
                spec['after_gpr'] = 'YALIUNK2 or YALI1D17462g'
            elif conflict == 'removal_scope':
                spec['remove_orphan_gene'] = 'YALI1D17462g'
            else:
                bad.reactions.R2176.gene_reaction_rule = 'YALI1D17462g or YALIUNK2'
            state = semantics(bad), annotations(bad)
            with self.subTest(conflict=conflict), self.assertRaises(ValueError):
                apply_r153_gpr_assignment(bad, spec)
            self.assertEqual((semantics(bad), annotations(bad)), state)


if __name__ == '__main__':
    unittest.main()
