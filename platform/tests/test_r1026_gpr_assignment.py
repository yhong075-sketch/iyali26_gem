"""Check the real source reaction, evidence roundtrip and fail-before-edit guards."""
import copy
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.patches import apply_r1026_gpr_assignment
from scripts.gem_annotate.metabolites import normalize_all_annotations
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics

ROOT = Path(__file__).resolve().parents[2]


from scripts.gem_annotate.model_layout import MODEL
class R1026AssignmentTests(unittest.TestCase):
    def test_assignment_roundtrip_and_conflicts(self):
        model = read_sbml_model(str(MODEL.candidate_file('model_metadata_trna.xml')))
        before = semantics(model)
        rx = model.reactions.R1026
        original_notes = dict(rx.notes)
        normalized = model.copy()
        normalize_all_annotations(normalized)
        self.assertEqual(apply_r1026_gpr_assignment(normalized)['status'], 'applied')
        self.assertEqual(apply_r1026_gpr_assignment(model)['status'], 'applied')
        self.assertEqual(rx.gene_reaction_rule, 'YALI1F28274g')
        self.assertEqual({k: rx.notes[k] for k in original_notes}, original_notes)
        with model:
            model.genes.YALI1F28274g.knock_out()
            self.assertEqual(model.reactions.R1026.bounds, (0, 0))
            self.assertEqual(model.reactions.R2202.bounds, (0, 0))
            self.assertEqual(model.reactions.R1025.bounds, (-1000, 1000))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'r1026.xml'
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual(loaded.reactions.R1026.notes, rx.notes)
        self.assertEqual(apply_r1026_gpr_assignment(loaded)['status'], 'already_correct')
        loaded.reactions.R1026.gene_reaction_rule = ''
        self.assertEqual(semantics(loaded), before)

        for conflict in ('gpr', 'counterpart', 'bound', 'chemistry', 'gene_identity', 'notes'):
            bad = model.copy()
            if conflict == 'gpr':
                bad.reactions.R1026.gene_reaction_rule = 'YALI1A08512g'
            elif conflict == 'counterpart':
                bad.reactions.R2202.gene_reaction_rule = ''
            elif conflict == 'bound':
                bad.reactions.R1026.upper_bound = 999
            elif conflict == 'chemistry':
                bad.metabolites.get_by_id('m199[C_cy]').charge = 0
            elif conflict == 'gene_identity':
                bad.genes.YALI1F28274g.annotation['uniprot'] = 'conflicting_accession'
            else:
                bad.reactions.R1026.notes['gpr_evidence_status'] = 'experimentally_verified'
            state = semantics(bad)
            notes = copy.deepcopy(bad.reactions.R1026.notes)
            with self.subTest(conflict=conflict), self.assertRaises(ValueError):
                apply_r1026_gpr_assignment(bad)
            self.assertEqual(semantics(bad), state)
            self.assertEqual(bad.reactions.R1026.notes, notes)


if __name__ == '__main__':
    unittest.main()
