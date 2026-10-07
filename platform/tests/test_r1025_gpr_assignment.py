"""Check nuclear candidate evidence survives export and closes its KO bypass."""
import copy
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.patches import apply_r1025_gpr_assignment, apply_r1026_gpr_assignment
from scripts.gem_annotate.metabolites import normalize_all_annotations
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics

ROOT = Path(__file__).resolve().parents[1]


class R1025AssignmentTests(unittest.TestCase):
    def test_nuclear_assignment_roundtrip_and_precondition_conflicts(self):
        model = read_sbml_model(str(ROOT/'model_metadata_trna.xml'))
        apply_r1026_gpr_assignment(model)
        before = semantics(model)
        normalized = model.copy()
        normalize_all_annotations(normalized)
        self.assertEqual(apply_r1025_gpr_assignment(normalized)['status'], 'applied')
        self.assertEqual(apply_r1025_gpr_assignment(model)['status'], 'applied')
        rx = model.reactions.R1025
        self.assertEqual(rx.gene_reaction_rule, 'YALI1F28274g')
        self.assertEqual(rx.notes['gpr_evidence_basis'], 'AlphaFold prediction')
        self.assertEqual(rx.notes['gpr_evidence_label'], 'AlphaFold 预测支持，仍需实验验证')
        self.assertEqual(rx.notes['gpr_experimental_confirmation'], 'required_not_confirmed')
        original_bounds = {r.id:r.bounds for r in model.reactions}
        with model:
            model.genes.YALI1F28274g.knock_out()
            self.assertEqual({r.id:r.bounds for r in model.reactions if r.bounds!=original_bounds[r.id]},
                             {'R1025':(0,0), 'R1026':(0,0), 'R2202':(0,0)})
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'r1025.xml'
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual(loaded.reactions.R1025.notes, rx.notes)
        self.assertEqual(apply_r1025_gpr_assignment(loaded)['status'], 'already_correct')
        loaded.reactions.R1025.gene_reaction_rule = ''
        self.assertEqual(semantics(loaded), before)

        for conflict in ('gpr', 'compartment', 'bound', 'gene_identity', 'evidence'):
            bad = model.copy()
            if conflict == 'gpr':
                bad.reactions.R1025.gene_reaction_rule = 'YALI1A08512g'
            elif conflict == 'compartment':
                bad.metabolites.get_by_id('m1175[C_nu]').compartment = 'C_cy'
            elif conflict == 'bound':
                bad.reactions.R1025.lower_bound = 0
            elif conflict == 'gene_identity':
                bad.genes.YALI1F28274g.annotation['uniprot'] = 'conflicting_accession'
            else:
                bad.reactions.R1025.notes['gpr_experimental_confirmation'] = 'experimentally_verified'
            state, notes = semantics(bad), copy.deepcopy(bad.reactions.R1025.notes)
            with self.subTest(conflict=conflict), self.assertRaises(ValueError):
                apply_r1025_gpr_assignment(bad)
            self.assertEqual(semantics(bad), state)
            self.assertEqual(bad.reactions.R1025.notes, notes)


if __name__ == '__main__':
    unittest.main()
