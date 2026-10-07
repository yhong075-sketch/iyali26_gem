"""R539 catalytic GPR/EC curation, conflict rejection and KO propagation; no solve."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model
from scripts.gem_annotate.patches import apply_r539_gpr_assignment, R539_ASSIGNMENT_PATH
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics, annotations

ROOT = Path(__file__).resolve().parents[1]


class R539AssignmentTests(unittest.TestCase):
    def test_catalytic_scope_conflicts_and_roundtrip(self):
        model = read_sbml_model(str(ROOT / 'model_metadata_trna_r1931_forward.xml'))
        spec = json.loads(R539_ASSIGNMENT_PATH.read_text())
        original = semantics(model)
        old_annotations = annotations(model)
        old_genes = {g.id: (g.name, copy.deepcopy(g.annotation)) for g in model.genes}
        r = model.reactions.R539
        with model:
            model.genes.YALI1E22262g.knock_out()
            self.assertTrue(r.functional)
            self.assertEqual(r.bounds, (-1000, 1000))

        # Every failed precondition must leave all scientific fields untouched.
        for field in ('ec', 'gpr', 'species', 'identity', 'notes'):
            bad = model.copy()
            rx = bad.reactions.R539
            if field == 'ec':
                rx.annotation = {**rx.annotation, 'ec-code': '2.3.1.38'}
            elif field == 'gpr':
                rx.gene_reaction_rule = 'YALI1F38317g'
            elif field == 'species':
                next(iter(rx.metabolites)).charge = 17
            elif field == 'identity':
                gene = bad.genes.YALI1E22262g
                gene.annotation = {**gene.annotation, 'refseq': 'conflicting_accession'}
            else:
                rx.notes = {**rx.notes, 'gpr_experimental_confirmation': 'experimentally_verified'}
            state = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaises(ValueError):
                apply_r539_gpr_assignment(bad)
            self.assertEqual((semantics(bad), annotations(bad)), state)

        self.assertEqual(apply_r539_gpr_assignment(model)['status'], 'applied')
        self.assertEqual(apply_r539_gpr_assignment(model)['status'], 'already_correct')
        self.assertEqual(r.gene_reaction_rule, 'YALI1E22262g')
        self.assertEqual(r.annotation['ec-code'], '2.3.1.39')
        self.assertEqual(r.notes, {**old_annotations['R539'][2], **spec['notes']})
        expected_annotation = {**old_annotations['R539'][1], 'ec-code': '2.3.1.39'}
        self.assertEqual(r.annotation, expected_annotation)
        self.assertEqual({g.id: (g.name, g.annotation) for g in model.genes}, old_genes)
        with model:
            r.gene_reaction_rule = spec['before_gpr']
            self.assertEqual(semantics(model), original)
        self.assertEqual({k: v for k, v in annotations(model).items() if k != 'R539'},
                         {k: v for k, v in old_annotations.items() if k != 'R539'})
        bounds = {rx.id: rx.bounds for rx in model.reactions}
        with model:
            model.genes.YALI1E22262g.knock_out()
            self.assertFalse(r.functional)
            self.assertEqual({rx.id: rx.bounds for rx in model.reactions if rx.bounds != bounds[rx.id]},
                             {'R539': (0, 0)})
            self.assertTrue(model.reactions.R78.functional)
        for gid in ('YALI1F38317g', 'YALI1A20089g', 'YALI1C26939g',
                    'YALI1D18037g', 'YALI1D32594g', 'YALI1F37498g'):
            with model:
                model.genes.get_by_id(gid).knock_out()
                self.assertTrue(r.functional)
                self.assertEqual(r.bounds, (-1000, 1000))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'r539.xml'
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual(semantics(loaded), semantics(model))
        self.assertEqual(annotations(loaded), annotations(model))
        self.assertEqual(apply_r539_gpr_assignment(loaded)['status'], 'already_correct')


if __name__ == '__main__':
    unittest.main()
