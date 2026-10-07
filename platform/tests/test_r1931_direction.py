"""Direction curation, conflict rejection and SBML roundtrip; no optimization."""
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model
from scripts.gem_annotate.patches import apply_r1931_direction, R1931_DIRECTION_PATH
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics, annotations

ROOT = Path(__file__).resolve().parents[2]


from scripts.gem_annotate.model_layout import MODEL
class R1931DirectionTests(unittest.TestCase):
    def test_direction_preservation_conflicts_and_roundtrip(self):
        path = MODEL.candidate_file('model_metadata_trna_r153_merged.xml')
        model = read_sbml_model(str(path))
        for field in ('bounds', 'gpr', 'stoichiometry', 'species', 'ec', 'notes'):
            bad = model.copy()
            r = bad.reactions.R1931
            if field == 'bounds':
                r.upper_bound = 999
            elif field == 'gpr':
                r.gene_reaction_rule = ''
            elif field == 'stoichiometry':
                r.add_metabolites({next(iter(r.metabolites)): 1})
            elif field == 'species':
                next(iter(r.metabolites)).charge = 17
            elif field == 'ec':
                r.annotation['ec-code'] = '1.2.1.41'
            else:
                r.notes['r1931_direction_curation'] = 'conflict'
            before = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaises(ValueError):
                apply_r1931_direction(bad)
            self.assertEqual((semantics(bad), annotations(bad)), before)
        self.assertEqual(apply_r1931_direction(model)['status'], 'applied')
        expected = read_sbml_model(str(path))
        expected.reactions.R1931.bounds = (0.0, 1000.0)
        expected.reactions.R1931.notes.update(json.loads(R1931_DIRECTION_PATH.read_text())['notes'])
        self.assertEqual(semantics(model), semantics(expected))
        self.assertEqual(annotations(model), annotations(expected))
        self.assertFalse(model.reactions.R1931.reversibility)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'forward.xml'
            write_deterministic_sbml_model(model, output)
            loaded = read_sbml_model(str(output))
        self.assertEqual(semantics(loaded), semantics(expected))
        self.assertEqual(annotations(loaded), annotations(expected))
        self.assertEqual(apply_r1931_direction(loaded)['status'], 'already_correct')


if __name__ == '__main__':
    unittest.main()
