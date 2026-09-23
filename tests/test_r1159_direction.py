"""Check signed leak bounds, preservation, rejection and SBML export without LPs."""

import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.r1159_direction import CURATION_PATH, apply_r1159_direction
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

ROOT = Path(__file__).resolve().parents[1]


class R1159DirectionTests(unittest.TestCase):
    def test_direction_preservation_rejection_and_roundtrip(self):
        model = read_sbml_model(str(ROOT / "model_metadata_trna.xml"))
        spec = json.loads(CURATION_PATH.read_text())
        for field in ("name", "bounds", "gpr", "stoichiometry", "species", "notes"):
            bad = model.copy()
            reaction = bad.reactions.R1159
            if field == "name":
                reaction.name = "different transport"
            elif field == "bounds":
                reaction.bounds = (0, 0)
            elif field == "gpr":
                reaction.gene_reaction_rule = model.genes[0].id
            elif field == "stoichiometry":
                reaction.add_metabolites({next(iter(reaction.metabolites)): 1})
            elif field == "species":
                next(iter(reaction.metabolites)).compartment = "C_va"
            else:
                reaction.notes["r1159_direction_curation"] = "conflicting review"
            before = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaises(ValueError):
                apply_r1159_direction(bad)
            self.assertEqual((semantics(bad), annotations(bad)), before)

        expected = read_sbml_model(str(ROOT / "model_metadata_trna.xml"))
        expected.reactions.R1159.bounds = (-1000.0, 0.0)
        expected.reactions.R1159.notes.update(spec["notes"])
        self.assertEqual(apply_r1159_direction(model)["status"], "applied")
        self.assertTrue((semantics(model), annotations(model)) ==
                        (semantics(expected), annotations(expected)), "Unexpected model change")
        self.assertEqual(model.reactions.R1159.bounds, (-1000, 0))
        self.assertFalse(model.reactions.R1159.reversibility)
        self.assertEqual(apply_r1159_direction(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "leak.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue((semantics(loaded), annotations(loaded)) ==
                        (semantics(expected), annotations(expected)), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r1159_direction(loaded)["status"], "already_correct")


if __name__ == "__main__":
    unittest.main()
