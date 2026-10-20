"""R2041 forward-only curation: preservation, rejection, SBML roundtrip and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r2041_direction import CURATION_PATH, apply_r2041_direction
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = "E5_coq9_resp_pkg_r132_20261012.xml"


class R2041DirectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        for field in ("name", "bounds", "forward_without_notes", "gpr", "ec", "patch_stoichiometry",
                      "species", "compartment", "notes"):
            bad = self.model.copy()
            reaction = bad.reactions.R2041
            if field == "name":
                reaction.name = "different kinase"
            elif field == "bounds":
                reaction.bounds = (0.0, 0.0)
            elif field == "forward_without_notes":
                reaction.bounds = (0.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = bad.genes[0].id
            elif field == "ec":
                reaction.annotation["ec-code"] = "2.7.1.55"
            elif field == "patch_stoichiometry":
                # The builder patch's 1 H+ equation (before the selection) is not this step's input.
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): -1})
            elif field == "species":
                # An ATP microspecies migration changes the balanced proton count.
                atp = bad.metabolites.get_by_id("m141[C_cy]")
                atp.formula, atp.charge = "C10H12N5O13P3", -4
            elif field == "compartment":
                bad.metabolites.get_by_id("m454[C_cy]").compartment = "C_mi"
            else:
                reaction.notes["r2041_direction_curation"] = "conflicting review"
            before = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaises(ValueError):
                apply_r2041_direction(bad)
            self.assertEqual((semantics(bad), annotations(bad)), before)

    def test_unexpected_scope_is_rejected(self):
        spec = copy.deepcopy(self.spec)
        spec["after_bounds"] = [-1000, 0]
        with self.assertRaises(ValueError):
            apply_r2041_direction(self.model.copy(), spec)

    def test_sets_only_the_lower_bound_and_notes_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected.reactions.R2041.bounds = (0.0, 1000.0)
        expected.reactions.R2041.notes.update(self.spec["notes"])
        self.assertEqual(apply_r2041_direction(model)["status"], "applied")
        self.assertTrue((semantics(model), annotations(model)) ==
                        (semantics(expected), annotations(expected)), "Unexpected model change")
        self.assertEqual(model.reactions.R2041.bounds, (0.0, 1000.0))
        self.assertEqual(model.reactions.R2041.metabolites[model.metabolites.get_by_id("m10[C_cy]")], 2.0)
        self.assertEqual(apply_r2041_direction(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r2041.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue((semantics(loaded), annotations(loaded)) ==
                        (semantics(expected), annotations(expected)), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r2041_direction(loaded)["status"], "already_correct")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r2041_forward)
        self.assertTrue(parse_args(["--r2041-forward"]).r2041_forward)


if __name__ == "__main__":
    unittest.main()
