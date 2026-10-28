"""R1726 opening curation: preservation, rejection, SBML roundtrip, the kept R1725 closure and the opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r1726_cpt_open import CURATION_PATH, apply_r1726_cpt_open
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_r2135_pool_20261027.xml")


class R1726CptOpenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        differs, kept = "R1726 opening precondition differs", "R1725 is not the closed reaction"
        expected = {"name": differs, "bounds": differs, "gpr": differs, "ec": differs, "stoichiometry": differs,
                    "species": differs, "compartment": differs,
                    "open_without_notes": "already open without this curation",
                    "notes": "Conflicting R1726 opening notes",
                    "r1725_open": kept, "r1725_gpr": kept, "r1725_name": kept, "r1725_stoichiometry": kept,
                    "r1725_missing": "R1726 or R1725 missing", "r1726_missing": "R1726 or R1725 missing",
                    "r1725_species": kept, "closed_with_notes": "closed but carries this curation's notes"}
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R1726
            if field == "name":
                reaction.name = "different transferase"
            elif field == "bounds":
                reaction.bounds = (0.0, 500.0)
            elif field == "open_without_notes":
                reaction.bounds = (0.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = "YALI1C15420g or YALI1E20159g"
            elif field == "ec":
                reaction.annotation["ec-code"] = "2.7.8.1"
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m377[C_em]"): -1})
            elif field == "species":
                # A lipid-unlump or microspecies change to the lumped PC would invalidate the review.
                bad.metabolites.get_by_id("m1632[C_em]").formula = "C44H84NO8P"
            elif field == "compartment":
                bad.metabolites.get_by_id("m1649[C_em]").compartment = "C_er"
            elif field == "notes":
                reaction.notes["r1726_cpt_open_curation"] = "conflicting review"
            elif field == "r1725_open":
                bad.reactions.R1725.bounds = (0.0, 1000.0)
            elif field == "r1725_gpr":
                bad.reactions.R1725.gene_reaction_rule = "YALI1E20159g"
            elif field == "r1725_name":
                bad.reactions.R1725.name = "ethanolaminephosphotransferase"
            elif field == "r1725_stoichiometry":
                # A rewritten but still closed R1725 is not the reaction this curation kept closed.
                bad.reactions.R1725.add_metabolites({bad.metabolites.get_by_id("m377[C_em]"): 1})
            elif field == "r1725_species":
                # m1633 (PE) is in R1725 only, so this hits the R1725 guard, not the R1726 one.
                bad.metabolites.get_by_id("m1633[C_em]").formula = "C41H78NO8P"
            elif field == "closed_with_notes":
                reaction.notes.update({"r1726_cpt_open_curation": self.spec["notes"]["r1726_cpt_open_curation"]})
            elif field == "r1725_missing":
                bad.remove_reactions([bad.reactions.R1725])
            else:
                bad.remove_reactions([bad.reactions.R1726])
            before = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r1726_cpt_open(bad)
            self.assertEqual((semantics(bad), annotations(bad)), before)

    def test_unexpected_scope_is_rejected(self):
        for path, value in ((("schema_version",), 2), (("curation_id",), "R1726-CPT-OPEN-20261028"),
                            (("status",), "proposed"), (("reaction_id",), "R1725"),
                            (("before_bounds",), [-1000, 1000]), (("after_bounds",), [-1000, 1000]),
                            (("kept_closed", "reaction_id"), "R1726"), (("kept_closed", "bounds"), [0, 1000])):
            spec = copy.deepcopy(self.spec)
            target = spec
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, "Unexpected R1726 opening curation scope"):
                apply_r1726_cpt_open(self.model.copy(), spec)

    def test_sets_only_the_r1726_bounds_and_notes_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected.reactions.R1726.bounds = (0.0, 1000.0)
        expected.reactions.R1726.notes.update(self.spec["notes"])
        self.assertEqual(apply_r1726_cpt_open(model)["status"], "applied")
        self.assertTrue((semantics(model), annotations(model)) ==
                        (semantics(expected), annotations(expected)), "Unexpected model change")
        self.assertEqual(model.reactions.R1726.bounds, (0.0, 1000.0))
        self.assertEqual(model.reactions.R1725.bounds, (0.0, 0.0))
        self.assertEqual(apply_r1726_cpt_open(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r1726.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue((semantics(loaded), annotations(loaded)) ==
                        (semantics(expected), annotations(expected)), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r1726_cpt_open(loaded)["status"], "already_correct")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r1726_cpt_open)
        self.assertTrue(parse_args(["--r1726-cpt-open"]).r1726_cpt_open)

    def test_flag_requires_the_r2135_acyl_coa_pool(self):
        with self.assertRaisesRegex(ValueError, "requires --r2135-acyl-coa-pool"):
            build_reference_chain(r1726_cpt_open=True)
        # The cofactor biomass-terms base alone is the one the first run stopped on (acyl leak).
        with self.assertRaisesRegex(ValueError, "requires --r2135-acyl-coa-pool"):
            build_reference_chain(r1726_cpt_open=True, lipoylation_module=True, cofactor_biomass_terms=True)


if __name__ == "__main__":
    unittest.main()
