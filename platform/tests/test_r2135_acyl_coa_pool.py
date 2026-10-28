"""R2135 acyl-CoA pool remap and R2145 rescale: preservation, rejection, SBML roundtrip and the opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r2135_acyl_coa_pool import CURATION_PATH, apply_r2135_acyl_coa_pool
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_20261027.xml")


class R2135AcylCoaPoolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def expected(self):
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        for rid, row in self.spec["reactions"].items():
            reaction = model.reactions.get_by_id(rid)
            for met in list(reaction.metabolites):
                reaction.subtract_metabolites({met: reaction.metabolites[met]})
            reaction.add_metabolites({model.metabolites.get_by_id(m): c for m, c in row["after_stoichiometry"].items()})
            reaction.notes.update(row["notes"])
        return model

    def test_stale_inputs_fail_closed_without_edits(self):
        species = "A species of the R2135 / R2145 curation differs"
        expected = {"r2135_name": "R2135 precondition differs", "r2145_gpr": "R2145 precondition differs",
                    "r2145_bounds": "R2145 precondition differs", "r2135_bounds": "R2135 precondition differs",
                    "r2135_stoichiometry": "R2135 stoichiometry is neither", "r2145_stoichiometry": "R2145 stoichiometry is neither",
                    "after_without_notes": "R2135 stoichiometry is neither",
                    "half_applied": "Only one of R2135 / R2145 carries this curation",
                    "notes": "Conflicting R2135 curation notes",
                    "acetyl_coa_formula": species, "pool_identity": species, "pool_compartment": species,
                    "phosphate_charge": species, "er_acyl_coa_identity": species, "er_acyl_coa_new_consumer": species,
                    "r2145_conflicting_notes": "Conflicting R2145 curation notes",
                    "reverse_half_applied": "Only one of R2135 / R2145 carries this curation",
                    "partial_notes": "R2135 stoichiometry is neither",
                    "notes_without_equation": "R2145 carries this curation's notes without its equation",
                    "missing_species": "or one of their species missing", "missing_reaction": "or one of their species missing",
                    "missing_r2135": "or one of their species missing"}
        for field, message in expected.items():
            bad = self.model.copy()
            r2135, r2145 = bad.reactions.R2135, bad.reactions.R2145
            if field == "r2135_name":
                r2135.name = "acetyl-CoA transport"
            elif field == "r2145_gpr":
                r2145.gene_reaction_rule = "YALI1D35593g or YALI1E20049g"
            elif field == "r2145_bounds":
                r2145.bounds = (-1000.0, 1000.0)
            elif field == "r2135_bounds":
                r2135.bounds = (0.0, 1000.0)
            elif field == "r2135_stoichiometry":
                r2135.add_metabolites({bad.metabolites.get_by_id("m567[C_em]"): -1.0})
            elif field == "r2145_stoichiometry":
                r2145.add_metabolites({bad.metabolites.get_by_id("m1993[C_cy]"): -0.5})
            elif field == "after_without_notes":
                # The reviewed after-equation without this curation's notes came from somewhere else.
                r2135.subtract_metabolites({bad.metabolites.get_by_id("m68[C_cy]"): -1.0})
                r2135.add_metabolites({bad.metabolites.get_by_id("m567[C_em]"): -1.0})
            elif field == "half_applied":
                row = self.spec["reactions"]["R2145"]
                r2145.add_metabolites({bad.metabolites.get_by_id(m): c for m, c in row["after_stoichiometry"].items()},
                                      combine=False)
                r2145.notes.update(row["notes"])
            elif field == "notes":
                r2135.notes["r2135_acyl_coa_pool_curation"] = "conflicting review"
            elif field == "acetyl_coa_formula":
                bad.metabolites.get_by_id("m68[C_cy]").formula = "C23H34N7O17P3S"
            elif field == "pool_identity":
                # Replace, not mutate: Model.copy() shares each metabolite's annotation dict with the fixture.
                met = bad.metabolites.get_by_id("m567[C_em]")
                met.annotation = {**met.annotation, "metanetx.chemical": "MNXM1104266"}
            elif field == "pool_compartment":
                bad.metabolites.get_by_id("m1984[C_er]").compartment = "C_em"
            elif field == "phosphate_charge":
                bad.metabolites.get_by_id("m35[C_cy]").charge = -2
            elif field == "er_acyl_coa_identity":
                met = bad.metabolites.get_by_id("m1984[C_er]")
                met.annotation = {**met.annotation, "metanetx.chemical": "MNXM1104266"}
            elif field == "er_acyl_coa_new_consumer":
                consumer = Reaction("R_new_acyl_coa_consumer", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([consumer])
                consumer.add_metabolites({bad.metabolites.get_by_id("m1984[C_er]"): -1.0})
            elif field == "r2145_conflicting_notes":
                r2145.notes["r2145_rescale_curation"] = "conflicting review"
            elif field in ("reverse_half_applied", "partial_notes"):
                row = self.spec["reactions"]["R2135"]
                r2135.subtract_metabolites({bad.metabolites.get_by_id("m68[C_cy]"): -1.0})
                r2135.add_metabolites({bad.metabolites.get_by_id("m567[C_em]"): -1.0})
                notes = row["notes"] if field == "reverse_half_applied" else dict(list(row["notes"].items())[:2])
                r2135.notes.update(notes)
            elif field == "notes_without_equation":
                r2145.notes.update(self.spec["reactions"]["R2145"]["notes"])
            elif field == "missing_species":
                bad.remove_metabolites([bad.metabolites.get_by_id("m567[C_em]")])
            elif field == "missing_r2135":
                bad.remove_reactions([r2135])
            else:
                bad.remove_reactions([r2145])
            before = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r2135_acyl_coa_pool(bad)
            self.assertEqual((semantics(bad), annotations(bad)), before)

    def test_unexpected_scope_is_rejected(self):
        swapped = copy.deepcopy(self.spec["reactions"])
        swapped["R2135"]["after_stoichiometry"] = {"m68[C_cy]": -1.0, "m567[C_em]": 1.0}
        unknown = copy.deepcopy(self.spec["reactions"])
        unknown["R2145"]["after_stoichiometry"]["m9999[C_cy]"] = 1.0
        species = {k: v for k, v in self.spec["species"].items() if k != "m35[C_cy]"}
        for key, value in (("schema_version", 2), ("curation_id", "R2135-ACYL-COA-POOL-20261028"),
                           ("status", "proposed"), ("reactions", {"R2135": self.spec["reactions"]["R2135"]}),
                           ("reactions", swapped), ("reactions", unknown), ("species", species)):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "Unexpected R2135 acyl-CoA pool curation scope"):
                apply_r2135_acyl_coa_pool(self.model.copy(), spec)

    def test_changes_only_the_two_reactions_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = self.expected()
        self.assertEqual(apply_r2135_acyl_coa_pool(model)["status"], "applied")
        self.assertTrue((semantics(model), annotations(model)) ==
                        (semantics(expected), annotations(expected)), "Unexpected model change")
        self.assertEqual(model.reactions.R2135.build_reaction_string(), "m567[C_em] <=> m1984[C_er]")
        self.assertEqual(model.reactions.R2145.build_reaction_string(),
                         "m1993[C_cy] + m32[C_cy] --> m1994[C_cy] + m35[C_cy]")
        self.assertEqual(apply_r2135_acyl_coa_pool(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r2135.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue((semantics(loaded), annotations(loaded)) ==
                        (semantics(expected), annotations(expected)), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r2135_acyl_coa_pool(loaded)["status"], "already_correct")

    def test_flag_is_opt_in_and_requires_the_cofactor_biomass_terms(self):
        self.assertFalse(parse_args([]).r2135_acyl_coa_pool)
        self.assertTrue(parse_args(["--r2135-acyl-coa-pool"]).r2135_acyl_coa_pool)
        with self.assertRaisesRegex(ValueError, "requires --cofactor-biomass-terms"):
            build_reference_chain(r2135_acyl_coa_pool=True)
        with self.assertRaisesRegex(ValueError, "requires --cofactor-biomass-terms"):
            build_reference_chain(r2135_acyl_coa_pool=True, lipoylation_module=True)


if __name__ == "__main__":
    unittest.main()
