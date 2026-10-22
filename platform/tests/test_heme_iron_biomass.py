"""Heme a + iron biomass terms: preservation, rejection, SBML roundtrip, exact edits and the flag requirement; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.heme_iron_biomass import CURATION_PATH, apply_heme_iron_biomass
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics
from tests.test_gpi_chain_links import group_table

BASELINE = "E5_coq9_resp_pkg_r132_combined_gpi_coproducts_20261021.xml"
HEME_A, IRON = "m1635[C_cy]", "m948[C_cy]"


def state(model):
    return semantics(model), annotations(model), group_table(model)


class HemeIronBiomassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        composition = "biomass_C composition differs from the reviewed baseline"
        expected = {"no_export": "R1709 missing",
                    "no_iron_exchange": "R1189 missing",
                    "altered_import": "R2193 differs from the reviewed supply route",
                    "closed_mito_import": "R863 differs from the reviewed supply route",
                    "heme_a_synthase": "R383 differs from the reviewed supply route",
                    "objective": "The objective is not biomass_C alone",
                    "species_identity": "m1635\\[C_cy\\] identity differs",
                    "iron_context": "m948\\[C_cy\\] reaction context differs",
                    "composition": composition, "bounds": composition,
                    "extra_note": "biomass_C notes outside this curation differ",
                    "conflicting_note": "Conflicting biomass_C heme/iron notes",
                    "coefficient": "Conflicting m948\\[C_cy\\] biomass coefficient",
                    "terms_without_notes": "already carries the heme/iron terms without this curation",
                    "notes_without_terms": "carries this curation's notes without the terms",
                    "partial": "Only one heme/iron biomass term is present"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            biomass = bad.reactions.biomass_C
            heme_a, iron = bad.metabolites.get_by_id(HEME_A), bad.metabolites.get_by_id(IRON)
            if field == "no_export":
                bad.remove_reactions([bad.reactions.R1709])
            elif field == "no_iron_exchange":
                bad.remove_reactions([bad.reactions.R1189])
            elif field == "altered_import":
                # Present but closed: iron would have no route into the cytosol.
                bad.reactions.R2193.bounds = (0.0, 0.0)
            elif field == "closed_mito_import":
                # Ferrochelatase would get no iron in the mitochondrion.
                bad.reactions.R863.bounds = (0.0, 0.0)
            elif field == "heme_a_synthase":
                bad.reactions.R383.add_metabolites({bad.metabolites.get_by_id("m28[C_mi]"): 1})
            elif field == "objective":
                bad.objective = "R1710"
            elif field == "species_identity":
                heme_a.formula = "C49H58FeN4O6"
            elif field == "iron_context":
                sink = Reaction("R_TEST_SINK", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([sink])
                sink.add_metabolites({iron: -1})
            elif field == "composition":
                # Another curation changed a biomass term after this review.
                biomass.add_metabolites({bad.metabolites.get_by_id("m468[C_mi]"): -1e-4})
            elif field == "bounds":
                biomass.bounds = (0.1, 1000.0)
            elif field == "extra_note":
                biomass.notes = {**biomass.notes, "unreviewed_note": "added elsewhere"}
            elif field == "conflicting_note":
                biomass.notes = {**biomass.notes, "heme_iron_biomass_status": "another review"}
            elif field == "coefficient":
                apply_heme_iron_biomass(bad)
                biomass.add_metabolites({iron: -1e-5})
            elif field == "terms_without_notes":
                biomass.add_metabolites({heme_a: -1e-6, iron: -3.04e-5})
            elif field == "notes_without_terms":
                apply_heme_iron_biomass(bad)
                biomass.subtract_metabolites({heme_a: -1e-6, iron: -3.04e-5})
            else:
                apply_heme_iron_biomass(bad)
                biomass.subtract_metabolites({heme_a: -1e-6})
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_heme_iron_biomass(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for key, value in (("status", "user_authorized"), ("curation_id", "HEME-IRON-BIOMASS"),
                           ("biomass_reaction", "R1710"), ("terms", {HEME_A: self.spec["terms"][HEME_A]}),
                           ("requires_reactions", {}), ("requires_objective", "R1710")):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_heme_iron_biomass(self.model.copy(), spec)

    def test_adds_only_the_two_terms_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        biomass, terms = expected.reactions.biomass_C, self.spec["terms"]
        biomass.add_metabolites({expected.metabolites.get_by_id(HEME_A): -1e-6,
                                 expected.metabolites.get_by_id(IRON): -3.04e-5})
        biomass.notes = {**biomass.notes, **{
            "heme_iron_biomass_curation_id": self.spec["curation_id"],
            "heme_iron_biomass_status": self.spec["status"],
            "heme_iron_biomass_authorization": self.spec["authorization"],
            "heme_iron_biomass_terms": "m1635[C_cy] 1e-06 mmol/gDW; m948[C_cy] 3.04e-05 mmol/gDW",
            "heme_iron_biomass_evidence": f"{terms[HEME_A]['source']} | {terms[IRON]['source']}",
            "heme_iron_biomass_evidence_tier": self.spec["evidence_tier"],
            "heme_iron_biomass_caveat": self.spec["organism_caveat"]}}
        result = apply_heme_iron_biomass(model)
        self.assertEqual((result["status"], result["changed_terms"]), ("applied", [HEME_A, IRON]))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        coefficients = {m.id: c for m, c in model.reactions.biomass_C.metabolites.items()}
        self.assertEqual((coefficients[HEME_A], coefficients[IRON]), (-1e-6, -3.04e-5))
        self.assertEqual(apply_heme_iron_biomass(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "heme_iron_biomass.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_heme_iron_biomass(loaded)["status"], "already_correct")

    def test_flag_requires_the_coq9_biomass_stage(self):
        self.assertFalse(parse_args([]).heme_iron_biomass)
        self.assertTrue(parse_args(["--heme-iron-biomass"]).heme_iron_biomass)
        with self.assertRaisesRegex(ValueError, "defined on the CoQ9 biomass candidate"):
            build_reference_chain(heme_iron_biomass=True)


if __name__ == "__main__":
    unittest.main()
