"""Cofactor biomass terms: preservation, rejection, SBML roundtrip, exact edits and the flag requirement; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.cofactor_biomass_terms import CURATION_PATH, TERM_IDS, apply_cofactor_biomass_terms
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics
from tests.test_gpi_chain_links import group_table

BASELINE = "E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml"
THF, FAD, RIBOFLAVIN, NAD, NADP, COA = "m479[C_cy]", "m528[C_cy]", "m926[C_cy]", "m122[C_cy]", "m41[C_cy]", "m69[C_cy]"
COEFFICIENTS = {THF: 6.34e-5, FAD: 1e-5, RIBOFLAVIN: 9.9e-4, NAD: 2.65e-3, NADP: 5.7e-4, COA: 1.9e-4}


def state(model):
    return semantics(model), annotations(model), group_table(model)


def add_terms(model, mids):
    model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id(mid): -COEFFICIENTS[mid] for mid in mids})


class CofactorBiomassTermsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_spec_matches_the_module_constants_and_the_task(self):
        self.assertEqual(set(self.spec["terms"]), TERM_IDS)
        self.assertEqual({mid: row["coefficient_mmol_per_gDW"] for mid, row in self.spec["terms"].items()}, COEFFICIENTS)
        self.assertFalse(TERM_IDS & set(self.spec["accepted_before_state"]["biomass_stoichiometry"]))
        self.assertEqual(len(self.spec["accepted_before_state"]["biomass_stoichiometry"]), 55)

    def test_stale_inputs_fail_closed_without_edits(self):
        composition = "biomass_C composition differs from the reviewed baseline"
        expected = {"objective": "The objective is not biomass_C alone",
                    "missing_species": "m479\\[C_cy\\] missing",
                    "species_identity": "m528\\[C_cy\\] identity differs",
                    "species_context": "m69\\[C_cy\\] reaction context differs",
                    "composition": composition, "bounds": composition,
                    "extra_note": "biomass_C notes outside this curation differ",
                    "conflicting_note": "Conflicting biomass_C cofactor notes",
                    "coefficient": "Conflicting m926\\[C_cy\\] biomass coefficient",
                    "terms_without_notes": "already carries the cofactor terms without this curation",
                    "notes_without_terms": "carries this curation's notes without the terms",
                    "partial": "Only some of the cofactor biomass terms are present"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            biomass = bad.reactions.biomass_C
            if field == "objective":
                bad.objective = "R1710"
            elif field == "missing_species":
                bad.remove_metabolites([bad.metabolites.get_by_id(THF)])
            elif field == "species_identity":
                bad.metabolites.get_by_id(FAD).formula = "C27H31N9O15P2"
            elif field == "species_context":
                sink = Reaction("R_TEST_SINK", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([sink])
                sink.add_metabolites({bad.metabolites.get_by_id(COA): -1})
            elif field == "composition":
                # Another curation changed a biomass term after this review.
                biomass.add_metabolites({bad.metabolites.get_by_id("m468[C_mi]"): -1e-4})
            elif field == "bounds":
                biomass.bounds = (0.1, 1000.0)
            elif field == "extra_note":
                biomass.notes = {**biomass.notes, "unreviewed_note": "added elsewhere"}
            elif field == "conflicting_note":
                biomass.notes = {**biomass.notes, "cofactor_biomass_terms_status": "another review"}
            elif field == "coefficient":
                apply_cofactor_biomass_terms(bad)
                biomass.add_metabolites({bad.metabolites.get_by_id(RIBOFLAVIN): -1e-5})
            elif field == "terms_without_notes":
                add_terms(bad, TERM_IDS)
            elif field == "notes_without_terms":
                apply_cofactor_biomass_terms(bad)
                biomass.subtract_metabolites({bad.metabolites.get_by_id(mid): -c for mid, c in COEFFICIENTS.items()})
            else:
                apply_cofactor_biomass_terms(bad)
                biomass.subtract_metabolites({bad.metabolites.get_by_id(NAD): -COEFFICIENTS[NAD]})
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_cofactor_biomass_terms(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        before_state = self.spec["accepted_before_state"]
        for key, value in (("status", "user_authorized"), ("curation_id", "COFACTOR-BIOMASS-TERMS"),
                           ("biomass_reaction", "R1710"), ("requires_objective", "R1710"),
                           ("terms", {mid: row for mid, row in self.spec["terms"].items() if mid != COA}),
                           ("terms", {**self.spec["terms"], THF: {**self.spec["terms"][THF], "coefficient_mmol_per_gDW": 0.0}}),
                           ("accepted_before_state", None),
                           ("accepted_before_state", {k: v for k, v in before_state.items() if k != "id"}),
                           ("accepted_before_state", {**before_state, "id": "another_reviewed_state"}),
                           ("accepted_before_state", {**before_state, "biomass_stoichiometry": {**before_state["biomass_stoichiometry"], NAD: -2.65e-3}}),
                           ("accepted_before_state", {**before_state, "biomass_note_keys": None})):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            model = self.model.copy()
            before = state(model)
            with self.subTest(key=key, value=str(value)[:40]), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_cofactor_biomass_terms(model, spec)
            self.assertEqual(state(model), before, "A rejected spec must leave the model untouched")

    def test_adds_only_the_six_terms_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        biomass, terms = expected.reactions.biomass_C, self.spec["terms"]
        add_terms(expected, TERM_IDS)
        biomass.notes = {**biomass.notes, **{
            "cofactor_biomass_terms_curation_id": self.spec["curation_id"],
            "cofactor_biomass_terms_status": self.spec["status"],
            "cofactor_biomass_terms_authorization": self.spec["authorization"],
            "cofactor_biomass_terms_terms": "; ".join(f"{mid} {terms[mid]['coefficient_mmol_per_gDW']:g} mmol/gDW" for mid in sorted(terms)),
            "cofactor_biomass_terms_evidence": " | ".join(terms[mid]["source"] for mid in sorted(terms)),
            "cofactor_biomass_terms_evidence_tier": self.spec["evidence_tier"],
            "cofactor_biomass_terms_caveat": self.spec["organism_caveat"]}}
        result = apply_cofactor_biomass_terms(model)
        self.assertEqual((result["status"], result["changed_terms"], result["before_state"]),
                         ("applied", sorted(TERM_IDS), "a1237821_biomass_C"))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        coefficients = {m.id: c for m, c in model.reactions.biomass_C.metabolites.items()}
        self.assertEqual({mid: coefficients[mid] for mid in TERM_IDS}, {mid: -c for mid, c in COEFFICIENTS.items()})
        self.assertEqual(len(coefficients), 61)
        self.assertEqual(apply_cofactor_biomass_terms(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cofactor_biomass_terms.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_cofactor_biomass_terms(loaded)["status"], "already_correct")

    def test_flag_requires_the_lipoylation_module(self):
        self.assertFalse(parse_args([]).cofactor_biomass_terms)
        self.assertTrue(parse_args(["--cofactor-biomass-terms"]).cofactor_biomass_terms)
        with self.assertRaisesRegex(ValueError, "requires --lipoylation-module"):
            build_reference_chain(cofactor_biomass_terms=True)
        with self.assertRaisesRegex(ValueError, "requires --lipoylation-module"):
            build_reference_chain(cofactor_biomass_terms=True, acyl_acp_chemistry=True)


if __name__ == "__main__":
    unittest.main()
