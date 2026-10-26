"""Heme a + iron biomass terms: preservation, rejection, SBML roundtrip, exact edits and the flag requirement; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.heme_iron_biomass import CURATION_PATH, apply_heme_iron_biomass
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.thiamine_phosphate_chemistry import CURATION_PATH as CHEMISTRY_PATH, apply_thiamine_phosphate_chemistry
from scripts.gem_annotate.thiamine_proton_lock import CURATION_PATH as LOCK_PATH, apply_thiamine_proton_lock
from scripts.gem_annotate.tpp_requirement import CURATION_PATH as TPP_CURATION_PATH, apply_tpp_requirement
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
                           ("requires_reactions", {}), ("requires_objective", "R1710"),
                           # The base must stay the before-state of tpp_requirement.json.
                           ("biomass_bounds_before", [0.0, 999.0]),
                           ("biomass_note_keys_before", self.spec["biomass_note_keys_before"][1:]),
                           ("biomass_stoichiometry_before", {**self.spec["biomass_stoichiometry_before"],
                                                             "m1000[C_cy]": -0.05})):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_heme_iron_biomass(self.model.copy(), spec)
        # A changed TPP curation fails every heme build closed, on the base state too.
        tpp_spec = json.loads(TPP_CURATION_PATH.read_text())
        tpp_spec["status"] = "rejected"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tpp_requirement.json"
            path.write_text(json.dumps(tpp_spec))
            with patch("scripts.gem_annotate.heme_iron_biomass.TPP_CURATION_PATH", path), \
                    self.assertRaisesRegex(ValueError, "TPP requirement scope/status"):
                apply_heme_iron_biomass(self.model.copy())
        states = self.spec["accepted_before_states"]
        for label, value in (("missing_tpp_state", states[:1]),
                             ("base_with_additions", [{**states[0], "added_stoichiometry": {"m975[C_cy]": -1.2e-6}}, states[1]]),
                             ("state_with_heme_term", [states[0], {**states[1], "added_stoichiometry": {HEME_A: -1e-6}}]),
                             # The TPP state must stay the after-state of tpp_requirement.json.
                             ("tpp_coefficient_copy", [states[0], {**states[1], "added_stoichiometry": {"m975[C_cy]": -2.4e-6}}]),
                             ("tpp_note_keys_copy", [states[0], {**states[1], "added_note_keys": states[1]["added_note_keys"][1:]}]),
                             ("repeats_base_key", [states[0], {**states[1], "added_stoichiometry": {
                                 "m975[C_cy]": -1.2e-6, "m1000[C_cy]": -0.05}}])):
            spec = copy.deepcopy(self.spec)
            spec["accepted_before_states"] = value
            with self.subTest(states=label), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_heme_iron_biomass(self.model.copy(), spec)

    def test_applies_after_the_full_reviewed_tpp_line_only(self):
        tpp_state = self.spec["accepted_before_states"][1]
        line = self.model.copy()
        apply_tpp_requirement(line)
        apply_thiamine_phosphate_chemistry(line)
        apply_thiamine_proton_lock(line)
        model = line.copy()
        result = apply_heme_iron_biomass(model)
        self.assertEqual((result["status"], result["before_state"]), ("applied", "after_TPP-REQUIREMENT-20261022"))
        coefficients = {m.id: c for m, c in model.reactions.biomass_C.metabolites.items()}
        self.assertEqual((coefficients[HEME_A], coefficients[IRON], coefficients["m975[C_cy]"]),
                         (-1e-6, -3.04e-5, tpp_state["added_stoichiometry"]["m975[C_cy]"]))
        self.assertEqual(apply_heme_iron_biomass(model)["status"], "already_correct")
        # The reverse order fails closed: the TPP guard rejects a biomass_C that carries the heme terms.
        reverse = self.model.copy()
        apply_heme_iron_biomass(reverse)
        before = state(reverse)
        with self.assertRaisesRegex(ValueError, "biomass_C precondition differs"):
            apply_tpp_requirement(reverse)
        self.assertEqual(state(reverse), before)
        partial = "combine only with the full reviewed TPP line"
        for label in ("tpp_only", "no_proton_lock", "hand_chemistry_bare_lock", "partial_proton_lock",
                      "r1342_reopened", "r1342_tpp_status", "r1342_chemistry_notes_removed",
                      "lock_curation_notes_removed", "bool_lock_note", "no_r750", "m975_chemistry_notes_removed",
                      "m643_chemistry_notes_removed", "m976_formula", "no_m643", "other_tpp_curation",
                      "rejected_tpp_status", "bool_tpp_note", "tpp_coefficient", "tpp_without_notes"):
            if label == "tpp_only":
                bad = self.model.copy()
                apply_tpp_requirement(bad)
                message = partial
            elif label == "no_proton_lock":
                bad = self.model.copy()
                apply_tpp_requirement(bad)
                apply_thiamine_phosphate_chemistry(bad)
                message = partial
            elif label == "hand_chemistry_bare_lock":
                # The reviewed TPP chemistry and the lock note, set by hand without either curation.
                bad = self.model.copy()
                apply_tpp_requirement(bad)
                tpp = bad.metabolites.get_by_id("m975[C_cy]")
                tpp.formula, tpp.charge = "C12H19N4O7P2S", 1
                for rid in ("R752", "R408", "R750"):
                    reaction = bad.reactions.get_by_id(rid)
                    reaction.notes = {**reaction.notes, "lock_proton_water_stoichiometry": "true"}
                message = f"{partial}.*R752.*m975\\[C_cy\\]"
            else:
                bad = line.copy()
                biomass = bad.reactions.biomass_C
                if label == "r1342_reopened":
                    bad.reactions.R1342.bounds = (0.0, 1000.0)
                    message = f"{partial}.*; R1342 not in"
                elif label == "r1342_tpp_status":
                    bad.reactions.R1342.notes = {**bad.reactions.R1342.notes, "tpp_requirement_status": "rejected"}
                    message = f"{partial}.*; R1342 not in"
                elif label == "r1342_chemistry_notes_removed":
                    bad.reactions.R1342.notes = {k: v for k, v in bad.reactions.R1342.notes.items()
                                                 if not k.startswith("thiamine_phosphate_chemistry_")}
                    message = f"{partial}.*; R1342 not in"
                elif label == "bool_lock_note":
                    # Non-string note values must fail with ValueError, not crash.
                    bad.reactions.R408.notes = {**bad.reactions.R408.notes, "lock_proton_water_stoichiometry": True}
                    message = f"{partial}.*; R408 not in"
                elif label == "no_r750":
                    bad.remove_reactions([bad.reactions.R750])
                    message = f"{partial}.*; R750 not in"
                elif label == "m976_formula":
                    # The chemistry notes are present, but the formula is not the reviewed one.
                    bad.metabolites.get_by_id("m976[C_cy]").formula = "C12H19N4O10P3S"
                    message = f"{partial}.*; m976\\[C_cy\\] not in"
                elif label == "no_m643":
                    bad.remove_metabolites([bad.metabolites.get_by_id("m643[C_cy]")])
                    message = f"{partial}.*; m643\\[C_cy\\] not in"
                elif label == "bool_tpp_note":
                    biomass.notes = {**biomass.notes, "tpp_requirement_status": True}
                    message = "notes outside this curation differ"
                elif label == "lock_curation_notes_removed":
                    for rid in ("R752", "R408", "R750"):
                        reaction = bad.reactions.get_by_id(rid)
                        reaction.notes = {k: v for k, v in reaction.notes.items()
                                          if not k.startswith("thiamine_proton_lock_")}
                    message = f"{partial}.*; R752, R408, R750 not in"
                elif label in ("m975_chemistry_notes_removed", "m643_chemistry_notes_removed"):
                    species = bad.metabolites.get_by_id(f"{label[:4]}[C_cy]")
                    species.notes = {k: v for k, v in species.notes.items()
                                     if not k.startswith("thiamine_phosphate_chemistry_")}
                    message = f"{partial}.*; {label[:4]}\\[C_cy\\] not in"
                elif label == "partial_proton_lock":
                    bad.reactions.R750.notes = {k: v for k, v in bad.reactions.R750.notes.items()
                                                if k != "lock_proton_water_stoichiometry"}
                    message = partial
                elif label == "other_tpp_curation":
                    # The TPP note keys are present, but they record another curation.
                    biomass.notes = {**biomass.notes, "tpp_requirement_curation_id": "SOME-OTHER-TPP-CURATION"}
                    message = "notes outside this curation differ"
                elif label == "rejected_tpp_status":
                    biomass.notes = {**biomass.notes, "tpp_requirement_status": "rejected"}
                    message = "notes outside this curation differ"
                elif label == "tpp_coefficient":
                    # A TPP term other than the reviewed one.
                    biomass.add_metabolites({bad.metabolites.get_by_id("m975[C_cy]"): -1.2e-6})
                    message = "composition differs from the reviewed baseline"
                else:
                    biomass.notes = {k: v for k, v in biomass.notes.items() if not k.startswith("tpp_requirement_")}
                    message = "notes outside this curation differ"
            before = state(bad)
            with self.subTest(label=label), self.assertRaisesRegex(ValueError, message):
                apply_heme_iron_biomass(bad)
            self.assertEqual(state(bad), before)
        # The chemistry and lock files are read only in the TPP state: a changed one fails that state closed and
        # leaves base-state builds alone.
        for name, source, message in (("CHEMISTRY_PATH", CHEMISTRY_PATH, "thiamine phosphate chemistry scope/status"),
                                      ("LOCK_PATH", LOCK_PATH, "thiamine proton lock scope/status")):
            changed = json.loads(source.read_text())
            changed["status"] = "rejected"
            with self.subTest(changed=name), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / source.name
                path.write_text(json.dumps(changed))
                with patch(f"scripts.gem_annotate.heme_iron_biomass.{name}", path):
                    with self.assertRaisesRegex(ValueError, message):
                        apply_heme_iron_biomass(line.copy())
                    self.assertEqual(apply_heme_iron_biomass(self.model.copy())["status"], "applied")

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
        # The TPP requirement alone (without its chemistry and proton lock) was never reviewed with these terms.
        line = dict(energy_candidate="E5", no_solve=True, allow_network=False, coq9_biomass_dilution=True,
                    gpi_chain_links=True, gpi_coproducts=True, tpp_requirement=True, heme_iron_biomass=True)
        for extra in ({}, {"thiamine_phosphate_chemistry": True}):
            with self.subTest(extra=extra), self.assertRaisesRegex(ValueError, "combine only with the full reviewed TPP line"):
                build_reference_chain(**line, **extra)


if __name__ == "__main__":
    unittest.main()
