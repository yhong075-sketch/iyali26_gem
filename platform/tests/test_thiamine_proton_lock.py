"""Thiamine proton lock: the gate honours the SBML-safe lock note, and the guarded step that sets it; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.microspecies import balance_protons_and_water
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.thiamine_proton_lock import CURATION_PATH, apply_thiamine_proton_lock
from tests.test_coq9_curation import annotations, semantics
from tests.test_gpi_chain_links import group_table

BASELINE = "E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_v2_20261022.xml"
TARGETS = ["R408", "R750", "R752"]


def state(model):
    return semantics(model), annotations(model), group_table(model)


def roundtrip(model):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "thiamine_proton_lock.xml"
        write_deterministic_sbml_model(model, path)
        return read_sbml_model(str(path))


class ProtonGateLockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))

    def gate(self, model):
        report = balance_protons_and_water(model, TARGETS)
        return sorted(c["reaction_id"] for c in report["changes"]), report["skipped_curated_lock_reaction_ids"]

    def test_gate_strips_the_h_without_a_lock(self):
        model = self.model.copy()
        self.assertEqual(self.gate(model), (TARGETS, []))
        self.assertNotIn(model.metabolites.get_by_id("m10[C_cy]"), model.reactions.R752.metabolites)

    def test_plain_lock_is_honoured_after_an_sbml_roundtrip(self):
        model = self.model.copy()
        for rid in TARGETS:
            model.reactions.get_by_id(rid).notes = {**model.reactions.get_by_id(rid).notes,
                                                    "lock_proton_water_stoichiometry": "true"}
        loaded = roundtrip(model)
        self.assertEqual(loaded.reactions.R752.notes["lock_proton_water_stoichiometry"], "true")
        stoichiometry = {rid: dict(loaded.reactions.get_by_id(rid).metabolites) for rid in TARGETS}
        self.assertEqual(self.gate(loaded), ([], TARGETS))
        self.assertEqual({rid: dict(loaded.reactions.get_by_id(rid).metabolites) for rid in TARGETS}, stoichiometry)

    def test_lock_value_is_case_and_space_insensitive_and_strict(self):
        for value, locked in ((" True ", True), ("TRUE", True), ("false", False), ("yes", False), ("", False)):
            model = self.model.copy()
            model.reactions.R408.notes = {**model.reactions.R408.notes, "lock_proton_water_stoichiometry": value}
            with self.subTest(value=value):
                changed, skipped = self.gate(model)
                self.assertEqual("R408" in skipped, locked)
                self.assertEqual("R408" in changed, not locked)

    def test_dict_lock_is_still_honoured_in_memory(self):
        model = self.model.copy()
        model.reactions.R750.notes = {**model.reactions.R750.notes, "curated_reaction_correction": {
            "lock_proton_water_stoichiometry": True}}
        self.assertEqual(self.gate(model), (["R408", "R752"], ["R750"]))


class ThiamineProtonLockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        expected = {"stoichiometry": "R752 precondition differs",
                    "gpr": "R408 precondition differs",
                    "residual": "R752 residual differs from the reviewed one H\\+",
                    "extra_note": "R750 notes outside this curation differ",
                    "conflicting_lock": "Conflicting R408 proton lock notes",
                    "lock_only": "R752 carries only part of this curation's notes",
                    "dict_note": "R752 notes outside this curation differ",
                    "bool_lock": "Conflicting R408 proton lock notes",
                    "partial": "Only part of the thiamine proton lock is present"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            if field == "stoichiometry":
                # The H+ removal the lock exists to prevent.
                bad.reactions.R752.subtract_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "gpr":
                bad.reactions.R408.gene_reaction_rule = "YALI1E25443g"
            elif field == "residual":
                # Without the thiamine phosphate chemistry the residual is not the reviewed one H+.
                bad.metabolites.get_by_id("m975[C_cy]").formula = None
            elif field == "extra_note":
                bad.reactions.R750.notes = {"unreviewed_note": "added elsewhere"}
            elif field == "conflicting_lock":
                bad.reactions.R408.notes = {**bad.reactions.R408.notes, "lock_proton_water_stoichiometry": "false"}
            elif field == "lock_only":
                bad.reactions.R752.notes = {**bad.reactions.R752.notes, "lock_proton_water_stoichiometry": "true"}
            elif field == "dict_note":
                # Non-string note values must fail with ValueError, not crash.
                bad.reactions.R752.notes = {**bad.reactions.R752.notes, "curated_reaction_correction": {
                    "lock_proton_water_stoichiometry": True}}
            elif field == "bool_lock":
                bad.reactions.R408.notes = {**bad.reactions.R408.notes, "lock_proton_water_stoichiometry": True}
            else:
                apply_thiamine_proton_lock(bad)
                bad.reactions.R750.notes = {}
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_thiamine_proton_lock(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_or_values_are_rejected(self):
        cases = (("status", "user_authorized", "scope/status"),
                 ("curation_id", "THIAMINE-PROTON-LOCK", "scope/status"),
                 ("schema_version", 1.0, "scope/status"),
                 ("reactions", {"R752": self.spec["reactions"]["R752"]}, "scope/status"),
                 ("lock_note", {"lock_proton_water_stoichiometry": "yes"}, "spec values differ"))
        for key, value, message in cases:
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, message):
                apply_thiamine_proton_lock(self.model.copy(), spec)
        spec = copy.deepcopy(self.spec)
        spec["reactions"]["R750"]["residual"] = {"H": 2.0, "charge": 2.0}
        with self.assertRaisesRegex(ValueError, "spec values differ"):
            apply_thiamine_proton_lock(self.model.copy(), spec)
        # An unsafe note text would make libSBML drop the whole notes block, lock included, on save.
        for key, value in (("rationale", "pH < 7.3 & neutral"), ("evidence", None)):
            spec = copy.deepcopy(self.spec)
            spec["reactions"]["R752"][key] = value
            model = self.model.copy()
            before = state(model)
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "note texts must be strings"):
                apply_thiamine_proton_lock(model, spec)
            self.assertEqual(state(model), before)

    def test_adds_only_the_notes_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        for rid, row in self.spec["reactions"].items():
            reaction = expected.reactions.get_by_id(rid)
            reaction.notes = {**reaction.notes, "lock_proton_water_stoichiometry": "true",
                              "thiamine_proton_lock_curation_id": self.spec["curation_id"],
                              "thiamine_proton_lock_status": self.spec["status"],
                              "thiamine_proton_lock_authorization": self.spec["authorization"],
                              "thiamine_proton_lock_evidence": row["evidence"],
                              "thiamine_proton_lock_evidence_tier": row["evidence_tier"],
                              "thiamine_proton_lock_rationale": row["rationale"],
                              "thiamine_proton_lock_previous_state": row["previous_state"]}
        result = apply_thiamine_proton_lock(model)
        self.assertEqual((result["status"], result["changed_reactions"]), ("applied", TARGETS))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual(apply_thiamine_proton_lock(model)["status"], "already_correct")
        loaded = roundtrip(model)
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_thiamine_proton_lock(loaded)["status"], "already_correct")
        report = balance_protons_and_water(loaded, TARGETS)
        self.assertEqual((report["changed_reactions"], report["skipped_curated_lock_reaction_ids"]), (0, TARGETS))

    def test_flag_requires_the_thiamine_phosphate_chemistry(self):
        self.assertFalse(parse_args([]).thiamine_proton_lock)
        self.assertTrue(parse_args(["--thiamine-proton-lock"]).thiamine_proton_lock)
        with self.assertRaisesRegex(ValueError, "defined on the thiamine phosphate chemistry candidate"):
            build_reference_chain(thiamine_proton_lock=True)


if __name__ == "__main__":
    unittest.main()
