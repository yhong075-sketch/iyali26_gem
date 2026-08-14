"""Regression tests for the machine-readable fatty acyl-CoA handoff."""

from __future__ import annotations

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model, write_sbml_model

import scripts.gem_annotate.patches as coa_patches
from scripts.fatty_acyl_coa_handoff import (
    HANDOFF_PATH,
    HandoffError,
    audit_handoff,
    load_handoff,
    _validate_output_path,
    validate_handoff_inputs,
)


REPOSITORY = Path(__file__).resolve().parents[1]
CURRENT_MAIN_MODEL = REPOSITORY.parent / "iyali26_gem" / "model.xml"
MODEL_PATH = CURRENT_MAIN_MODEL
CURRENT_MAIN_REPORT = REPOSITORY / "artifacts" / "fatty_acyl_coa_handoff_20260806.json"


class FattyAcylCoAHandoffTests(unittest.TestCase):
    def test_handoff_matches_pipeline_inputs_and_complete_inventory(self) -> None:
        handoff = load_handoff()
        validate_handoff_inputs(handoff)
        self.assertEqual(len(handoff["concrete_acyl_coa_groups"]), 7)
        self.assertEqual(
            len([copy for group in handoff["concrete_acyl_coa_groups"] for copy in group["copies"]]),
            36,
        )
        self.assertEqual(len(handoff["generic_acyl_coa_pools"]), 3)
        self.assertFalse(handoff["activation"]["ready_for_activation"])
        review = handoff["authoritative_inputs"]["r1521_current_snapshot_handoff"]
        self.assertEqual(review["activation_state"], "blocked")
        self.assertFalse(review["ready_for_activation"])
        r1521_blocker = next(
            row for row in handoff["remaining_activation_blockers"]
            if row["reaction_id"] == "R1521"
        )
        self.assertEqual(
            r1521_blocker["status"],
            "connected_component_migration_blocked",
        )

    def test_generic_pool_cannot_become_a_chemical_tuple(self) -> None:
        changed = copy.deepcopy(load_handoff())
        changed["generic_acyl_coa_pools"][0]["formula"] = "C1"
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temp_dir:
            path = Path(temp_dir) / "handoff.json"
            path.write_text(json.dumps(changed), encoding="utf-8")
            with self.assertRaises(HandoffError):
                load_handoff(path)

    def test_malformed_activation_is_a_handoff_error(self) -> None:
        changed = copy.deepcopy(load_handoff())
        changed["activation"] = None
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temp_dir:
            path = Path(temp_dir) / "handoff.json"
            path.write_text(json.dumps(changed), encoding="utf-8")
            with self.assertRaises(HandoffError):
                load_handoff(path)

    def test_report_output_cannot_alias_an_authoritative_input(self) -> None:
        handoff = load_handoff()
        protected = [
            MODEL_PATH,
            CURRENT_MAIN_MODEL,
            HANDOFF_PATH,
            REPOSITORY / handoff["authoritative_inputs"]["coa_protonation_curation"]["path"],
            REPOSITORY / handoff["authoritative_inputs"]["lipid_moiety_ledger_spec"]["path"],
            REPOSITORY / handoff["authoritative_inputs"]["r1521_current_snapshot_handoff"]["path"],
            Path(handoff["authoritative_inputs"]["human_review_table"]["path"]),
        ]
        for path in protected:
            with self.subTest(path=path):
                with self.assertRaises(HandoffError):
                    _validate_output_path(path, MODEL_PATH)

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_report_is_read_only_and_blocked_on_its_declared_source_snapshot(self) -> None:
        report = audit_handoff(MODEL_PATH)
        self.assertTrue(report["current_model_matches_declared_source"])
        self.assertFalse(report["ready_for_activation"])
        self.assertEqual(
            [row["reaction_id"] for row in report["generic_pool_policy"]],
            ["xPOOL_AC_EM", "xPOOL_AC_LP", "xPOOL_AC_MM"],
        )
        self.assertTrue(all(row["formula_is_absent"] for row in report["generic_pool_policy"]))
        self.assertTrue(all(row["handoff_tuple_is_null"] for row in report["generic_pool_policy"]))
        self.assertTrue(all(row["charge_is_nonchemical_placeholder"] for row in report["generic_pool_policy"]))
        self.assertTrue(all(row["chebi_is_absent"] for row in report["generic_pool_policy"]))
        self.assertTrue(all(row["no_complete_chemical_tuple"] for row in report["generic_pool_policy"]))
        self.assertTrue(all(row["product_is_the_only_positive_member"] for row in report["generic_pool_policy"]))
        self.assertTrue(all(row["inventory_linked"] for row in report["generic_pool_policy"]))
        self.assertTrue(all(len(row["actual_member_ids"]) == 7 for row in report["generic_pool_policy"]))

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_candidate_survives_sbml_roundtrip_without_relaxing_blocked_state(self) -> None:
        model = read_sbml_model(str(MODEL_PATH))
        curation = coa_patches.load_coa_protonation_curation()
        coa_patches._apply_coa_protonation_curation(model, curation)
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temp_dir:
            path = Path(temp_dir) / "candidate.xml"
            write_sbml_model(model, str(path))
            loaded = read_sbml_model(str(path))
        validation = coa_patches.validate_coa_protonation_curation(loaded, curation)
        self.assertTrue(validation["valid"], validation["errors"])
        self.assertEqual(
            (loaded.metabolites.get_by_id("m1437[C_em]").formula, loaded.metabolites.get_by_id("m1437[C_em]").charge),
            ("C47H80N7O18P3S", -4),
        )
        self.assertEqual(loaded.reactions.get_by_id("R1412").check_mass_balance(), {})
        self.assertEqual(curation["activation_state"], "blocked")

    def test_current_main_artifact_is_fresh(self) -> None:
        self.assertEqual(
            json.loads(CURRENT_MAIN_REPORT.read_text(encoding="utf-8")),
            audit_handoff(CURRENT_MAIN_MODEL),
        )


if __name__ == "__main__":
    unittest.main()
