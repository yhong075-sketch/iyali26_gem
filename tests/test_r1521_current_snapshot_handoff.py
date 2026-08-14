"""Regression tests for the blocked, current-main R1521 handoff."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra import Metabolite
from cobra.io import read_sbml_model, write_sbml_model

from scripts.gem_annotate.patches import (
    CoAProtonationActivationBlocked,
    _model_snapshot_fingerprint,
)
from scripts.r1521_current_snapshot_handoff import (
    ER_EVIDENCE_PATH,
    HANDOFF_PATH,
    R1521CurrentSnapshotError,
    _local_counterfactual_on_copy,
    _validate_output_path,
    apply_r1521_current_snapshot_handoff,
    audit_r1521_current_snapshot,
    load_r1521_current_snapshot_handoff,
    validate_r1521_current_snapshot_model,
)


REPOSITORY = Path(__file__).resolve().parents[1]
CURRENT_MAIN_MODEL = REPOSITORY.parent / "iyali26_gem" / "model.xml"
CURRENT_MAIN_AUDIT = REPOSITORY / "artifacts" / "r1521_current_snapshot_handoff_20260806.json"


@unittest.skipUnless(CURRENT_MAIN_MODEL.exists(), "requires the declared current-main model")
class R1521CurrentSnapshotHandoffTests(unittest.TestCase):
    def setUp(self) -> None:
        self.handoff = load_r1521_current_snapshot_handoff()
        self.model = read_sbml_model(str(CURRENT_MAIN_MODEL))

    def test_current_sha_fingerprint_legacy_tuple_and_copy_inventory_are_exact(self) -> None:
        report = audit_r1521_current_snapshot()
        self.assertTrue(report["source_snapshot"]["sha256_verified"])
        self.assertTrue(report["source_snapshot"]["fingerprint_verified"])
        self.assertTrue(report["validation"]["valid"], report["validation"]["errors"])
        self.assertEqual(
            self.handoff["legacy_metabolite_contract"]["copies"][0]["legacy_tuple"],
            {"formula": None, "charge": 0},
        )
        self.assertFalse(report["ready_for_activation"])
        self.assertEqual(report["activation"]["state"], "blocked")
        candidate = self.handoff["local_counterfactual"]
        self.assertIn("MFE2", candidate["evidence"])
        self.assertIn("3R", candidate["evidence"])
        self.assertIn("3S", candidate["supersedes"])
        self.assertEqual(candidate["rhea"]["hydratase_direction"], "39213")
        self.assertEqual(candidate["rhea"]["dehydrogenase_bidirectional"], "78638")
        self.assertIn("stale metadata", self.handoff["blockers"][2])

    def test_unknown_third_tuple_and_extra_same_identity_copy_fail_closed(self) -> None:
        changed = self.model.copy()
        changed.metabolites.get_by_id("m1546[C_pe]").formula = "C47H86N7O18P3S"
        errors = validate_r1521_current_snapshot_model(changed, self.handoff)["errors"]
        self.assertTrue(any("third tuple" in error for error in errors), errors)

        changed = self.model.copy()
        changed.metabolites.get_by_id("m226[C_pe]").charge = -4
        errors = validate_r1521_current_snapshot_model(changed, self.handoff)["errors"]
        self.assertTrue(any("m226[C_pe]: third tuple" in error for error in errors), errors)

        changed = self.model.copy()
        changed.add_metabolites(
            [
                Metabolite(
                    "m1546_extra[C_pe]",
                    name="(R)-3-hydroxyhexacosanoyl-CoA_",
                    formula=None,
                    charge=0,
                    compartment="C_pe",
                )
            ]
        )
        errors = validate_r1521_current_snapshot_model(changed, self.handoff)["errors"]
        self.assertTrue(any("copy completeness" in error for error in errors), errors)

        changed = self.model.copy()
        changed.metabolites.get_by_id("m27[C_mi]").charge = -1
        errors = validate_r1521_current_snapshot_model(changed, self.handoff)["errors"]
        self.assertTrue(any("nad_plus: m27[C_mi] has third tuple" in error for error in errors), errors)

        changed = self.model.copy()
        changed.add_metabolites(
            [
                Metabolite(
                    "m27_extra[C_mi]",
                    name="NAD_C21H27N7O14P2",
                    formula="C21H27N7O14P2",
                    charge=0,
                    compartment="C_mi",
                )
            ]
        )
        errors = validate_r1521_current_snapshot_model(changed, self.handoff)["errors"]
        self.assertTrue(any("nad_plus: copy completeness" in error for error in errors), errors)

    def test_source_sha_or_fingerprint_drift_fails_before_counterfactual(self) -> None:
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temp_dir:
            path = Path(temp_dir) / "drifted.xml"
            write_sbml_model(self.model, str(path))
            with self.assertRaises(R1521CurrentSnapshotError):
                audit_r1521_current_snapshot(path)

    def test_direct_script_entry_resolves_to_the_lipid_worktree(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/r1521_current_snapshot_handoff.py", "--help"],
            cwd=REPOSITORY,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--output", completed.stdout)

    def test_output_cannot_alias_any_authoritative_input(self) -> None:
        for protected in (CURRENT_MAIN_MODEL, HANDOFF_PATH, ER_EVIDENCE_PATH):
            with self.subTest(protected=protected):
                with self.assertRaises(R1521CurrentSnapshotError):
                    _validate_output_path(protected, CURRENT_MAIN_MODEL)

    def test_global_nad_gate_drift_aborts_the_audit(self) -> None:
        drifted = {"matches_expected": False, "checks": {"newly_unbalanced_count": False}}
        with patch(
            "scripts.r1521_current_snapshot_handoff._global_nad_migration_audit",
            return_value=drifted,
        ):
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "global NAD"):
                audit_r1521_current_snapshot()

    def test_adjacent_contracts_preserve_the_exact_current_absolute_balances(self) -> None:
        report = audit_r1521_current_snapshot()
        balances = report["adjacent_absolute_balance"]
        self.assertEqual(set(balances), {"R1504", "R1521", "R80"})
        self.assertTrue(all(row["matches_exactly"] for row in balances.values()))
        self.assertFalse(balances["R1504"]["is_balanced"])
        self.assertFalse(balances["R1521"]["is_balanced"])
        self.assertTrue(balances["R80"]["is_balanced"])
        self.assertTrue(report["candidate_closes_adjacent_reactions"])
        self.assertEqual(report["candidate_adjacent_absolute_balance"], {"R1504": {}, "R1521": {}, "R80": {}})
        global_gate = report["global_nad_migration_gate"]
        self.assertEqual(self.handoff["global_nad_migration_gate"]["baseline"], "source_snapshot_only; no CoA/local transaction pre-applied")
        self.assertTrue(global_gate["matches_expected"], global_gate)
        self.assertEqual(
            {key: global_gate["actual"][key] for key in ("incident_reaction_count", "changed_residual_count", "newly_unbalanced_count", "repaired_count")},
            {"incident_reaction_count": 121, "changed_residual_count": 114, "newly_unbalanced_count": 73, "repaired_count": 2},
        )

    def test_manifest_drift_and_reaction_drift_are_rejected(self) -> None:
        changed = copy.deepcopy(self.handoff)
        changed["activation"]["state"] = "approved"
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temp_dir:
            path = Path(temp_dir) / "handoff.json"
            path.write_text(json.dumps(changed), encoding="utf-8")
            with self.assertRaises(R1521CurrentSnapshotError):
                load_r1521_current_snapshot_handoff(path)

        changed_model = self.model.copy()
        changed_model.reactions.get_by_id("R1521").add_metabolites(
            {changed_model.metabolites.get_by_id("m118[C_pe]"): 1}
        )
        errors = validate_r1521_current_snapshot_model(changed_model, self.handoff)["errors"]
        self.assertTrue(any("R1521: stoichiometry drifted" in error for error in errors), errors)

        changed_model = self.model.copy()
        changed_model.reactions.get_by_id("R1521").upper_bound = 1000.0000005
        errors = validate_r1521_current_snapshot_model(changed_model, self.handoff)["errors"]
        self.assertTrue(any("R1521: bounds or reversibility drifted" in error for error in errors), errors)

        changed_model = self.model.copy()
        changed_model.add_metabolites(
            [
                Metabolite(
                    "m1546_3r_extra[C_pe]",
                    name="(3R)-3-hydroxyhexacosanoyl-CoA",
                    formula=None,
                    charge=0,
                    compartment="C_pe",
                )
            ]
        )
        errors = validate_r1521_current_snapshot_model(changed_model, self.handoff)["errors"]
        self.assertTrue(any("m1546[C_pe]: copy completeness" in error for error in errors), errors)

    def test_er_evidence_dependency_drift_is_rejected(self) -> None:
        with patch(
            "scripts.r1521_current_snapshot_handoff._canonical_json_file_digest",
            return_value="0" * 64,
        ):
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                load_r1521_current_snapshot_handoff()

    def test_caller_supplied_handoff_cannot_bypass_er_evidence_check(self) -> None:
        with patch(
            "scripts.r1521_current_snapshot_handoff._canonical_json_file_digest",
            return_value="0" * 64,
        ):
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                validate_r1521_current_snapshot_model(self.model, self.handoff)
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                audit_r1521_current_snapshot(handoff=self.handoff)
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                apply_r1521_current_snapshot_handoff(self.model, self.handoff)

    def test_blocked_apply_rolls_back_before_any_write_and_is_idempotent(self) -> None:
        before = _model_snapshot_fingerprint(self.model)
        for _ in range(2):
            with self.assertRaises(CoAProtonationActivationBlocked):
                apply_r1521_current_snapshot_handoff(self.model, self.handoff)
            self.assertEqual(_model_snapshot_fingerprint(self.model), before)

    def test_read_only_audit_preserves_fba_objective_and_repeats_deterministically(self) -> None:
        first = audit_r1521_current_snapshot()
        second = audit_r1521_current_snapshot()
        self.assertEqual(first, second)
        self.assertTrue(first["objective_unchanged"])
        self.assertAlmostEqual(float(first["objective_before"]), float(first["objective_after"]))

    def test_candidate_sbml_round_trip_keeps_3r_identity_and_blocked_handoff(self) -> None:
        candidate = _local_counterfactual_on_copy(self.model, self.handoff)
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temp_dir:
            path = Path(temp_dir) / "r1521-candidate-roundtrip.xml"
            write_sbml_model(candidate, str(path))
            reloaded = read_sbml_model(str(path))
        hydroxy_c26 = reloaded.metabolites.get_by_id("m1546[C_pe]")
        self.assertEqual(hydroxy_c26.name, "(3R)-3-hydroxyhexacosanoyl-CoA")
        self.assertEqual((hydroxy_c26.formula, hydroxy_c26.charge), ("C47H82N7O18P3S", -4))
        self.assertEqual(hydroxy_c26.annotation.get("chebi"), "CHEBI:76378")
        self.assertEqual(
            reloaded.reactions.get_by_id("R80").check_mass_balance(),
            {},
        )
        self.assertEqual(reloaded.reactions.get_by_id("R1504").check_mass_balance(), {})
        self.assertEqual(reloaded.reactions.get_by_id("R1521").check_mass_balance(), {})
        self.assertEqual(self.handoff["activation"]["state"], "blocked")

    def test_checked_in_audit_is_fresh_and_records_evidence_provenance(self) -> None:
        expected = audit_r1521_current_snapshot()
        observed = json.loads(CURRENT_MAIN_AUDIT.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(
            observed["handoff_contract_sha256"],
            self.handoff["target_contract_sha256"],
        )
        self.assertEqual(observed["evidence_dependencies"], self.handoff["evidence_dependencies"])


if __name__ == "__main__":
    unittest.main()
