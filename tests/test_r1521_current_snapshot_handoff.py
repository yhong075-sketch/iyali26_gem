"""Regression tests for the blocked, current-main R1521 handoff."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model, write_sbml_model

import scripts.r1521_current_snapshot_handoff as r1521_handoff
from scripts.gem_annotate.patches import (
    CoAProtonationActivationBlocked,
    _model_snapshot_fingerprint,
)
from scripts.r1521_current_snapshot_handoff import (
    COA_CURATION_PATH,
    HANDOFF_PATH,
    R1521CurrentSnapshotError,
    _complete_mass_balance,
    _local_counterfactual_on_copy,
    _validate_output_path,
    apply_r1521_current_snapshot_handoff,
    audit_r1521_current_snapshot,
    load_r1521_current_snapshot_handoff,
    validate_r1521_current_snapshot_model,
)


REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE_MODEL_ENV = os.environ.get("IYALI26_SOURCE_MODEL")
CURRENT_MAIN_MODEL = Path(SOURCE_MODEL_ENV) if SOURCE_MODEL_ENV else None
CURRENT_MAIN_AUDIT = REPOSITORY / "artifacts" / "r1521_current_snapshot_handoff_20260818.json"


class R1521EvidenceContractTests(unittest.TestCase):
    def test_frontier_sources_are_complete_disjoint_and_fail_closed(self) -> None:
        handoff = load_r1521_current_snapshot_handoff()
        dependencies = handoff["evidence_dependencies"]
        self.assertEqual(len(dependencies), 4)
        documents = (
            json.loads(
                (REPOSITORY / dependencies["r1521_rhea_frontier_source_audit"]["path"])
                .read_text(encoding="utf-8")
            ),
            json.loads(
                (REPOSITORY / dependencies["r1521_kegg_mnx_frontier_source_audit"]["path"])
                .read_text(encoding="utf-8")
            ),
            json.loads(
                (REPOSITORY / dependencies["r1521_unresolved_frontier_source_audit"]["path"])
                .read_text(encoding="utf-8")
            ),
        )
        for document in documents:
            self.assertEqual(document["activation_state"], "blocked")
            self.assertFalse(document["ready_for_activation"])
            self.assertTrue(document["production_apply_forbidden"])
            self.assertFalse(document["isolated_proton_edits_authorized"])
        batches = (documents[0]["rows"], documents[1]["rows"], documents[2]["reactions"])
        self.assertTrue(all(not row["h_edit_authorized"] for rows in batches for row in rows))
        id_sets = [
            {row["reaction_id"] for row in rows}
            for rows in batches
        ]
        self.assertEqual([len(ids) for ids in id_sets], [53, 22, 21])
        self.assertFalse(id_sets[0] & id_sets[1])
        self.assertFalse(id_sets[0] & id_sets[2])
        self.assertFalse(id_sets[1] & id_sets[2])
        audited_ids = set().union(*id_sets)
        self.assertEqual(len(audited_ids), 96)
        coverage = handoff["frontier_closure_gate"]["source_audit_coverage"]
        self.assertEqual(
            hashlib.sha256(
                json.dumps(sorted(audited_ids), separators=(",", ":")).encode()
            ).hexdigest(),
            coverage["audited_regression_reaction_ids_sha256"],
        )
        self.assertEqual(coverage["global_incident_not_source_audited_count"], 48)
        self.assertEqual(coverage["atomic_newly_unbalanced_not_audited"], 51)
        self.assertEqual(
            len(coverage["global_incident_not_source_audited_reaction_ids"]), 48
        )
        self.assertEqual(
            len(coverage["atomic_newly_unbalanced_not_audited_reaction_ids"]), 51
        )
        self.assertTrue(
            {"R613", "R2076", "R1708", "R2004"}
            <= set(coverage["atomic_newly_unbalanced_not_audited_reaction_ids"])
        )
        self.assertFalse(coverage["all_atomic_frontier_identities_closed"])
        self.assertFalse(coverage["isolated_proton_edits_authorized"])

        gpr = handoff["local_counterfactual"]["gpr_evidence"]
        self.assertEqual(
            (gpr["systematic_id"], gpr["established_name"], gpr["legacy_locus_tag"]),
            ("YALI1E18441g", "MFE2", "YALI0E15378g"),
        )
        self.assertEqual(gpr["protein_accessions"], ["Q9P4D9", "F2Z6I5"])
        self.assertEqual(gpr["evidence_status"]["exact_c26_substrate_turnover"], "unverified")
        r364 = handoff["frontier_closure_gate"]["r364_symbolic_contract"]
        self.assertEqual((r364["canonical_h_coefficient"], r364["model_h_coefficient"]), (1, 1))
        self.assertEqual(r364["balance_assessment_status"], "unassessable_symbolic_moiety")
        self.assertFalse(r364["h_edit_authorized"])
        self.assertFalse(handoff["activation"]["ready_for_activation"])

    def test_any_evidence_dependency_drift_is_rejected(self) -> None:
        handoff = load_r1521_current_snapshot_handoff()
        original = r1521_handoff._canonical_json_file_digest
        for dependency in handoff["evidence_dependencies"].values():
            bad_path = (REPOSITORY / dependency["path"]).resolve()
            with self.subTest(path=bad_path), patch(
                "scripts.r1521_current_snapshot_handoff._canonical_json_file_digest",
                side_effect=lambda path, bad=bad_path: (
                    "0" * 64 if path.resolve() == bad else original(path)
                ),
            ):
                with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                    load_r1521_current_snapshot_handoff()

    def test_exact_balance_assessor_rejects_incomplete_or_invalid_tuples(self) -> None:
        model = Model("balance-assessor")
        left = Metabolite("left", formula="C1H2", charge=0)
        right = Metabolite("right", formula="C1H2", charge=0)
        reaction = Reaction("test")
        reaction.add_metabolites({left: -1, right: 1})
        model.add_reactions([reaction])
        self.assertEqual(_complete_mass_balance(reaction), {})
        left.formula = None
        self.assertIsNone(_complete_mass_balance(reaction))
        left.formula = "C6H5O6*3"
        self.assertIsNone(_complete_mass_balance(reaction))
        left.formula = "C6H5O6foo"
        self.assertIsNone(_complete_mass_balance(reaction))
        left.formula = "C6H5O6?"
        self.assertIsNone(_complete_mass_balance(reaction))
        left.formula = "C1H2"
        left.charge = None
        self.assertIsNone(_complete_mass_balance(reaction))

    def test_output_protects_every_authoritative_runtime_input(self) -> None:
        handoff = load_r1521_current_snapshot_handoff()
        source = Path("/tmp/explicit-r1521-source.xml")
        protected = [
            source,
            HANDOFF_PATH,
            COA_CURATION_PATH,
            *(REPOSITORY / row["path"] for row in handoff["evidence_dependencies"].values()),
        ]
        for path in protected:
            with self.subTest(path=path), self.assertRaises(R1521CurrentSnapshotError):
                _validate_output_path(path, source, handoff)


@unittest.skipUnless(SOURCE_MODEL_ENV, "requires explicit IYALI26_SOURCE_MODEL")
class R1521CurrentSnapshotHandoffTests(unittest.TestCase):
    def setUp(self) -> None:
        self.handoff = load_r1521_current_snapshot_handoff()
        self.model = read_sbml_model(str(CURRENT_MAIN_MODEL))

    def test_current_sha_fingerprint_legacy_tuple_and_copy_inventory_are_exact(self) -> None:
        report = audit_r1521_current_snapshot(CURRENT_MAIN_MODEL)
        self.assertTrue(report["source_snapshot"]["sha256_verified"])
        self.assertTrue(report["source_snapshot"]["fingerprint_verified"])
        self.assertEqual(
            report["source_snapshot"]["declared_model_locator"],
            self.handoff["source_model"],
        )
        self.assertTrue(
            report["source_snapshot"]
            ["physical_input_locator_recorded_in_execution_metadata"]
        )
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
        self.assertEqual(candidate["gpr_evidence"]["systematic_id"], "YALI1E18441g")
        self.assertTrue(any("stale" in blocker.lower() for blocker in self.handoff["blockers"]))

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
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "drifted.xml"
            write_sbml_model(self.model, str(path))
            with self.assertRaises(R1521CurrentSnapshotError):
                audit_r1521_current_snapshot(path)

    def test_audit_requires_an_explicit_source_model_path(self) -> None:
        with self.assertRaisesRegex(R1521CurrentSnapshotError, "explicit R1521 source model"):
            audit_r1521_current_snapshot()

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
        dependency_paths = [
            REPOSITORY / entry["path"]
            for entry in self.handoff["evidence_dependencies"].values()
        ]
        for protected in (CURRENT_MAIN_MODEL, HANDOFF_PATH, COA_CURATION_PATH, *dependency_paths):
            with self.subTest(protected=protected):
                with self.assertRaises(R1521CurrentSnapshotError):
                    _validate_output_path(protected, CURRENT_MAIN_MODEL)

    def test_global_nad_gate_drift_aborts_the_audit(self) -> None:
        drifted = {
            "matches_expected": False,
            "mismatched_fields": ["newly_unbalanced_count"],
        }
        with patch(
            "scripts.r1521_current_snapshot_handoff._global_nad_migration_audit",
            return_value=drifted,
        ):
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "global NAD"):
                audit_r1521_current_snapshot(CURRENT_MAIN_MODEL)

    def test_full_frontiers_fail_closed_on_regressions_and_incomplete_tuples(self) -> None:
        report = audit_r1521_current_snapshot(CURRENT_MAIN_MODEL)
        balances = report["source_partial_known_formula_residual"]
        self.assertEqual(set(balances), {"R1504", "R1521", "R80"})
        self.assertTrue(all(row["matches_contract"] for row in balances.values()))
        self.assertFalse(balances["R1504"]["exact_balance_assessable"])
        self.assertFalse(balances["R1521"]["exact_balance_assessable"])
        self.assertTrue(balances["R80"]["exact_balance_assessable"])
        self.assertTrue(report["candidate_balances_named_adjacent_reactions_exactly"])
        self.assertEqual(report["candidate_adjacent_exact_balance"], {"R1504": {}, "R1521": {}, "R80": {}})
        local_gate = report["local_counterfactual_full_model_gate"]
        self.assertTrue(local_gate["matches_expected"], local_gate)
        self.assertFalse(local_gate["zero_new_regressions"])
        self.assertFalse(local_gate["frontier_closed"])
        self.assertEqual(
            {key: local_gate["actual"][key] for key in ("migrated_metabolite_count", "incident_reaction_count", "newly_unbalanced_count", "repaired_count")},
            {"migrated_metabolite_count": 8, "incident_reaction_count": 70, "newly_unbalanced_count": 30, "repaired_count": 0},
        )
        self.assertEqual(local_gate["actual"]["unassessable_reaction_count"], 23)
        global_gate = report["global_nad_migration_gate"]
        self.assertEqual(self.handoff["global_nad_migration_gate"]["baseline"], "source_snapshot_only; no CoA/local transaction pre-applied")
        self.assertTrue(global_gate["matches_expected"], global_gate)
        self.assertFalse(global_gate["frontier_closed"])
        self.assertEqual(
            {key: global_gate["actual"][key] for key in ("migrated_metabolite_count", "incident_reaction_count", "changed_residual_count", "newly_unbalanced_count", "repaired_count", "unassessable_reaction_count")},
            {"migrated_metabolite_count": 11, "incident_reaction_count": 121, "changed_residual_count": 83, "newly_unbalanced_count": 73, "repaired_count": 2, "unassessable_reaction_count": 32},
        )
        self.assertEqual(
            self.handoff["global_nad_migration_gate"]
            ["historical_partial_formula_baseline"]["changed_residual_count"],
            114,
        )
        atomic_gate = report["atomic_identity_frontier_gate"]
        self.assertTrue(atomic_gate["matches_expected"], atomic_gate)
        self.assertFalse(atomic_gate["frontier_closed"])
        self.assertEqual(
            {key: atomic_gate["actual"][key] for key in ("migrated_metabolite_count", "incident_reaction_count", "newly_unbalanced_count", "repaired_count", "unassessable_reaction_count")},
            {"migrated_metabolite_count": 36, "incident_reaction_count": 303, "newly_unbalanced_count": 146, "repaired_count": 2, "unassessable_reaction_count": 75},
        )
        self.assertIn("R364", atomic_gate["actual"]["unassessable_reaction_ids"])
        self.assertEqual(
            report["runtime_inputs"]["atomic_identity_curation"]["declared_locator"],
            "data/coa_protonation_curation.json",
        )
        self.assertFalse(report["frontier_closed"])
        self.assertFalse(report["production_gate_passed"])

    def test_manifest_drift_and_reaction_drift_are_rejected(self) -> None:
        changed = copy.deepcopy(self.handoff)
        changed["activation"]["state"] = "approved"
        with tempfile.TemporaryDirectory() as temp_dir:
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

    def test_caller_supplied_handoff_cannot_bypass_evidence_checks(self) -> None:
        with patch(
            "scripts.r1521_current_snapshot_handoff._canonical_json_file_digest",
            return_value="0" * 64,
        ):
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                validate_r1521_current_snapshot_model(self.model, self.handoff)
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                audit_r1521_current_snapshot(
                    CURRENT_MAIN_MODEL, handoff=self.handoff
                )
            with self.assertRaisesRegex(R1521CurrentSnapshotError, "evidence dependency"):
                apply_r1521_current_snapshot_handoff(self.model, self.handoff)

    def test_blocked_apply_rolls_back_before_any_write_and_is_idempotent(self) -> None:
        before = _model_snapshot_fingerprint(self.model)
        for _ in range(2):
            with self.assertRaises(CoAProtonationActivationBlocked):
                apply_r1521_current_snapshot_handoff(self.model, self.handoff)
            self.assertEqual(_model_snapshot_fingerprint(self.model), before)

    def test_read_only_audit_preserves_fba_objective_and_repeats_deterministically(self) -> None:
        first = audit_r1521_current_snapshot(CURRENT_MAIN_MODEL)
        second = audit_r1521_current_snapshot(CURRENT_MAIN_MODEL)
        self.assertEqual(first, second)
        self.assertTrue(first["objective_unchanged"])
        self.assertAlmostEqual(float(first["objective_before"]), float(first["objective_after"]))

    def test_candidate_sbml_round_trip_keeps_3r_identity_and_blocked_handoff(self) -> None:
        candidate = _local_counterfactual_on_copy(self.model, self.handoff)
        with tempfile.TemporaryDirectory() as temp_dir:
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
        expected = audit_r1521_current_snapshot(CURRENT_MAIN_MODEL)
        observed = json.loads(CURRENT_MAIN_AUDIT.read_text(encoding="utf-8"))
        self.assertEqual(observed, expected)
        self.assertEqual(
            observed["evidence_role"], "local_cross_check_not_hpcc_acceptance"
        )
        self.assertEqual(
            observed["handoff_contract_sha256"],
            self.handoff["target_contract_sha256"],
        )
        self.assertEqual(observed["evidence_dependencies"], self.handoff["evidence_dependencies"])


if __name__ == "__main__":
    unittest.main()
