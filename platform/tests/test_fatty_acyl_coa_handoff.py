"""Regression tests for the machine-readable fatty acyl-CoA handoff."""

from __future__ import annotations

import copy
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

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
SOURCE_MODEL_ENV = os.environ.get("IYALI26_SOURCE_MODEL")
CURRENT_MAIN_MODEL = Path(SOURCE_MODEL_ENV) if SOURCE_MODEL_ENV else None
MODEL_PATH = CURRENT_MAIN_MODEL
CURRENT_MAIN_REPORT = REPOSITORY / "artifacts" / "fatty_acyl_coa_handoff_20260818.json"


class FattyAcylCoAHandoffTests(unittest.TestCase):
    @unittest.skipUnless(SOURCE_MODEL_ENV, "requires explicit IYALI26_SOURCE_MODEL")
    def test_handoff_matches_pipeline_inputs_and_complete_inventory(self) -> None:
        handoff = load_handoff()
        validate_handoff_inputs(handoff, source_model_path=MODEL_PATH)
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
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "handoff.json"
            path.write_text(json.dumps(changed), encoding="utf-8")
            with self.assertRaises(HandoffError):
                load_handoff(path)

    def test_d_review_inventory_ledgers_and_science_are_fail_closed(self) -> None:
        review = load_handoff()["d_connected_component_review"]
        self.assertEqual(
            [row["id"] for row in review["group_inventory"]],
            [
                "acetate", "succinate", "acetyl_coa", "succinyl_coa",
                "malonyl_coa", "pyruvate", "oxaloacetate",
            ],
        )
        self.assertEqual(
            len({metabolite_id for row in review["group_inventory"] for metabolite_id in row["expected_ids"]}),
            25,
        )
        ledgers = review["components"]
        for name, counts, digest in (
            ("r1708", (15, 100, 84, 16, 76, 64, 3, 9, 8, 0), "076b4c3f9cdcb80d7dc60ed8be54335280ec2c4722d73d805461c1243acc79ed"),
            ("r2004", (14, 119, 100, 19, 91, 78, 0, 13, 9, 0), "8c3499d76f0999ce6b98c819ac8e6a7ff1a5329b785c5cef29dff0b8dcdd287a"),
            ("union", (25, 155, 133, 22, 117, 95, 3, 19, 16, 0), "731707e29d8385cdccb8188e77f8778f2435cdd6618da1fd339d94da92e28839"),
        ):
            expected = ledgers[name]["expected"]
            self.assertEqual(
                tuple(expected[field] for field in (
                    "migrated_metabolite_count", "incident_reaction_count",
                    "evaluable_reaction_count", "unassessable_reaction_count",
                    "changed_residual_count", "newly_unbalanced_count",
                    "repaired_count", "still_unbalanced_changed_count",
                    "balanced_unchanged_count", "unbalanced_unchanged_count",
                )),
                counts,
            )
            self.assertEqual(expected["ledger_sha256"], digest)
        union = ledgers["union"]
        self.assertEqual(len(union["newly_unbalanced_reaction_ids"]), 95)
        self.assertEqual(len(union["still_unbalanced_changed_reaction_ids"]), 19)
        self.assertEqual(len(union["expected"]["unassessable_reaction_ids"]), 22)
        self.assertEqual(union["expected"]["repaired_reaction_ids"], ["R2077", "R2269", "R_CYSS_m"])

        r1708 = review["reaction_reviews"]["R1708"]
        self.assertEqual(r1708["ph_7_3_aggregate_each_side"], {"formula": "C27H38N7O21P3S", "charge": -6})
        self.assertEqual(r1708["explicit_proton_coefficient"], 0)
        gene = r1708["gpr_evidence"]
        self.assertEqual(
            (gene["systematic_identifier"], gene["established_name"], gene["legacy_systematic_identifier"], gene["uniprot_accession"], gene["protein_function"]),
            ("YALI1E36437g", "ACH1", "YALI0E30965g", "Q6C3Z9", "acetate:succinate CoA-transferase"),
        )
        self.assertIn("experimentally_supported_function", gene["evidence_status"])
        self.assertIn("unverified", gene["localization_evidence"])
        self.assertIn("homology", gene["database_conflict"])
        self.assertEqual(
            r1708["preserve_database_reaction_ids"]["rhea"],
            {
                "master": "35711", "left_to_right": "35712",
                "right_to_left": "35713", "bidirectional": "35714",
            },
        )
        provenance = r1708["provenance"]
        self.assertEqual(provenance["accessed_on"], "2026-08-18")
        self.assertEqual(
            (
                provenance["identifier_mapping"]["model_identifier"],
                provenance["identifier_mapping"]["s2_identifier"],
                provenance["identifier_mapping"]["legacy_identifier"],
                provenance["identifier_mapping"]["ncbi_gene_id"],
            ),
            ("YALI1E36437g", "YALI1_E36437g", "YALI0E30965g", "2911920"),
        )
        self.assertIn("partial", provenance["identifier_mapping"]["ncbi_feature_annotation_status"])
        self.assertEqual(
            provenance["functional_evidence"]["evidence_status"],
            "experimentally_supported_condition_specific",
        )
        self.assertEqual(provenance["proteomics_evidence"]["matches"], 35)
        self.assertEqual(provenance["proteomics_evidence"]["coverage_percent"], 66)
        self.assertEqual(provenance["proteomics_evidence"]["supplementary_table"], "Table S2")
        self.assertIn(
            "not_function_or_organelle_localization",
            provenance["proteomics_evidence"]["evidence_status"],
        )
        self.assertIsNone(r1708["second_isoenzyme_assignment"])
        self.assertEqual(r1708["r89_adjacent_candidate"]["status"], "blocked_not_in_reaction_corrections")
        self.assertFalse(r1708["production_authorized"])
        self.assertFalse(r1708["production_decision_applied"])
        self.assertFalse(r1708["r89_adjacent_candidate"]["production_authorized"])
        self.assertFalse(r1708["r89_adjacent_candidate"]["production_decision_applied"])

        r2004 = review["reaction_reviews"]["R2004"]
        self.assertEqual(r2004["current_malonyl_coa_metabolite_id"], "m200[C_cy]")
        self.assertEqual(r2004["forbidden_historical_metabolite_id"], "m1855[C_cy]")
        self.assertEqual(r2004["ph_7_3_aggregate_each_side"], {"formula": "C27H36N7O22P3S", "charge": -6})
        self.assertEqual(r2004["explicit_proton_coefficient"], 0)
        self.assertTrue(all(r2004["gpr_evidence"][field] is None for field in ("systematic_identifier", "established_name", "protein_function", "model_gpr")))
        self.assertEqual(r2004["gpr_evidence"]["source_model_gene_reaction_rule"], "")
        provenance = r2004["gpr_evidence"]["search_provenance"]
        self.assertEqual(provenance["accessed_on"], "2026-08-18")
        self.assertEqual(len(provenance["zero_hit_searches"]), 6)
        self.assertTrue(all(row["hit_count"] == 0 for row in provenance["zero_hit_searches"]))
        self.assertEqual(len(provenance["positive_controls"]), 5)
        self.assertTrue(all(row["hit_count"] > 0 for row in provenance["positive_controls"]))
        self.assertEqual(
            r2004["non_yarrowia_biochemical_evidence"]["evidence_status"],
            "partial_non_yarrowia_not_gpr_evidence",
        )
        self.assertIn("cannot justify replacing", r2004["ec_scope_limitation"])
        self.assertEqual(r2004["preferred_candidate_decision"], "remove")
        self.assertFalse(r2004["production_authorized"])
        self.assertFalse(r2004["production_decision_applied"])
        self.assertEqual(
            {row["decision"]: row["status"] for row in r2004["candidates"]},
            {
                "retain_reversible_without_gpr": "allowed_for_sensitivity_only",
                "remove": "recommended_fail_closed_candidate",
                "replace_with_R00930": "forbidden_wrong_substrates",
            },
        )
        self.assertFalse(review["frontier_closed"])
        self.assertFalse(review["production_gate_passed"])
        self.assertFalse(review["ready_for_activation"])
        self.assertTrue(review["human_gate_required"])

    def test_d_nested_contract_drift_fails_closed(self) -> None:
        paths_and_values = (
            (("status",), "approved"),
            (("reaction_reviews", "R1708"), None),
            (("reaction_reviews", "R1708", "production_authorized"), True),
            (("reaction_reviews", "R1708", "gpr_evidence", "systematic_identifier"), "WRONG"),
            (("reaction_reviews", "R1708", "provenance", "proteomics_evidence", "supplementary_table"), "Table S1"),
            (("reaction_reviews", "R1708", "second_isoenzyme_assignment"), "guessed_gene"),
            (("reaction_reviews", "R1708", "r89_adjacent_candidate", "production_decision_applied"), True),
            (("reaction_reviews", "R1708", "explicit_proton_coefficient"), 1),
            (("reaction_reviews", "R1708", "provenance", "accessed_on"), None),
            (("reaction_reviews", "R2004", "production_decision_applied"), True),
            (("reaction_reviews", "R2004", "current_malonyl_coa_metabolite_id"), "m1855[C_cy]"),
            (("reaction_reviews", "R2004", "preferred_candidate_decision"), "retain_reversible_without_gpr"),
            (("reaction_reviews", "R2004", "gpr_evidence", "systematic_identifier"), "guessed_gene"),
            (("reaction_reviews", "R2004", "gpr_evidence", "evidence_status"), "experimentally_verified"),
            (("reaction_reviews", "R2004", "candidates", 0, "production_authorized"), True),
            (("reaction_reviews", "R9999"), {"production_authorized": True, "production_decision_applied": True}),
        )
        for path, value in paths_and_values:
            with self.subTest(path=path):
                changed = copy.deepcopy(load_handoff())
                target = changed["d_connected_component_review"]
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                with tempfile.TemporaryDirectory() as temp_dir:
                    changed_path = Path(temp_dir) / "handoff.json"
                    changed_path.write_text(json.dumps(changed), encoding="utf-8")
                    with self.assertRaises(HandoffError):
                        load_handoff(changed_path)

    @unittest.skipUnless(SOURCE_MODEL_ENV, "requires explicit IYALI26_SOURCE_MODEL")
    def test_d_focal_activation_blockers_cannot_be_removed(self) -> None:
        handoff = load_handoff()
        curation = coa_patches.load_coa_protonation_curation(
            source_model_path=MODEL_PATH
        )
        for reaction_id in ("R1708", "R2004"):
            with self.subTest(reaction_id=reaction_id):
                changed = copy.deepcopy(curation)
                changed["activation_blockers"] = [
                    row for row in changed["activation_blockers"]
                    if row["reaction_id"] != reaction_id
                ]
                changed["target_balance_reaction_ids"].remove(reaction_id)
                with patch(
                    "scripts.fatty_acyl_coa_handoff.load_coa_protonation_curation",
                    return_value=changed,
                ):
                    with self.assertRaises(HandoffError):
                        validate_handoff_inputs(handoff, source_model_path=MODEL_PATH)
        missing_residual = copy.deepcopy(curation)
        next(
            row for row in missing_residual["activation_blockers"]
            if row["reaction_id"] == "R1708"
        ).pop("residual_after_curated_tuples")
        with patch(
            "scripts.fatty_acyl_coa_handoff.load_coa_protonation_curation",
            return_value=missing_residual,
        ):
            with self.assertRaises(HandoffError):
                validate_handoff_inputs(handoff, source_model_path=MODEL_PATH)

    def test_d_union_frontier_id_drift_fails_closed(self) -> None:
        changed = copy.deepcopy(load_handoff())
        changed["d_connected_component_review"]["components"]["union"]["newly_unbalanced_reaction_ids"].pop()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "handoff.json"
            path.write_text(json.dumps(changed), encoding="utf-8")
            with self.assertRaises(HandoffError):
                load_handoff(path)

    def test_malformed_activation_is_a_handoff_error(self) -> None:
        changed = copy.deepcopy(load_handoff())
        changed["activation"] = None
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "handoff.json"
            path.write_text(json.dumps(changed), encoding="utf-8")
            with self.assertRaises(HandoffError):
                load_handoff(path)

    def test_malformed_d_structures_are_handoff_errors(self) -> None:
        cases = (
            (("d_connected_component_review", "group_inventory"), [None], False),
            (("d_connected_component_review", "components", "r1708"), None, False),
            (("d_connected_component_review", "components", "r1708", "expected", "changed_residual_count"), None, True),
            (("d_connected_component_review", "components", "r2004", "expected", "repaired_count"), False, False),
            (("d_connected_component_review", "components", "union", "expected", "newly_unbalanced_count"), None, True),
            (("authoritative_inputs",), None, False),
        )
        for path, value, delete in cases:
            with self.subTest(path=path):
                changed = copy.deepcopy(load_handoff())
                target = changed
                for key in path[:-1]:
                    target = target[key]
                if delete:
                    del target[path[-1]]
                else:
                    target[path[-1]] = value
                with tempfile.TemporaryDirectory() as temp_dir:
                    changed_path = Path(temp_dir) / "handoff.json"
                    changed_path.write_text(json.dumps(changed), encoding="utf-8")
                    with self.assertRaises(HandoffError):
                        load_handoff(changed_path)

    def test_report_output_cannot_alias_an_authoritative_input(self) -> None:
        handoff = load_handoff()
        r1521 = json.loads(
            (REPOSITORY / "data" / "r1521_current_snapshot_handoff.json")
            .read_text(encoding="utf-8")
        )
        model_path = CURRENT_MAIN_MODEL or REPOSITORY / "explicit-source-model.xml"
        protected = [
            model_path,
            HANDOFF_PATH,
            REPOSITORY / handoff["authoritative_inputs"]["coa_protonation_curation"]["path"],
            REPOSITORY / handoff["authoritative_inputs"]["lipid_moiety_ledger_spec"]["path"],
            REPOSITORY / handoff["authoritative_inputs"]["r1521_current_snapshot_handoff"]["path"],
            Path(handoff["authoritative_inputs"]["human_review_table"]["path"]),
            *(
                REPOSITORY / dependency["path"]
                for dependency in r1521["evidence_dependencies"].values()
            ),
        ]
        for path in protected:
            with self.subTest(path=path):
                with self.assertRaises(HandoffError):
                    _validate_output_path(path, model_path)

    @unittest.skipUnless(SOURCE_MODEL_ENV, "requires explicit IYALI26_SOURCE_MODEL")
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
        d_review = report["current_snapshot_component_reviews"]["D_R1708_R2004"]
        declared = load_handoff()["d_connected_component_review"]
        self.assertTrue(d_review["deterministic_recomputation"])
        for name in ("r1708", "r2004", "union"):
            self.assertTrue(d_review["ledgers"][name]["matches_expected"])
            self.assertEqual(
                d_review["ledgers"][name]["actual"],
                declared["components"][name]["expected"],
            )
        self.assertEqual(
            d_review["d_only_union_focal_balances"],
            {"R1708": {}, "R2004": {}, "R613": {"H": -5.0, "charge": -5.0}, "R2076": {}},
        )
        self.assertEqual(
            d_review["complete_coa_candidate_focal_balances"],
            {"R1708": {}, "R2004": {}, "R613": {"H": -1.0, "charge": -1.0}, "R2076": {"H": 4.0, "charge": 4.0}},
        )
        self.assertEqual(len(d_review["union_frontier_reaction_ids"]["newly_unbalanced_reaction_ids"]), 95)
        self.assertEqual(len(d_review["union_frontier_reaction_ids"]["still_unbalanced_changed_reaction_ids"]), 19)
        self.assertFalse(d_review["frontier_closed"])
        self.assertFalse(d_review["production_gate_passed"])
        self.assertFalse(d_review["ready_for_activation"])
        self.assertTrue(d_review["human_gate_required"])

    @unittest.skipUnless(SOURCE_MODEL_ENV, "requires explicit IYALI26_SOURCE_MODEL")
    def test_candidate_survives_sbml_roundtrip_without_relaxing_blocked_state(self) -> None:
        model = read_sbml_model(str(MODEL_PATH))
        curation = coa_patches.load_coa_protonation_curation(
            source_model_path=MODEL_PATH
        )
        coa_patches._apply_coa_protonation_curation(model, curation)
        with tempfile.TemporaryDirectory() as temp_dir:
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

    @unittest.skipUnless(SOURCE_MODEL_ENV, "requires explicit IYALI26_SOURCE_MODEL")
    def test_current_main_artifact_is_fresh(self) -> None:
        observed = json.loads(CURRENT_MAIN_REPORT.read_text(encoding="utf-8"))
        expected = audit_handoff(CURRENT_MAIN_MODEL)
        self.assertEqual(observed, expected)
        self.assertEqual(observed["evidence_role"], "local_cross_check_not_hpcc_acceptance")


if __name__ == "__main__":
    unittest.main()
