"""Strict, inactive CoA biochemical-pH protonation curation tests."""

from __future__ import annotations

import copy
import hashlib
import inspect
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model

import scripts.gem_annotate.patches as coa_patches

from scripts.gem_annotate.patches import (
    CoAProtonationActivationBlocked,
    CoAProtonationCurationError,
    apply_all_patches,
    audit_coa_protonation_curation,
    normalize_coa_protonation,
    validate_coa_protonation_curation,
)


REPOSITORY = Path(__file__).resolve().parents[1]
MODEL_PATH = Path(os.environ.get(
    "IYALI26_SOURCE_MODEL", REPOSITORY.parent / "iyali26_gem" / "model.xml"
))
COA_CURATION_PATH = REPOSITORY / "data" / "coa_protonation_curation.json"
LEDGER_SPEC_PATH = REPOSITORY / "data" / "lipid_moiety_ledger_spec.json"


def load_coa_protonation_curation(curation_path: Path | None = None) -> dict:
    return coa_patches.load_coa_protonation_curation(
        curation_path, source_model_path=MODEL_PATH
    )


CORE_TARGETS = {
    "coenzyme_a": (9, "C21H32N7O16P3S", -4),
    "lauroyl_coa": (5, "C33H54N7O17P3S", -4),
    "myristoyl_coa": (5, "C35H58N7O17P3S", -4),
    "palmitoyl_coa": (6, "C37H62N7O17P3S", -4),
    "palmitoleoyl_coa": (5, "C37H60N7O17P3S", -4),
    "stearoyl_coa": (5, "C39H66N7O17P3S", -4),
    "oleoyl_coa": (5, "C39H64N7O17P3S", -4),
    "linoleoyl_coa": (5, "C39H62N7O17P3S", -4),
}


def group_by_id(curation: dict, group_id: str) -> dict:
    return next(group for group in curation["groups"] if group["id"] == group_id)


def approved_minimal_curation() -> dict:
    """A self-contained approved fixture for exercising the gate, not production data."""
    curation = {
        "schema_version": 1,
        "activation_state": "approved",
        "source_stage": "in_memory_test_fixture",
        "activation_blockers": [],
        "target_balance_reaction_ids": ["R2076"],
        "groups": [
            {
                "id": "coenzyme_a",
                "identity_names": ["coenzyme a"],
                "expected_ids": ["coa_c"],
                "expected_copy_count": 1,
                "legacy_tuples": [
                    {
                        "ids": ["coa_c"],
                        "formula": "C21H36N7O16P3S",
                        "charge": 0,
                    }
                ],
                "target_tuple": {"formula": "C21H32N7O16P3S", "charge": -4},
            }
        ],
        "reaction_corrections": [
            {
                "reaction_id": "R2076",
                "metabolite_id": "h_c",
                "legacy_coefficient": 3,
                "target_coefficient": -1,
            }
        ],
    }
    source = build_balanced_fixture()
    curation["source_model_fingerprint"] = coa_patches._model_snapshot_fingerprint(source)
    target = source.copy()
    coa_patches._apply_coa_protonation_curation(target, curation)
    curation["target_model_fingerprint"] = coa_patches._model_snapshot_fingerprint(
        target
    )
    return curation


def build_balanced_fixture() -> Model:
    """Build a model where the exact R2076 correction preserves all gates."""
    model = Model("coa_curation_fixture")
    coa = Metabolite(
        "coa_c",
        name="coenzyme A_C21H36N7O16P3S",
        formula="C21H36N7O16P3S",
        charge=0,
        compartment="c",
    )
    methylmalonate = Metabolite("methylmalonate_c", formula="C4H4O4", charge=-2, compartment="c")
    methylmalonyl_coa = Metabolite(
        "methylmalonyl_coa_c",
        formula="C25H35N7O19P3S",
        charge=-5,
        compartment="c",
    )
    proton = Metabolite("h_c", formula="H", charge=1, compartment="c")
    water = Metabolite("h2o_c", formula="H2O", charge=0, compartment="c")
    r2076 = Reaction("R2076")
    r2076.add_metabolites(
        {
            methylmalonate: -1,
            coa: -1,
            proton: 3,
            methylmalonyl_coa: 1,
            water: 1,
        }
    )

    feed = Metabolite("feed_c", formula="H", charge=0, compartment="c")
    product = Metabolite("product_c", formula="H", charge=0, compartment="c")
    growth = Reaction("growth")
    growth.add_metabolites({feed: -1, product: 1})
    exchange = Reaction("EX_feed")
    exchange.add_metabolites({feed: -1})
    exchange.bounds = (-10, 1000)
    demand = Reaction("DM_product")
    demand.add_metabolites({product: -1})
    demand.bounds = (0, 1000)

    model.add_reactions([r2076, growth, exchange, demand])
    model.objective = demand
    return model


class CoAProtonationCurationDataTests(unittest.TestCase):
    def test_shared_core_tuples_and_identity_mappings_match_ledger_contract(self) -> None:
        """The duplicated CoA chemistry contracts must not drift independently."""

        curation = json.loads(COA_CURATION_PATH.read_text(encoding="utf-8"))
        ledger_spec = json.loads(LEDGER_SPEC_PATH.read_text(encoding="utf-8"))
        curation_groups = {group["id"]: group for group in curation["groups"]}
        convention = ledger_spec["source_chemical_convention"]
        ledger_tuples = {
            "coenzyme_a": convention["same_compartment_coa"],
            **{
                f"{chain_id}_coa": entry
                for chain_id, entry in convention["acyl_coas"].items()
            },
        }
        self.assertEqual(
            set(ledger_tuples),
            {
                "coenzyme_a",
                "lauroyl_coa",
                "myristoyl_coa",
                "palmitoyl_coa",
                "palmitoleoyl_coa",
                "stearoyl_coa",
                "oleoyl_coa",
                "linoleoyl_coa",
            },
        )
        for group_id, ledger_tuple in ledger_tuples.items():
            with self.subTest(group=group_id):
                curation_tuple = curation_groups[group_id]["target_tuple"]
                self.assertEqual(set(curation_tuple), {"formula", "charge"})
                for field in ("formula", "charge"):
                    self.assertIs(type(curation_tuple[field]), type(ledger_tuple[field]))
                    self.assertEqual(curation_tuple[field], ledger_tuple[field])

        for chain_id in ("palmitoleoyl", "linoleoyl"):
            group_id = f"{chain_id}_coa"
            with self.subTest(annotation_group=group_id):
                mapping = convention["acyl_coas"][chain_id]["canonical_mapping"]
                self.assertEqual(
                    curation_groups[group_id]["annotation_target"],
                    {
                        "bigg.metabolite": mapping["bigg"],
                        "kegg.compound": mapping["kegg"],
                        "metanetx.chemical": mapping["mnxm"],
                        "chebi": mapping["chebi"],
                    },
                )

    def test_production_manifest_is_locked_to_the_post_annotation_model_snapshot(self) -> None:
        curation = load_coa_protonation_curation()
        self.assertEqual(curation["source_model"], "../iyali26_gem/model.xml")
        self.assertEqual(curation["source_stage"], "current_main_model_xml")
        self.assertEqual(
            curation["source_sha256"],
            hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest(),
        )
        self.assertRegex(curation["source_model_fingerprint"], r"^[0-9a-f]{64}$")
        self.assertTrue({"R1419", "R1504", "R1521", "R1708", "R2004", "R613", "R614"}.issubset(curation["target_balance_reaction_ids"]))

    @unittest.skipUnless(os.environ.get("IYALI26_SOURCE_MODEL"), "requires explicit IYALI26_SOURCE_MODEL")
    def test_r613_identity_is_resolved_but_its_full_copy_frontier_stays_blocked(self) -> None:
        curation = load_coa_protonation_curation()
        blocker = next(row for row in curation["activation_blockers"] if row["reaction_id"] == "R613")
        resolution = blocker["identity_resolution"]
        frontier = resolution["frontier"]
        self.assertEqual(resolution["status"], "identity_resolved_frontier_unclosed")
        self.assertFalse(resolution["ready_for_activation"])
        self.assertEqual(resolution["biological_identity"]["ph_7_3_tuple"], {"formula": "C12H20NO4S2", "charge": -1})
        self.assertEqual(resolution["namespace_resolution"]["excluded_s6_metanetx_ids"], ["MNXM735016", "MNXM735017"])

        model = read_sbml_model(str(MODEL_PATH))
        source_fingerprint = coa_patches._model_snapshot_fingerprint(model)

        def annotation_values(metabolite, key: str) -> set[str]:
            value = metabolite.annotation.get(key)
            return {value} if isinstance(value, str) else set(value or [])

        same_identity_ids = {
            metabolite.id
            for metabolite in model.metabolites
            if coa_patches._coa_identity_key(metabolite.name) == "s(8)-succinyldihydrolipoamide"
            or annotation_values(metabolite, "bigg.metabolite") & {"HC00695", "sdhlam"}
            or annotation_values(metabolite, "metanetx.chemical") & {"MNXM735015", "MNXM735016", "MNXM735017", "MNXM735018"}
        }
        self.assertEqual(same_identity_ids, {"m853[C_mi]"})
        m853 = model.metabolites.get_by_id("m853[C_mi]")
        self.assertEqual((m853.formula, m853.charge), ("C12H21NO4S2", 0))
        self.assertEqual(annotation_values(m853, "bigg.metabolite"), {"HC00695"})
        self.assertIn("MNXM735015", annotation_values(m853, "metanetx.chemical"))
        self.assertTrue(
            annotation_values(m853, "metanetx.chemical").isdisjoint(
                resolution["namespace_resolution"]["excluded_s6_metanetx_ids"]
            )
        )
        self.assertEqual({reaction.id for reaction in m853.reactions}, set(resolution["current_snapshot"]["incident_reaction_ids"]))
        self.assertEqual(model.reactions.get_by_id("R613").check_mass_balance(), {})
        self.assertEqual(model.reactions.get_by_id("R614").check_mass_balance(), {})

        two_oxoglutarate_ids = set(frontier["required_copy_ids"])
        self.assertEqual(
            {metabolite.id for metabolite in model.metabolites if coa_patches._coa_identity_key(metabolite.name) == "2-oxoglutarate"},
            two_oxoglutarate_ids,
        )
        self.assertEqual(
            len({reaction.id for metabolite_id in two_oxoglutarate_ids for reaction in model.metabolites.get_by_id(metabolite_id).reactions}),
            frontier["incident_reaction_count"],
        )

        candidate = model.copy()
        coa_patches._apply_coa_protonation_curation(candidate, curation)
        balances_before = coa_patches._all_mass_balances(candidate)
        candidate_m853 = candidate.metabolites.get_by_id("m853[C_mi]")
        candidate_m853.formula = resolution["biological_identity"]["ph_7_3_tuple"]["formula"]
        candidate_m853.charge = resolution["biological_identity"]["ph_7_3_tuple"]["charge"]
        for metabolite_id in two_oxoglutarate_ids:
            metabolite = candidate.metabolites.get_by_id(metabolite_id)
            metabolite.formula = frontier["target_tuple"]["formula"]
            metabolite.charge = frontier["target_tuple"]["charge"]
        proton_contract = resolution["proton_contract"]
        candidate.reactions.get_by_id("R614").add_metabolites({
            candidate.metabolites.get_by_id(proton_contract["R614"]["metabolite_id"]): proton_contract["R614"]["target_coefficient"]
        })
        balances_after = coa_patches._all_mass_balances(candidate)
        self.assertEqual(
            {reaction_id for reaction_id in balances_before if balances_before[reaction_id] == {} and balances_after[reaction_id] not in ({}, None)},
            set(frontier["incremental_newly_unbalanced_reaction_ids"]),
        )
        self.assertEqual(
            {reaction_id for reaction_id in balances_before if balances_before[reaction_id] not in ({}, None) and balances_after[reaction_id] == {}},
            set(frontier["incremental_repaired_reaction_ids"]),
        )
        self.assertEqual(balances_after["R613"], {})
        self.assertEqual(balances_after["R614"], {})
        self.assertEqual(candidate.reactions.get_by_id("R613").metabolites.get(candidate.metabolites.get_by_id("m28[C_mi]"), 0), 0)
        self.assertEqual(candidate.reactions.get_by_id("R614").metabolites[candidate.metabolites.get_by_id("m28[C_mi]")], -1)
        self.assertNotIn("m853[C_mi]", {metabolite_id for group in curation["groups"] for metabolite_id in group["expected_ids"]})
        self.assertNotIn("R614", {correction["reaction_id"] for correction in curation["reaction_corrections"]})
        self.assertEqual(coa_patches._model_snapshot_fingerprint(model), source_fingerprint)

    def test_loader_rejects_a_drifted_source_snapshot_digest(self) -> None:
        curation = json.loads(COA_CURATION_PATH.read_text(encoding="utf-8"))
        curation["source_sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "curation.json"
            path.write_text(json.dumps(curation), encoding="utf-8")
            with self.assertRaises(CoAProtonationCurationError):
                load_coa_protonation_curation(path)

    def test_explicit_source_model_path_overrides_manifest_location(self) -> None:
        curation = {
            "source_model": "missing.xml",
            "source_sha256": hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest(),
        }
        self.assertEqual(
            coa_patches._validate_curation_source_file(curation, MODEL_PATH),
            MODEL_PATH.resolve(),
        )
        with self.assertRaisesRegex(CoAProtonationCurationError, "unavailable"):
            coa_patches._validate_curation_source_file(curation)

    def test_loader_rejects_r1521_handoff_tuple_drift(self) -> None:
        curation = json.loads(COA_CURATION_PATH.read_text(encoding="utf-8"))
        group_by_id(curation, "nad_plus")["target_tuple"]["charge"] = 0
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "curation.json"
            path.write_text(json.dumps(curation), encoding="utf-8")
            with self.assertRaisesRegex(CoAProtonationCurationError, "R1521 tuple contract"):
                load_coa_protonation_curation(path)

    def test_core_45_tuples_and_identity_annotation_targets_are_explicit(self) -> None:
        curation = load_coa_protonation_curation()
        observed_copies = 0
        for group_id, (copies, formula, charge) in CORE_TARGETS.items():
            group = group_by_id(curation, group_id)
            self.assertEqual(group["expected_copy_count"], copies)
            self.assertEqual(len(group["expected_ids"]), copies)
            self.assertEqual(
                group["target_tuple"], {"formula": formula, "charge": charge}
            )
            observed_copies += copies
        self.assertEqual(observed_copies, 45)

        linoleoyl = group_by_id(curation, "linoleoyl_coa")
        self.assertEqual(
            linoleoyl["annotation_target"],
            {
                "bigg.metabolite": "lnlccoa",
                "kegg.compound": "C02050",
                "metanetx.chemical": "MNXM638",
                "chebi": "CHEBI:57383",
            },
        )
        palmitoleoyl = group_by_id(curation, "palmitoleoyl_coa")
        self.assertEqual(
            palmitoleoyl["annotation_target"],
            {
                "bigg.metabolite": "hdcoa",
                "kegg.compound": "C21072",
                "metanetx.chemical": "MNXM781",
                "chebi": "CHEBI:61540",
            },
        )
        for group in (linoleoyl, palmitoleoyl):
            self.assertIn("inchi", group["replace_annotation_keys"])
            self.assertIn("inchikey", group["replace_annotation_keys"])

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_all_expected_copies_match_exact_legacy_or_target_states(self) -> None:
        model = read_sbml_model(str(MODEL_PATH))
        validation = validate_coa_protonation_curation(
            model, source_model_path=MODEL_PATH
        )
        self.assertTrue(validation["valid"], validation["errors"])

    def test_r2076_is_the_exact_enumerated_correction(self) -> None:
        curation = load_coa_protonation_curation()
        corrections = {entry["reaction_id"]: entry for entry in curation["reaction_corrections"]}
        self.assertEqual(corrections["R2076"]["metabolite_id"], "m10[C_cy]")
        self.assertEqual(corrections["R2076"]["legacy_coefficient"], 3)
        self.assertEqual(corrections["R2076"]["target_coefficient"], -1)
        self.assertEqual(set(corrections), {"R74", "R2076", "R1412"})
        self.assertEqual(
            (corrections["R1412"]["metabolite_id"], corrections["R1412"]["legacy_coefficient"], corrections["R1412"]["target_coefficient"]),
            ("m377[C_em]", 0, -1),
        )

    def test_c24_c26_closure_members_are_exact_and_source_tuples_cover_every_copy(self) -> None:
        curation = load_coa_protonation_curation()
        tetracosanoyl = group_by_id(curation, "tetracosanoyl_coa")
        oxohexacosanoyl = group_by_id(curation, "3_oxohexacosanoyl_coa")
        self.assertEqual(
            set(tetracosanoyl["expected_ids"]),
            {"m186[C_pe]", "m394[C_er]", "m1374[C_cy]", "m1436[C_em]", "m1502[C_en]", "m1511[C_lp]"},
        )
        self.assertEqual(tetracosanoyl["target_tuple"], {"formula": "C45H78N7O17P3S", "charge": -4})
        self.assertEqual(
            oxohexacosanoyl["expected_ids"], ["m1437[C_em]", "m185[C_pe]"]
        )
        self.assertEqual(oxohexacosanoyl["target_tuple"], {"formula": "C47H80N7O18P3S", "charge": -4})
        for group, chebi in (
            (tetracosanoyl, "CHEBI:65052"),
            (oxohexacosanoyl, "CHEBI:73980"),
        ):
            self.assertNotIn("source_tuple_fingerprint", group)
            self.assertEqual(
                {metabolite_id for legacy in group["legacy_tuples"] for metabolite_id in legacy["ids"]},
                set(group["expected_ids"]),
            )
            self.assertEqual(group["annotation_target"], {"chebi": chebi})
            self.assertEqual(
                set(group["replace_annotation_keys"]),
                coa_patches._COA_CHEMICAL_IDENTITY_KEYS,
            )

    def test_legacy_tuple_group_cannot_omit_a_copy(self) -> None:
        curation = json.loads(COA_CURATION_PATH.read_text(encoding="utf-8"))
        group_by_id(curation, "3_oxohexacosanoyl_coa")["legacy_tuples"][0]["ids"].pop()
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "curation.json"
            path.write_text(json.dumps(curation), encoding="utf-8")
            with self.assertRaises(CoAProtonationCurationError):
                load_coa_protonation_curation(path)

    def test_redox_and_focal_acid_groups_are_pinned(self) -> None:
        curation = load_coa_protonation_curation()
        for group_id, copies, target in (
            ("nad_plus", 6, ("C21H26N7O14P2", -1)),
            ("nadh", 5, ("C21H27N7O14P2", -2)),
            ("nadp_plus", 6, ("C21H25N7O17P3", -3)),
            ("nadph", 6, ("C21H26N7O17P3", -4)),
            ("acetate", 6, ("C2H3O2", -1)),
            ("succinate", 3, ("C4H4O4", -2)),
            ("pyruvate", 3, ("C3H3O3", -1)),
            ("oxaloacetate", 4, ("C4H2O5", -2)),
        ):
            group = group_by_id(curation, group_id)
            self.assertEqual(group["expected_copy_count"], copies)
            self.assertEqual((group["target_tuple"]["formula"], group["target_tuple"]["charge"]), target)
        hydroxy = group_by_id(curation, "er_vlcfa_3r_hydroxyhexacosanoyl_coa")
        self.assertEqual(hydroxy["expected_ids"], ["m1446[C_em]", "m1546[C_pe]"])
        self.assertIn({"ids": ["m1546[C_pe]"], "formula": None, "charge": 0}, hydroxy["legacy_tuples"])

    def test_blocked_curation_is_not_wired_into_the_build_pipeline(self) -> None:
        self.assertNotIn("normalize_coa_protonation", inspect.getsource(apply_all_patches))
        main_source = (REPOSITORY / "scripts" / "gem_annotate" / "main.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("normalize_coa_protonation", main_source)


class CoAProtonationGateTests(unittest.TestCase):
    def test_approved_curation_requires_a_target_model_fingerprint(self) -> None:
        curation = approved_minimal_curation()
        curation.pop("target_model_fingerprint")
        with self.assertRaises(CoAProtonationCurationError):
            audit_coa_protonation_curation(build_balanced_fixture(), curation)

    def test_public_gate_rejects_curation_without_source_provenance(self) -> None:
        curation = approved_minimal_curation()
        curation.pop("source_stage")
        curation.pop("source_model_fingerprint")
        with self.assertRaises(CoAProtonationCurationError):
            audit_coa_protonation_curation(build_balanced_fixture(), curation)

    def test_fractional_boolean_and_nonfinite_live_charges_are_rejected(self) -> None:
        curation = approved_minimal_curation()
        for label, charge in (
            ("fractional", -0.5),
            ("boolean", True),
            ("nonfinite", float("inf")),
        ):
            with self.subTest(label=label):
                model = build_balanced_fixture()
                model.metabolites.get_by_id("coa_c").charge = charge
                with self.assertRaises(CoAProtonationCurationError):
                    validate_coa_protonation_curation(model, curation)

    def test_unknown_prestate_missing_copy_and_new_copy_fail_closed(self) -> None:
        curation = approved_minimal_curation()

        unknown = build_balanced_fixture()
        unknown.metabolites.get_by_id("coa_c").formula = "C21H31N7O16P3S"
        errors = validate_coa_protonation_curation(unknown, curation)["errors"]
        self.assertTrue(any("third tuple" in error for error in errors), errors)

        missing_curation = copy.deepcopy(curation)
        missing_curation["groups"][0]["expected_ids"] = ["missing_coa"]
        missing_curation["groups"][0]["legacy_tuples"][0]["ids"] = ["missing_coa"]
        errors = validate_coa_protonation_curation(build_balanced_fixture(), missing_curation)["errors"]
        self.assertTrue(any("missing expected" in error for error in errors), errors)

        extra = build_balanced_fixture()
        extra.add_metabolites(
            [
                Metabolite(
                    "coa_extra",
                    name="coenzyme A_C21H36N7O16P3S",
                    formula="C21H36N7O16P3S",
                    charge=0,
                    compartment="c",
                )
            ]
        )
        errors = validate_coa_protonation_curation(extra, curation)["errors"]
        self.assertTrue(any("unexpected same-identity" in error for error in errors), errors)

    def test_unavailable_objective_fails_closed(self) -> None:
        curation = approved_minimal_curation()
        for label, replacement in (
            ("none", {"return_value": None}),
            ("exception", {"side_effect": RuntimeError("simulated solver failure")}),
        ):
            with self.subTest(label=label):
                with patch.object(Model, "slim_optimize", **replacement):
                    report = audit_coa_protonation_curation(
                        build_balanced_fixture(), curation
                    )
                self.assertFalse(report["objective_unchanged"])
                self.assertFalse(report["gate_passed"])
                self.assertFalse(report["ready_for_activation"])

    def test_r2076_zero_newly_unbalanced_gate_idempotence_and_fba(self) -> None:
        model = build_balanced_fixture()
        curation = approved_minimal_curation()
        before = model.slim_optimize()
        self.assertEqual(model.reactions.get_by_id("R2076").check_mass_balance(), {})

        counts = normalize_coa_protonation(model, curation)
        self.assertEqual(
            counts,
            {
                "metabolites": 1,
                "annotations": 0,
                "reaction_corrections": 1,
                "identity_corrections": 0,
            },
        )
        self.assertEqual(model.metabolites.get_by_id("coa_c").formula, "C21H32N7O16P3S")
        self.assertEqual(model.metabolites.get_by_id("coa_c").charge, -4)
        self.assertEqual(model.reactions.get_by_id("R2076").metabolites[model.metabolites.get_by_id("h_c")], -1)
        self.assertEqual(model.reactions.get_by_id("R2076").check_mass_balance(), {})
        self.assertAlmostEqual(model.slim_optimize(), before)

        self.assertEqual(
            normalize_coa_protonation(model, curation),
            {
                "metabolites": 0,
                "annotations": 0,
                "reaction_corrections": 0,
                "identity_corrections": 0,
            },
        )

    def test_production_curation_is_blocked_before_any_mutation(self) -> None:
        model = build_balanced_fixture()
        curation = approved_minimal_curation()
        curation["activation_state"] = "blocked"
        before = (
            model.metabolites.get_by_id("coa_c").formula,
            model.metabolites.get_by_id("coa_c").charge,
            model.reactions.get_by_id("R2076").metabolites[model.metabolites.get_by_id("h_c")],
        )
        with self.assertRaises(CoAProtonationActivationBlocked):
            normalize_coa_protonation(model, curation)
        after = (
            model.metabolites.get_by_id("coa_c").formula,
            model.metabolites.get_by_id("coa_c").charge,
            model.reactions.get_by_id("R2076").metabolites[model.metabolites.get_by_id("h_c")],
        )
        self.assertEqual(after, before)

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_preapplied_target_model_cannot_bypass_blockers_or_absolute_balance(self) -> None:
        """A target-state model cannot turn known closure failures into a green gate."""

        model = read_sbml_model(str(MODEL_PATH))
        curation = load_coa_protonation_curation()
        coa_patches._apply_coa_protonation_curation(model, curation)
        approved = copy.deepcopy(curation)
        approved["activation_state"] = "approved"
        approved["target_model_fingerprint"] = coa_patches._model_snapshot_fingerprint(
            model
        )

        report = audit_coa_protonation_curation(
            model, approved, source_model_path=MODEL_PATH
        )
        self.assertFalse(report["activation_blockers_clear"])
        self.assertTrue(report["source_snapshot"]["verified"])
        self.assertTrue(report["source_snapshot"]["target_model_fingerprint_verified"])
        self.assertEqual(set(report["target_reactions_unbalanced"]), {"R613"})
        self.assertFalse(report["gate_passed"])
        self.assertFalse(report["ready_for_activation"])

        # Even an attempted deletion of the declared blockers cannot hide the
        # target-model residuals from the absolute-balance safety gate.
        without_blockers = copy.deepcopy(approved)
        without_blockers["activation_blockers"] = []
        report = audit_coa_protonation_curation(
            model, without_blockers, source_model_path=MODEL_PATH
        )
        self.assertTrue(report["activation_blockers_clear"])
        self.assertEqual(set(report["target_reactions_unbalanced"]), {"R613"})
        self.assertFalse(report["gate_passed"])
        self.assertFalse(report["ready_for_activation"])

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_er_vlcfa_3r_tuple_closes_r1426_without_a_proton_repair(self) -> None:
        """The newly curated ER tuple resolves R1426, not by adding H+."""

        model = read_sbml_model(str(MODEL_PATH))
        curation = load_coa_protonation_curation()
        coa_patches._apply_coa_protonation_curation(model, curation)

        for metabolite_id in ("m1446[C_em]", "m1546[C_pe]"):
            target = model.metabolites.get_by_id(metabolite_id)
            self.assertEqual(target.name, "(3R)-3-hydroxyhexacosanoyl-CoA")
            self.assertEqual((target.formula, target.charge), ("C47H82N7O18P3S", -4))
            self.assertEqual(target.annotation["chebi"], "CHEBI:76378")
        self.assertEqual(model.reactions.get_by_id("R1504").check_mass_balance(), {})
        self.assertEqual(model.reactions.get_by_id("R1521").check_mass_balance(), {})
        self.assertEqual(model.reactions.get_by_id("R1419").check_mass_balance(), {})
        self.assertEqual(model.reactions.get_by_id("R1426").check_mass_balance(), {})
        self.assertEqual(
            model.reactions.get_by_id("R1419").name,
            "B-ketoacyl-CoA reductase ((3R)-3-hydroxyhexacosanoyl-CoA)",
        )
        self.assertEqual(model.reactions.get_by_id("R1426").annotation["rhea"], ["39211", "39214"])

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_c24_c26_candidate_closes_authoritative_reaction_set_and_is_idempotent(self) -> None:
        model = read_sbml_model(str(MODEL_PATH))
        curation = load_coa_protonation_curation()
        coa_patches._apply_coa_protonation_curation(model, curation)

        for reaction_id in ("R1412", "R1419", "R1504", "R1521", "R1708", "R2004", "R80", "R202", "R204"):
            self.assertEqual(model.reactions.get_by_id(reaction_id).check_mass_balance(), {})
        for group_id, chebi in (
            ("tetracosanoyl_coa", "CHEBI:65052"),
            ("3_oxohexacosanoyl_coa", "CHEBI:73980"),
        ):
            for metabolite_id in group_by_id(curation, group_id)["expected_ids"]:
                annotation = model.metabolites.get_by_id(metabolite_id).annotation
                self.assertEqual(annotation.get("chebi"), chebi)
                self.assertTrue(
                    all(key == "chebi" or key not in annotation for key in coa_patches._COA_CHEMICAL_IDENTITY_KEYS)
                )
        self.assertEqual(
            model.reactions.get_by_id("R1412").annotation["rhea"], ["36515", "36516"]
        )
        self.assertEqual(
            model.reactions.get_by_id("R1419").annotation["rhea"], ["39207", "39208"]
        )
        self.assertEqual(
            model.reactions.get_by_id("R80").annotation["rhea"], ["78639", "78641"]
        )
        self.assertEqual(
            coa_patches._apply_coa_protonation_curation(model, curation),
            {"metabolites": 0, "annotations": 0, "reaction_corrections": 0, "identity_corrections": 0},
        )

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_candidate_transaction_rolls_back_if_a_prewrite_contract_fails(self) -> None:
        model = read_sbml_model(str(MODEL_PATH))
        curation = load_coa_protonation_curation()
        before = coa_patches._model_snapshot_fingerprint(model)
        with patch.object(coa_patches, "_reaction_matches_contract", return_value=False):
            with self.assertRaises(CoAProtonationCurationError):
                coa_patches._apply_coa_protonation_curation(model, curation)
        self.assertEqual(coa_patches._model_snapshot_fingerprint(model), before)

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_coa_tuple_curation_does_not_mark_ledger_source_normalized(self) -> None:
        """Tuple fixes alone leave the ledger's independent annotation gate open."""

        from scripts.lipid_moiety_ledger import (
            POOL_REACTION_ID,
            _source_chemical_verification,
            _source_pool_bindings,
            load_spec,
        )

        model = read_sbml_model(str(MODEL_PATH))
        curation = load_coa_protonation_curation()
        coa_patches._apply_coa_protonation_curation(model, curation)
        convention = load_spec(LEDGER_SPEC_PATH)["source_chemical_convention"]
        report = _source_chemical_verification(
            _source_pool_bindings(model, convention, POOL_REACTION_ID), convention
        )
        self.assertEqual(report["status"], "not_normalized")
        self.assertEqual(len(report["errors"]), 24)

    @unittest.skipUnless(MODEL_PATH.exists(), "requires repository model.xml")
    def test_production_audit_reports_the_recursive_closure_blockers(self) -> None:
        report = audit_coa_protonation_curation(
            read_sbml_model(str(MODEL_PATH)), source_model_path=MODEL_PATH
        )
        self.assertEqual(report["activation_state"], "blocked")
        self.assertFalse(report["ready_for_activation"])
        self.assertFalse(report["gate_passed"])
        self.assertGreater(len(report["newly_unbalanced"]), 100)
        self.assertTrue({"R344", "R613", "R742"}.issubset(report["newly_unbalanced"]))
        self.assertEqual(set(report["target_reactions_unbalanced"]), {"R613"})
        self.assertEqual(report["non_h_charge_deltas"], {})
        self.assertTrue(report["source_snapshot"]["verified"])
        self.assertFalse(report["activation_blockers_clear"])
        self.assertTrue(report["objective_unchanged"])


if __name__ == "__main__":
    unittest.main()
