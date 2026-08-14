"""Regression tests for fail-closed ER VLCFA (3R) curation."""

from __future__ import annotations

import copy
import importlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra.io import read_sbml_model

from scripts.gem_annotate.vlcfa_stereochemistry import (
    ERVLCFAStereochemistryError,
    correct_er_vlcfa_3r_stereochemistry,
    load_er_vlcfa_stereochemistry_curation,
    verify_er_vlcfa_3r_stereochemistry_target,
)
from scripts.gem_annotate import vlcfa_stereochemistry


REPOSITORY = Path(__file__).resolve().parents[1]
RAW_MODEL = REPOSITORY / "data" / "iyli21.xml"
CURATION = REPOSITORY / "data" / "er_vlcfa_3r_stereochemistry.json"


@unittest.skipUnless(RAW_MODEL.exists(), "requires source SBML")
class ERVLCFAStereochemistryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.model = read_sbml_model(str(RAW_MODEL))
        self.curation = load_er_vlcfa_stereochemistry_curation()

    @staticmethod
    def _target_state(model):
        metabolite = model.metabolites.get_by_id("m1446[C_em]")
        reactions = [model.reactions.get_by_id(reaction_id) for reaction_id in ("R1419", "R1426")]
        return (
            metabolite.name,
            metabolite.formula,
            metabolite.charge,
            copy.deepcopy(metabolite.annotation),
            tuple((reaction.name, copy.deepcopy(reaction.annotation)) for reaction in reactions),
        )

    def test_raw_source_is_corrected_offline_and_peroxisomal_copy_is_untouched(self) -> None:
        peroxisomal = self.model.metabolites.get_by_id("m1546[C_pe]")
        untouched = (
            peroxisomal.name,
            peroxisomal.formula,
            peroxisomal.charge,
            copy.deepcopy(peroxisomal.annotation),
        )

        changes = correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)

        target = self.model.metabolites.get_by_id("m1446[C_em]")
        self.assertEqual(changes, {"metabolites": 2, "reactions": 1, "annotations": 3})
        self.assertEqual(target.name, "(3R)-3-hydroxyhexacosanoyl-CoA")
        self.assertEqual((target.formula, target.charge), ("C47H86N7O18P3S", 0))
        self.assertEqual(target.annotation["chebi"], "CHEBI:76465")
        self.assertEqual(
            self.model.reactions.get_by_id("R1419").name,
            "B-ketoacyl-CoA reductase ((3R)-3-hydroxyhexacosanoyl-CoA)",
        )
        self.assertEqual(self.model.reactions.get_by_id("R1419").annotation["rhea"], ["39207", "39208"])
        self.assertEqual(self.model.reactions.get_by_id("R1426").annotation["rhea"], ["39211", "39214"])
        self.assertEqual(
            (
                peroxisomal.name,
                peroxisomal.formula,
                peroxisomal.charge,
                copy.deepcopy(peroxisomal.annotation),
            ),
            untouched,
        )

    def test_correction_is_idempotent(self) -> None:
        correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)
        self.assertEqual(
            correct_er_vlcfa_3r_stereochemistry(self.model, self.curation),
            {"metabolites": 0, "reactions": 0, "annotations": 0},
        )

    def test_loader_rejects_source_sha_and_target_contract_drift(self) -> None:
        original = json.loads(CURATION.read_text(encoding="utf-8"))
        for field, value in (
            ("source_sha256", "0" * 64),
            ("target_contract_sha256", "0" * 64),
        ):
            with self.subTest(field=field), tempfile.TemporaryDirectory(dir="/private/tmp") as temp_dir:
                changed = copy.deepcopy(original)
                changed[field] = value
                path = Path(temp_dir) / "curation.json"
                path.write_text(json.dumps(changed), encoding="utf-8")
                with self.assertRaises(ERVLCFAStereochemistryError):
                    load_er_vlcfa_stereochemistry_curation(path)

    def test_direct_curation_cannot_bypass_source_sha_and_leaves_model_unchanged(self) -> None:
        changed = copy.deepcopy(self.curation)
        changed["source_sha256"] = "0" * 64
        before = self._target_state(self.model)

        with self.assertRaises(ERVLCFAStereochemistryError):
            correct_er_vlcfa_3r_stereochemistry(self.model, changed)

        self.assertEqual(self._target_state(self.model), before)

    def test_direct_curation_rechecks_the_actual_source_file_sha(self) -> None:
        before = self._target_state(self.model)

        with patch.object(vlcfa_stereochemistry, "_sha256_file", return_value="0" * 64):
            with self.assertRaises(ERVLCFAStereochemistryError):
                correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)

        self.assertEqual(self._target_state(self.model), before)

    def test_recomputed_curation_digest_cannot_authorize_an_impossible_target_tuple(self) -> None:
        changed = copy.deepcopy(self.curation)
        changed["neutral_patch"]["metabolite"]["target_tuple"] = {"formula": "C1", "charge": 99}
        before = self._target_state(self.model)

        # The checked-in SHA is the authority.  Even if a caller supplies a
        # freshly computed digest, no field inside JSON can replace that anchor.
        self.assertNotEqual(
            vlcfa_stereochemistry._target_contract_digest(changed),
            vlcfa_stereochemistry._TRUSTED_CURATION_SHA256,
        )
        with self.assertRaises(ERVLCFAStereochemistryError):
            correct_er_vlcfa_3r_stereochemistry(self.model, changed)
        self.assertEqual(self._target_state(self.model), before)

        # Defense in depth: the deep tuple contract independently rejects C1/+99
        # if the digest verifier itself is intentionally mocked in a test.
        with patch.object(
            vlcfa_stereochemistry,
            "_TRUSTED_CURATION_SHA256",
            vlcfa_stereochemistry._target_contract_digest(changed),
        ):
            with self.assertRaises(ERVLCFAStereochemistryError):
                correct_er_vlcfa_3r_stereochemistry(self.model, changed)
        self.assertEqual(self._target_state(self.model), before)

    def test_float_schema_version_is_rejected_even_with_a_matching_mock_digest(self) -> None:
        changed = copy.deepcopy(self.curation)
        changed["schema_version"] = 2.0
        before = self._target_state(self.model)

        with patch.object(
            vlcfa_stereochemistry,
            "_TRUSTED_CURATION_SHA256",
            vlcfa_stereochemistry._target_contract_digest(changed),
        ):
            with self.assertRaises(ERVLCFAStereochemistryError):
                correct_er_vlcfa_3r_stereochemistry(self.model, changed)

        self.assertEqual(self._target_state(self.model), before)

    def test_missing_excluded_peroxisomal_copy_fails_closed(self) -> None:
        peroxisomal = self.model.metabolites.get_by_id("m1546[C_pe]")
        self.model.remove_metabolites([peroxisomal])

        with self.assertRaises(ERVLCFAStereochemistryError):
            correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)

    def test_malformed_or_unknown_curation_fails_before_mutation(self) -> None:
        for label, mutate in (
            (
                "malformed annotation target",
                lambda value: value["neutral_patch"]["metabolite"].__setitem__("annotation_target", {"chebi": ["CHEBI:76465"]}),
            ),
            ("unknown key", lambda value: value.__setitem__("unreviewed_override", True)),
        ):
            with self.subTest(label=label):
                changed = copy.deepcopy(self.curation)
                mutate(changed)
                before = self._target_state(self.model)
                with self.assertRaises(ERVLCFAStereochemistryError):
                    correct_er_vlcfa_3r_stereochemistry(self.model, changed)
                self.assertEqual(self._target_state(self.model), before)
                # This makes the test exercise strict schema validation rather
                # than merely the outer immutable-digest guard.
                with patch.object(
                    vlcfa_stereochemistry,
                    "_TRUSTED_CURATION_SHA256",
                    vlcfa_stereochemistry._target_contract_digest(changed),
                ):
                    with self.assertRaises(ERVLCFAStereochemistryError):
                        correct_er_vlcfa_3r_stereochemistry(self.model, changed)
                self.assertEqual(self._target_state(self.model), before)

    def test_reaction_bounds_direction_and_stoichiometry_drift_fail_closed(self) -> None:
        for label, mutate in (
            ("bounds", lambda m: setattr(m.reactions.get_by_id("R1419"), "upper_bound", 999.0)),
            ("direction", lambda m: setattr(m.reactions.get_by_id("R1426"), "lower_bound", 0.0)),
            ("stoichiometry", lambda m: m.reactions.get_by_id("R1426").add_metabolites({m.metabolites.get_by_id("m1447[C_em]"): 1.0})),
        ):
            with self.subTest(label=label):
                model = read_sbml_model(str(RAW_MODEL))
                mutate(model)
                with self.assertRaises(ERVLCFAStereochemistryError):
                    correct_er_vlcfa_3r_stereochemistry(model, self.curation)

    def test_main_is_offline_and_legacy_entry_forwards_to_unified_builder(self) -> None:
        main_source = (REPOSITORY / "scripts" / "gem_annotate" / "main.py").read_text(encoding="utf-8")
        call = main_source.index("n_vlcfa_stereo = correct_er_vlcfa_3r_stereochemistry(model)")
        offline_branch = main_source.index("MetaNetX files not found")
        self.assertGreater(call, offline_branch)

        update_model = importlib.import_module("scripts.update_model")
        with patch("scripts.gem_annotate.main.main", return_value="unified") as unified:
            self.assertEqual(update_model.legacy_main(), "unified")
            self.assertEqual(update_model.main(), "unified")
        self.assertEqual(unified.call_count, 2)

        final_gate = main_source.index("verify_er_vlcfa_3r_stereochemistry_target(model)")
        writer = main_source.index("write_sbml_model(model, str(OUTPUT_MODEL_PATH))")
        self.assertLess(final_gate, writer)

    def test_direct_script_entry_resolves_without_running_the_builder(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/update_model.py", "--resolve-entry-smoke"],
            cwd=REPOSITORY,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("update_model entry resolution: OK", completed.stdout)

    def test_final_target_gate_detects_postcondition_drift_without_mutating_model(self) -> None:
        correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)
        self.model.metabolites.get_by_id("m1446[C_em]").annotation["chebi"] = "CHEBI:bad"
        before = self._target_state(self.model)

        with self.assertRaises(ERVLCFAStereochemistryError):
            verify_er_vlcfa_3r_stereochemistry_target(self.model)

        self.assertEqual(self._target_state(self.model), before)

    def test_final_target_gate_rejects_uncurated_live_closure_tuple(self) -> None:
        correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)
        partner = self.model.metabolites.get_by_id("m1437[C_em]")
        partner.formula, partner.charge = "C1", 99
        before = (partner.formula, partner.charge)

        with self.assertRaises(ERVLCFAStereochemistryError):
            verify_er_vlcfa_3r_stereochemistry_target(self.model)

        self.assertEqual((partner.formula, partner.charge), before)

    def test_final_target_gate_rejects_reinserted_chemical_identity_annotation(self) -> None:
        correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)
        target = self.model.metabolites.get_by_id("m1446[C_em]")
        target.annotation["inchikey"] = "OLD-IDENTITY"
        before = self._target_state(self.model)

        with self.assertRaises(ERVLCFAStereochemistryError):
            verify_er_vlcfa_3r_stereochemistry_target(self.model)

        self.assertEqual(self._target_state(self.model), before)

    def test_final_target_gate_rejects_excluded_peroxisomal_identity_drift(self) -> None:
        correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)
        peroxisomal = self.model.metabolites.get_by_id("m1546[C_pe]")
        peroxisomal.name = "UNREVIEWED"

        with self.assertRaises(ERVLCFAStereochemistryError):
            verify_er_vlcfa_3r_stereochemistry_target(self.model)

    def test_postcondition_failure_rolls_back_every_live_edit(self) -> None:
        before = self._target_state(self.model)
        with patch.object(
            vlcfa_stereochemistry,
            "_verify_live_target_state",
            side_effect=ERVLCFAStereochemistryError("forced final-gate failure"),
        ):
            with self.assertRaises(ERVLCFAStereochemistryError):
                correct_er_vlcfa_3r_stereochemistry(self.model, self.curation)
        self.assertEqual(self._target_state(self.model), before)

    def test_pH_tuple_is_recorded_but_not_part_of_the_neutral_patch(self) -> None:
        pH_identity = self.curation["biochemical_ph_identity"]
        self.assertEqual(pH_identity["chebi"], "CHEBI:76378")
        self.assertEqual(pH_identity["biochemical_ph_tuple"], {"formula": "C47H82N7O18P3S", "charge": -4})
        self.assertEqual(self.curation["scope"]["excluded_metabolite_ids"], ["m1546[C_pe]"])


if __name__ == "__main__":
    unittest.main()
