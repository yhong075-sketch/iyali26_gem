"""Schema-v2 contract tests for the read-only aggregate acyl-moiety ledger."""

from __future__ import annotations

import copy
import csv
import hashlib
import importlib
import inspect
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
SPEC_PATH = REPOSITORY / "data" / "lipid_moiety_ledger_spec.json"
CURATION_PATH = REPOSITORY / "data" / "lipid_combo_curation.csv"
DOC_PATH = REPOSITORY / "docs" / "lipid_moiety_ledger.md"
CORE_PATH = REPOSITORY / "scripts" / "lipid_moiety_ledger.py"
PLANNER_PATH = REPOSITORY / "scripts" / "plan_lipid_moiety_ledger.py"
MODEL_PATH = REPOSITORY / "model.xml"

if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))


EXPECTED_TUPLES = {
    "same_compartment_coa": {
        "chemical_identity": "coenzyme A",
        "formula": "C21H32N7O16P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "coa",
            "kegg": "C00010",
            "mnxm": "MNXM12",
            "chebi": "CHEBI:57287",
        },
    },
    "lauroyl": {
        "chemical_identity": "lauroyl-CoA",
        "formula": "C33H54N7O17P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "ddcacoa",
            "kegg": "C01832",
            "mnxm": "MNXM363",
            "chebi": "CHEBI:57375",
        },
    },
    "myristoyl": {
        "chemical_identity": "myristoyl-CoA",
        "formula": "C35H58N7O17P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "tdcoa",
            "kegg": "C02593",
            "mnxm": "MNXM224",
            "chebi": "CHEBI:57385",
        },
    },
    "palmitoyl": {
        "chemical_identity": "palmitoyl-CoA",
        "formula": "C37H62N7O17P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "pmtcoa",
            "kegg": "C00154",
            "mnxm": "MNXM88",
            "chebi": "CHEBI:57379",
        },
    },
    "palmitoleoyl": {
        "chemical_identity": "palmitoleoyl-CoA",
        "formula": "C37H60N7O17P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "hdcoa",
            "kegg": "C21072",
            "mnxm": "MNXM781",
            "chebi": "CHEBI:61540",
        },
    },
    "stearoyl": {
        "chemical_identity": "stearoyl-CoA",
        "formula": "C39H66N7O17P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "stcoa",
            "kegg": "C00412",
            "mnxm": "MNXM272",
            "chebi": "CHEBI:57394",
        },
    },
    "oleoyl": {
        "chemical_identity": "oleoyl-CoA",
        "formula": "C39H64N7O17P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "ocdce9coa",
            "kegg": "C00510",
            "mnxm": "MNXM686",
            "chebi": "CHEBI:57387",
        },
    },
    "linoleoyl": {
        "chemical_identity": "linoleoyl-CoA",
        "formula": "C39H62N7O17P3S",
        "charge": -4,
        "canonical_mapping": {
            "bigg": "lnlccoa",
            "kegg": "C02050",
            "mnxm": "MNXM638",
            "chebi": "CHEBI:57383",
        },
    },
}


def load_spec() -> dict[str, Any]:
    return json.loads(SPEC_PATH.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def immutable_input_paths() -> dict[str, Path]:
    return {
        "model.xml": MODEL_PATH,
        "data/lipid_combo_curation.csv": CURATION_PATH,
        "data/lipid_moiety_ledger_spec.json": SPEC_PATH,
    }


def immutable_input_sha256() -> dict[str, str]:
    return {name: sha256(path) for name, path in immutable_input_paths().items()}


def parse_residue_formula(formula: str) -> tuple[int, int, int]:
    match = re.fullmatch(r"C(\d+)H(\d+)O(\d*)", formula)
    if not match:
        raise AssertionError(f"not an explicit C/H/O residue formula: {formula!r}")
    carbon, hydrogen, oxygen = match.groups()
    return int(carbon), int(hydrogen), int(oxygen or 1)


def parse_canonical_manifest(payload: bytes) -> dict[str, Any]:
    assert payload.endswith(b"\n")
    decoded = json.loads(payload.decode("utf-8"))
    expected = (
        json.dumps(decoded, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
        + b"\n"
    )
    assert payload == expected
    return decoded


def canonical_annotation(mapping: Mapping[str, str]) -> dict[str, object]:
    """Exercise scalar and list annotation forms accepted by schema v2."""

    return {
        "bigg.metabolite": mapping["bigg"],
        "kegg.compound": [mapping["kegg"]],
        "metanetx.chemical": mapping["mnxm"],
        "chebi": [mapping["chebi"]],
    }


def normalized_source_fixture(core: Any) -> tuple[Any, Mapping[str, Any]]:
    """Normalize only an in-memory source copy; never repair model.xml."""

    spec = core.load_spec(SPEC_PATH)
    convention = spec["source_chemical_convention"]
    source = core._read_model(MODEL_PATH).copy()
    bindings = core._source_pool_bindings(source, convention, core.POOL_REACTION_ID)
    coa = bindings[0][3]
    coa_expected = convention["same_compartment_coa"]
    coa.name = coa_expected["chemical_identity"]
    coa.formula = coa_expected["formula"]
    coa.charge = coa_expected["charge"]
    coa.annotation = canonical_annotation(coa_expected["canonical_mapping"])
    for chain_id, _member, metabolite, _coa in bindings:
        expected = convention["acyl_coas"][chain_id]
        metabolite.name = expected["chemical_identity"]
        metabolite.formula = expected["formula"]
        metabolite.charge = expected["charge"]
        metabolite.annotation = canonical_annotation(expected["canonical_mapping"])
    return source, convention


def run_cli(*arguments: str, hash_seed: str | None = None) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    if hash_seed is not None:
        environment["PYTHONHASHSEED"] = hash_seed
    return subprocess.run(
        [sys.executable, "-B", str(PLANNER_PATH), *arguments],
        cwd=REPOSITORY,
        text=True,
        capture_output=True,
        check=False,
        env=environment,
    )


class LipidMoietyLedgerDataContractTests(unittest.TestCase):
    def test_schema_v2_is_read_only_and_nonactivating(self) -> None:
        spec = load_spec()
        self.assertEqual(spec["schema_version"], 2)
        self.assertFalse(spec["activation_ready"])
        self.assertEqual(spec["execution_mode"]["build_scope"], "in_memory_only")
        self.assertTrue(spec["execution_mode"]["source_model_immutable"])
        self.assertTrue(spec["execution_mode"]["xml_output_forbidden"])
        self.assertTrue(spec["execution_mode"]["model_activation_forbidden"])
        self.assertTrue(spec["execution_mode"]["immutable_input_alias_output_forbidden"])
        self.assertEqual(spec["public_compile_contract"]["entrypoint"], "compile_dry_run")
        self.assertIn(
            "activation_opt_in", spec["public_compile_contract"]["forbidden_parameters"]
        )

    def test_biochemical_ph_anionic_tuples_and_canonical_mappings_are_complete(self) -> None:
        convention = load_spec()["source_chemical_convention"]
        self.assertEqual(convention["name"], "biochemical_pH_anionic")
        self.assertTrue(convention["annotation_policy"]["wrong_same_namespace_forbidden"])
        self.assertEqual(
            convention["annotation_policy"]["required_canonical_fields"],
            ["bigg", "kegg", "mnxm", "chebi"],
        )
        self.assertEqual(convention["same_compartment_coa"], EXPECTED_TUPLES["same_compartment_coa"])
        self.assertEqual(set(convention["acyl_coas"]), set(EXPECTED_TUPLES) - {"same_compartment_coa"})
        for chain_id, expected in EXPECTED_TUPLES.items():
            if chain_id == "same_compartment_coa":
                continue
            with self.subTest(chain=chain_id):
                self.assertEqual(convention["acyl_coas"][chain_id], expected)

    def test_residues_curated_tag_budget_and_architecture_are_closed(self) -> None:
        spec = load_spec()
        chains = spec["chains"]
        self.assertEqual(len(chains), 7)
        self.assertEqual({chain["id"] for chain in chains}, set(EXPECTED_TUPLES) - {"same_compartment_coa"})
        for chain in chains:
            carbon, hydrogen, oxygen = parse_residue_formula(chain["acyl_residue_formula"])
            self.assertEqual(carbon, chain["carbon_atoms"])
            self.assertEqual(oxygen, 1)
            self.assertEqual(hydrogen, 2 * carbon - 2 * chain["double_bonds"] - 2)
            self.assertEqual(chain["curation_member"], f"{chain['id']}-CoA")
        with CURATION_PATH.open(newline="", encoding="utf-8") as handle:
            kept = [
                row
                for row in csv.DictReader(handle)
                if row["layer"] == "tri" and row["verdict"] == "keep"
            ]
        self.assertEqual(len(kept), 35)
        self.assertEqual(
            sum(len(set(row["member_chains"].split(";"))) for row in kept), 81
        )
        architecture = spec["architecture"]
        self.assertEqual(architecture["first_acylation_reactions"], 7)
        self.assertEqual(architecture["second_acylation_reactions"], 7)
        self.assertEqual(architecture["terminal_tag_routes"], 81)
        self.assertEqual(architecture["biochemical_reactions_exact"], 95)
        self.assertEqual(architecture["generated_metabolites_exact"], 44)
        self.assertEqual(architecture["generated_metabolites_maximum"], 44)

    def test_manifest_contract_and_remaining_activation_blockers_are_explicit(self) -> None:
        spec = load_spec()
        manifest = spec["required_manifest"]
        self.assertIn("source_chemical_convention", manifest["must_include"])
        self.assertIn("source_chemical_convention", manifest["validation_required_fields"])
        self.assertIn(
            "source_coa_acyl_formula_charge_verified", manifest["validation_flat_fields"]
        )
        blocker_ids = {blocker["id"] for blocker in spec["required_blockers"]}
        self.assertEqual(
            blocker_ids,
            {
                "read_only_manifest_only_prototype",
                "candidate_gpr_compartment_evidence_unresolved",
                "aggregate_ledger_lacks_molecular_species_provenance",
                "chemistry_source_not_normalized",
            },
        )
        self.assertNotIn("non_uniform_source_acyl_coa_charge_convention", blocker_ids)

    def test_capability_boundary_retains_probabilities_only_for_explicit_soft_prior(self) -> None:
        boundary = load_spec()["capability_boundary"]
        self.assertEqual(
            boundary["curation_inputs_retained"],
            ["tri", "keep", "member_chains", "prob", "cumulative_coverage", "xPOOL abundance weights"],
        )
        self.assertEqual(boundary["curation_inputs_not_retained"], [])
        self.assertEqual(boundary["toy_supply_upper_bound"], 3)
        self.assertIn(
            "quantitative TAG composition inference without an explicit soft-prior probe", boundary["unsupported_uses"]
        )
        self.assertIn("quantitative TAG yield inference", boundary["unsupported_uses"])

    def test_documentation_states_chemical_and_capability_boundaries(self) -> None:
        text = DOC_PATH.read_text(encoding="utf-8").lower()
        for required in (
            "biochemical-ph anionic",
            "source_coa_acyl_formula_charge_verified",
            "chemistry_source_not_normalized",
            "aggregate chain histogram",
            "dag",
            "sn position",
            "cross-compartment",
            "gpr",
            "activation_ready: false",
            "not an integrated regression",
            "hard-link output",
            "atomically",
        ):
            with self.subTest(required=required):
                self.assertIn(required, text)


@unittest.skipUnless(CORE_PATH.exists(), "ledger module is unavailable")
class LipidMoietyLedgerFixtureV2Tests(unittest.TestCase):
    """Strict source checks against a canonicalized in-memory fixture only."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.core = importlib.import_module("scripts.lipid_moiety_ledger")
        cls.source, cls.convention = normalized_source_fixture(cls.core)
        cls.bindings = cls.core._source_pool_bindings(
            cls.source, cls.convention, cls.core.POOL_REACTION_ID
        )
        cls.source_report = cls.core._source_chemical_verification(
            cls.bindings, cls.convention
        )
        cls.chains = cls.core.load_chains(cls.source, cls.convention)
        cls.curated_tags = cls.core.load_curated_tag_records(CURATION_PATH, cls.chains)
        cls.tags = tuple(item["composition"] for item in cls.curated_tags)
        cls.core._validate_spec_against_plan(
            cls.core.load_spec(SPEC_PATH), cls.chains, cls.tags
        )
        (
            cls.model,
            cls.classes,
            cls.transitions,
            cls.ledger_ids,
            cls.tag_records,
            cls.source_proxies,
        ) = cls.core._build_ledger(cls.chains, cls.curated_tags)
        cls.validation = {
            "mass_balance": cls.core._validate_mass_balance(cls.model, cls.transitions),
            "occupancy": cls.core._validate_occupancy(cls.transitions, cls.ledger_ids),
            "chain_identity": cls.core._validate_chain_identity(
                cls.transitions, cls.chains, cls.ledger_ids, cls.tag_records
            ),
            "ids": cls.core._validate_ids(cls.model, cls.transitions),
            "budget": cls.core._validate_budget(
                cls.classes, cls.transitions, cls.tag_records
            ),
            "fba_probe": cls.core._run_fba_probes(
                cls.model, cls.chains, cls.ledger_ids, cls.tag_records
            ),
            "abundance": cls.core._abundance_report(
                cls.tag_records,
                cls.core._source_pool_abundance_weights(cls.source, cls.chains),
            ),
            "source_chemical_convention": cls.source_report,
            "independent_source_baseline": {"passed": True},
            "source_model_baseline": {"passed": True},
            "source_model_unchanged": True,
            "curation_unchanged": True,
            "spec_unchanged": True,
            "determinism": {"passed": True, "status": "fixture"},
        }
        cls.manifest = cls.core._manifest(
            input_sha256={"fixture": "0" * 64},
            chains=cls.chains,
            classes=cls.classes,
            transitions=cls.transitions,
            ledger_ids=cls.ledger_ids,
            tags=cls.tag_records,
            source_proxies=cls.source_proxies,
            validation=cls.validation,
        )
        cls.manifest["gates"]["determinism"] = True
        cls.core._update_activation_status(cls.manifest)

    def test_canonical_identifier_normalization_is_exact(self) -> None:
        self.assertEqual(
            self.core._normalize_annotation_value("chebi", "CHEBI:57287"),
            "CHEBI:57287",
        )
        self.assertEqual(
            self.core._normalize_annotation_value("chebi", "57287"),
            "CHEBI:57287",
        )

    def test_strict_loader_accepts_exact_tuples_and_scalar_or_list_annotations(self) -> None:
        self.assertTrue(self.source_report["passed"])
        self.assertEqual(self.source_report["status"], "verified")
        self.assertEqual(self.source_report["tuple_origin"], "source_model")
        self.assertEqual(len(self.chains), 7)
        self.assertEqual(self.core.load_chains(self.source), self.chains)
        self.assertTrue(all(chain.coa_formula == "C21H32N7O16P3S" for chain in self.chains))
        self.assertTrue(all(chain.coa_charge == -4 for chain in self.chains))
        self.assertTrue(all(chain.acyl_coa_charge == -4 for chain in self.chains))
        self.assertTrue(all(chain.residue_charge == 0 for chain in self.chains))
        for chain in self.chains:
            with self.subTest(chain=chain.id):
                expected = EXPECTED_TUPLES[chain.id]
                self.assertEqual(chain.acyl_coa_formula, expected["formula"])
                self.assertEqual(chain.residue_formula, next(
                    item["acyl_residue_formula"]
                    for item in load_spec()["chains"]
                    if item["id"] == chain.id
                ))

    def test_fixture_reactions_are_charge_balanced_and_all_generated_lipids_are_neutral(self) -> None:
        self.assertEqual(len(self.transitions), 95)
        self.assertEqual(self.validation["mass_balance"]["unbalanced_reactions"], {})
        for transition in self.transitions:
            with self.subTest(reaction=transition.id):
                reaction = self.model.reactions.get_by_id(transition.id)
                self.assertEqual(reaction.check_mass_balance(), {})
        self.assertTrue(all(item.charge == 0 for item in self.classes if item.generated))
        self.assertTrue(all(tag["charge"] == 0 for tag in self.tag_records))
        self.assertTrue(self.validation["fba_probe"]["passed"])

    def test_candidate_gprs_cover_all_generated_reactions_without_clearing_blocker(self) -> None:
        self.assertFalse(self.model.genes)
        candidate = self.core.build_candidate_gpr_evaluation_model(self.model)
        rules = [
            candidate.reactions.get_by_id(transition.id).gene_reaction_rule
            for transition in self.transitions
        ]
        self.assertEqual(len(rules), 95)
        self.assertTrue(all(rules))
        self.assertNotIn("YALI1F13550g", " ".join(rules))
        self.assertNotIn("YALI1D21511g", " ".join(rules))
        self.assertIn("YALI1C00230g", rules)
        self.assertIn("YALI1E38810g", rules)
        self.assertIn(
            "candidate_gpr_compartment_evidence_unresolved",
            {item["id"] for item in self.manifest["blockers"]},
        )

    def test_verified_chemistry_gate_replaces_the_old_charge_blocker(self) -> None:
        self.assertTrue(self.manifest["gates"]["source_coa_acyl_formula_charge_verified"])
        self.assertTrue(self.manifest["validation"]["source_coa_acyl_formula_charge_verified"])
        self.assertTrue(self.manifest["gates"]["all_required_gates"])
        self.assertFalse(self.manifest["activation_ready"])
        blocker_ids = {item["id"] for item in self.manifest["blockers"]}
        self.assertNotIn("chemistry_source_not_normalized", blocker_ids)
        self.assertNotIn("non_uniform_source_acyl_coa_charge_convention", blocker_ids)
        self.assertEqual(
            blocker_ids,
            {
                "read_only_manifest_only_prototype",
                "candidate_gpr_compartment_evidence_unresolved",
                "aggregate_ledger_lacks_molecular_species_provenance",
            },
        )

    def _mutated_fixture(self) -> tuple[Any, Mapping[str, Any]]:
        # COBRA copies can retain annotation dictionaries by reference, so make
        # a fresh normalized fixture before every destructive in-memory test.
        return normalized_source_fixture(self.core)

    def test_formula_and_charge_drift_fail_closed(self) -> None:
        for chain_id, field, value in (
            ("lauroyl", "formula", "C33H55N7O17P3S"),
            ("myristoyl", "charge", -3),
        ):
            source, convention = self._mutated_fixture()
            binding = next(
                item
                for item in self.core._source_pool_bindings(
                    source, convention, self.core.POOL_REACTION_ID
                )
                if item[0] == chain_id
            )
            setattr(binding[2], field, value)
            with self.subTest(chain=chain_id, field=field):
                with self.assertRaises(self.core.LedgerError):
                    self.core.load_chains(source, convention)

    def test_coa_and_identity_annotation_drift_fail_closed(self) -> None:
        cases = (
            ("coa_formula", "same_compartment_coa", "formula", "C21H36N7O16P3S"),
            ("palmitoleoyl_chebi", "palmitoleoyl", "chebi", "CHEBI:52381"),
            ("linoleoyl_bigg", "linoleoyl", "bigg.metabolite", "stcoa"),
        )
        for label, chain_id, field, value in cases:
            source, convention = self._mutated_fixture()
            bindings = self.core._source_pool_bindings(
                source, convention, self.core.POOL_REACTION_ID
            )
            metabolite = bindings[0][3] if chain_id == "same_compartment_coa" else next(
                binding[2] for binding in bindings if binding[0] == chain_id
            )
            if field == "formula":
                metabolite.formula = value
            else:
                metabolite.annotation[field] = [value, value]
            with self.subTest(case=label):
                with self.assertRaises(self.core.LedgerError):
                    self.core.load_chains(source, convention)

    def test_extra_wrong_value_in_a_canonical_namespace_is_rejected(self) -> None:
        source, convention = self._mutated_fixture()
        linoleoyl = next(
            binding[2]
            for binding in self.core._source_pool_bindings(
                source, convention, self.core.POOL_REACTION_ID
            )
            if binding[0] == "linoleoyl"
        )
        linoleoyl.annotation["bigg.metabolite"] = ["lnlccoa", "stcoa"]
        with self.assertRaises(self.core.LedgerError):
            self.core.load_chains(source, convention)


@unittest.skipUnless(CORE_PATH.exists() and PLANNER_PATH.exists(), "ledger modules unavailable")
class LipidMoietyLedgerIntegrationTests(unittest.TestCase):
    """Public compiler contracts, including the temporary source-mismatch path."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.core = importlib.import_module("scripts.lipid_moiety_ledger")
        cls.model_sha_before = sha256(MODEL_PATH)
        cls.result = cls.core.compile_dry_run(
            MODEL_PATH, CURATION_PATH, spec_path=SPEC_PATH
        )
        cls.payload = cls.result.canonical_json
        cls.manifest = parse_canonical_manifest(cls.payload)
        cls.model_sha_after = sha256(MODEL_PATH)

    def test_public_compile_contract_and_inputs_remain_immutable(self) -> None:
        self.assertEqual(self.model_sha_before, self.model_sha_after)
        self.assertEqual(
            self.manifest["input_sha256"],
            {
                "model.xml": sha256(MODEL_PATH),
                "data/lipid_combo_curation.csv": sha256(CURATION_PATH),
                "data/lipid_moiety_ledger_spec.json": sha256(SPEC_PATH),
            },
        )
        for name in ("compile_dry_run", "build_dry_run", "plan_lipid_moiety_ledger"):
            with self.subTest(api=name):
                self.assertNotIn(
                    "activation_opt_in", inspect.signature(getattr(self.core, name)).parameters
                )
        validation = self.manifest["validation"]
        self.assertTrue(validation["source_model_unchanged"])
        self.assertTrue(validation["curation_unchanged"])
        self.assertTrue(validation["spec_unchanged"])

    def test_manifest_preserves_architecture_neutral_lipids_and_95_balanced_reactions(self) -> None:
        self.assertEqual(self.manifest["schema_version"], 2)
        self.assertEqual(self.manifest["counts"]["chains"], 7)
        self.assertEqual(self.manifest["counts"]["tag_outputs"], 35)
        self.assertEqual(self.manifest["counts"]["terminal_tag_routes"], 81)
        self.assertEqual(self.manifest["counts"]["biochemical_reactions"], 95)
        self.assertEqual(self.manifest["counts"]["generated_metabolites"], 44)
        self.assertTrue(self.manifest["validation"]["mass_balance"]["passed"])
        self.assertTrue(all(chain["residue_charge"] == 0 for chain in self.manifest["chains"]))
        self.assertTrue(all(tag["charge"] == 0 for tag in self.manifest["tag_outputs"]))
        for transition in self.result.transitions:
            with self.subTest(reaction=transition.id):
                self.assertEqual(
                    self.result.model.reactions.get_by_id(transition.id).check_mass_balance(),
                    {},
                )

    def test_probability_metadata_and_opt_in_soft_prior_are_well_formed(self) -> None:
        abundance = self.manifest["abundance"]
        self.assertAlmostEqual(abundance["kept_probability_coverage"], 0.99083121, places=7)
        self.assertAlmostEqual(abundance["normalized_probability_sum"], 1.0, places=8)
        self.assertEqual(abundance["solver_use"], "opt_in_soft_prior_only")
        self.assertEqual(abundance["cumulative_coverage_use"], "curation_qc_only")
        self.assertTrue(all(tag["prob"] > 0 for tag in self.manifest["tag_outputs"]))
        self.assertAlmostEqual(
            sum(tag["normalized_prob"] for tag in self.manifest["tag_outputs"]), 1.0, places=8
        )
        probe = self.core.run_composition_soft_prior(
            self.result.model, self.manifest["tag_outputs"]
        )
        self.assertEqual(probe["fit_status"], "optimal")
        self.assertGreater(probe["primary_tag_maximum"], 0.0)
        self.assertLess(probe["l1_deviation"], 1e-7)

    def test_manifest_reports_candidate_only_gpr_coverage(self) -> None:
        candidate = self.manifest["candidate_gpr"]
        self.assertEqual(candidate["reaction_count"], 95)
        self.assertEqual(candidate["assigned_reaction_count"], 95)
        self.assertEqual(candidate["activation_status"], "blocked_pending_compartment_and_substrate_evidence")
        self.assertTrue(all(item["candidate_gpr"] for item in self.manifest["generated_reactions"]))

    @unittest.skipUnless(
        importlib.util.find_spec("memote") is not None,
        "memote is an optional targeted validation dependency",
    )
    def test_targeted_memote_candidate_evaluation_mass_charge_and_gpr_checks(self) -> None:
        from memote.suite.tests.test_basic import test_gene_protein_reaction_rule_presence
        from memote.suite.tests.test_consistency import (
            test_reaction_charge_balance,
            test_reaction_mass_balance,
        )

        candidate = self.core.build_candidate_gpr_evaluation_model(self.result.model)
        test_reaction_mass_balance(candidate)
        test_reaction_charge_balance(candidate)
        test_gene_protein_reaction_rule_presence(candidate)

    def test_source_chemistry_is_never_silently_repaired(self) -> None:
        report = self.manifest["source_chemical_convention"]
        blocker_ids = {item["id"] for item in self.manifest["blockers"]}
        self.assertEqual(
            report["passed"],
            self.manifest["gates"]["source_coa_acyl_formula_charge_verified"],
        )
        if report["passed"]:
            self.assertEqual(report["status"], "verified")
            self.assertEqual(report["tuple_origin"], "source_model")
            self.assertNotIn("chemistry_source_not_normalized", blocker_ids)
        else:
            self.assertEqual(report["status"], "not_normalized")
            self.assertEqual(report["tuple_origin"], "explicit_canonical_curation_tuple")
            self.assertIn("chemistry_source_not_normalized", blocker_ids)
            self.assertTrue(report["errors"])
        self.assertFalse(self.manifest["activation_ready"])

    def test_nonchemical_activation_blockers_and_capability_boundary_remain_open(self) -> None:
        blocker_ids = {item["id"] for item in self.manifest["blockers"]}
        self.assertTrue(
            {
                "read_only_manifest_only_prototype",
                "candidate_gpr_compartment_evidence_unresolved",
                "aggregate_ledger_lacks_molecular_species_provenance",
            }.issubset(blocker_ids)
        )
        scopes = {item["id"] for item in self.manifest["capability_loss"]["blocked_scopes"]}
        self.assertEqual(
            scopes,
            {
                "unsupported_transport_scope",
                "unsupported_remodeling_scope",
                "unsupported_topology_scope",
            },
        )
        self.assertTrue(self.manifest["gates"]["activation_blockers_open"])
        self.assertFalse(self.manifest["activation_ready"])

    def test_fba_pfba_supply_cut_and_determinism_contracts_hold(self) -> None:
        validation = self.manifest["validation"]
        self.assertTrue(validation["fba_probe"]["passed"])
        self.assertTrue(validation["source_model_baseline"]["passed"])
        self.assertEqual(
            validation["source_model_baseline"]["objective_reaction"], "biomass_C"
        )
        self.assertEqual(
            validation["source_model_baseline"]["objective_direction"], "max"
        )
        self.assertEqual(
            validation["source_model_baseline"]["objective_expression"],
            {"biomass_C": 1.0},
        )
        self.assertEqual(
            validation["source_model_baseline"]["scope"],
            "independent_source_baseline_not_integrated_regression",
        )
        self.assertTrue(validation["determinism"]["passed"])
        self.assertEqual(len(validation["fba_probe"]["positive_controls"]), 35)
        self.assertEqual(len(validation["fba_probe"]["required_chain_cuts"]), 81)
        self.assertEqual(
            validation["fba_probe"]["all_acyl_closed_internal_flux_blocked"][
                "reaction_count"
            ],
            95,
        )
        second = self.core.compile_dry_run(MODEL_PATH, CURATION_PATH, spec_path=SPEC_PATH)
        self.assertEqual(second.canonical_json, self.payload)

    def test_strict_real_source_loader_is_kept_for_the_canonical_model(self) -> None:
        report = self.manifest["source_chemical_convention"]
        if not report["passed"]:
            self.skipTest("canonical source tuples/annotations have not yet been rebuilt")
        source = self.core._read_model(MODEL_PATH)
        convention = self.core.load_spec(SPEC_PATH)["source_chemical_convention"]
        chains = self.core.load_chains(source, convention)
        self.assertEqual(len(chains), 7)

    def test_source_baseline_rejects_objective_drift_and_tiny_secondary_terms(self) -> None:
        alternatives = [
            reaction.id
            for reaction in self.core._read_model(MODEL_PATH).reactions
            if reaction.id != "biomass_C"
        ]
        self.assertTrue(alternatives)
        alternate = alternatives[0]

        wrong_objective = self.core._read_model(MODEL_PATH)
        wrong_objective.objective = alternate
        with self.assertRaises(self.core.LedgerError):
            self.core._independent_source_baseline(wrong_objective)

        tiny_secondary = self.core._read_model(MODEL_PATH)
        tiny_secondary.objective = "biomass_C"
        tiny_secondary.reactions.get_by_id(alternate).objective_coefficient = 1e-9
        with self.assertRaises(self.core.LedgerError):
            self.core._independent_source_baseline(tiny_secondary)

    def test_cli_hash_seed_determinism_and_atomic_output_alias_guards(self) -> None:
        before = immutable_input_sha256()
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temporary_directory:
            directory = Path(temporary_directory)
            payloads: list[bytes] = []
            for seed in ("1", "777"):
                output = directory / f"seed-{seed}.json"
                completed = run_cli("--output", str(output), hash_seed=seed)
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                payloads.append(output.read_bytes())
                parse_canonical_manifest(payloads[-1])
            self.assertEqual(payloads[0], payloads[1])

            payload = b'{"atomic_test":true}\n'
            symlink = directory / "model-symlink.json"
            os.symlink(MODEL_PATH, symlink)
            with self.assertRaises(self.core.LedgerError):
                self.core.atomic_write_manifest(
                    symlink, payload, immutable_input_paths(), before
                )
            output = directory / "regular-output.json"
            self.core.atomic_write_manifest(output, payload, immutable_input_paths(), before)
            self.assertEqual(output.read_bytes(), payload)
        self.assertEqual(immutable_input_sha256(), before)

    def test_strict_spec_drift_fails_closed(self) -> None:
        base = load_spec()
        mutations = (
            lambda spec: spec.__setitem__("schema_version", 1),
            lambda spec: spec["source_chemical_convention"]["same_compartment_coa"].__setitem__(
                "charge", 0
            ),
            lambda spec: spec["source_chemical_convention"]["acyl_coas"]["linoleoyl"][
                "canonical_mapping"
            ].__setitem__("bigg", "stcoa"),
            lambda spec: spec["required_blockers"].pop(),
            lambda spec: spec["validation_gates"].__setitem__(
                "source_coa_acyl_formula_charge_verified", "not checked"
            ),
        )
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temporary_directory:
            directory = Path(temporary_directory)
            for index, mutate in enumerate(mutations):
                candidate = copy.deepcopy(base)
                mutate(candidate)
                path = directory / f"drift-{index}.json"
                path.write_text(json.dumps(candidate), encoding="utf-8")
                with self.subTest(index=index):
                    with self.assertRaises(self.core.LedgerError):
                        self.core.load_spec(path)

    def test_schema_v2_rejects_exact_type_confusion(self) -> None:
        """JSON numbers/bools must not exploit Python equality coercion."""

        base = load_spec()
        mutations = {
            "schema_version_float": lambda spec: spec.__setitem__(
                "schema_version", 2.0
            ),
            "execution_bool_as_int": lambda spec: spec["execution_mode"].__setitem__(
                "source_model_immutable", 1
            ),
            "declared_charge_float": lambda spec: spec[
                "source_chemical_convention"
            ]["same_compartment_coa"].__setitem__("charge", -4.0),
            "reaction_count_float": lambda spec: spec["architecture"].__setitem__(
                "biochemical_reactions_exact", 95.0
            ),
            "nested_chain_count_float": lambda spec: spec["validation_gates"][
                "chain_identity_conserved"
            ].__setitem__("chains", 7.0),
        }
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temporary_directory:
            directory = Path(temporary_directory)
            for label, mutate in mutations.items():
                candidate = copy.deepcopy(base)
                mutate(candidate)
                path = directory / f"type-confusion-{label}.json"
                path.write_text(json.dumps(candidate), encoding="utf-8")
                with self.subTest(label=label):
                    with self.assertRaises(self.core.LedgerError):
                        self.core.load_spec(path)

    def test_load_spec_binds_chain_ids_and_residues_to_declared_acyl_coas(self) -> None:
        base = load_spec()
        mutations = {
            "renamed_self_consistent_decanoyl": lambda chain: chain.update(
                {
                    "id": "decanoyl",
                    "curation_member": "decanoyl-CoA",
                    "acyl_residue_formula": "C10H18O",
                    "carbon_atoms": 10,
                    "double_bonds": 0,
                }
            ),
            "lauroyl_with_self_consistent_wrong_residue": lambda chain: chain.update(
                {
                    "acyl_residue_formula": "C10H18O",
                    "carbon_atoms": 10,
                    "double_bonds": 0,
                }
            ),
        }
        with tempfile.TemporaryDirectory(dir="/private/tmp") as temporary_directory:
            directory = Path(temporary_directory)
            for label, mutate in mutations.items():
                candidate = copy.deepcopy(base)
                lauroyl = next(
                    chain for chain in candidate["chains"] if chain["id"] == "lauroyl"
                )
                mutate(lauroyl)
                path = directory / f"chain-contract-{label}.json"
                path.write_text(json.dumps(candidate), encoding="utf-8")
                with self.subTest(label=label):
                    with self.assertRaises(self.core.LedgerError):
                        self.core.load_spec(path)


if __name__ == "__main__":
    unittest.main()
