#!/usr/bin/env python3
"""Audit the current-main R1521 handoff without proposing a model patch."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import tempfile
import warnings
from math import isclose, isfinite
from pathlib import Path

from cobra.core.formula import elements_and_molecular_weights
from cobra.io import read_sbml_model

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.gem_annotate.patches import (
    CoAProtonationActivationBlocked,
    _model_snapshot_fingerprint,
    _safe_mass_balance,
    load_coa_protonation_curation,
)


REPOSITORY = Path(__file__).resolve().parents[1]
HANDOFF_PATH = REPOSITORY / "data" / "r1521_current_snapshot_handoff.json"
ER_EVIDENCE_PATH = REPOSITORY / "data" / "er_vlcfa_3r_stereochemistry.json"
COA_CURATION_PATH = REPOSITORY / "data" / "coa_protonation_curation.json"
_TRUSTED_CONTRACT_SHA256 = "18d60ef4e583311fd95581dcce2a96f4af6ac6c1d0dc19bc36799c163334d678"
_BALANCE_EXPECTATION_KEYS = frozenset(
    {
        "migrated_metabolite_count", "incident_reaction_count",
        "evaluable_reaction_count", "unassessable_reaction_count",
        "changed_residual_count", "newly_unbalanced_count", "repaired_count",
        "still_unbalanced_changed_count", "balanced_unchanged_count",
        "unbalanced_unchanged_count", "migrated_metabolite_ids_sha256",
        "incident_reaction_ids_sha256", "changed_reaction_ids_sha256",
        "newly_unbalanced_ids_sha256", "repaired_ids_sha256",
        "still_unbalanced_changed_ids_sha256", "unassessable_ids_sha256",
        "ledger_sha256", "unassessable_reaction_ids", "repaired_reaction_ids",
    }
)
_IDENTITY_KEYS = frozenset(
    {
        "chebi", "hmdb", "inchi", "inchikey", "lipidmaps", "lipidmapsm",
        "metacyc.compound", "metacycm", "metanetx.chemical", "seed.compound", "seedm",
    }
)
_R1521_IDENTITY_NAMES = frozenset(
    {"(r)-3-hydroxyhexacosanoyl-coa", "(3r)-3-hydroxyhexacosanoyl-coa"}
)
_COMPLETE_FORMULA = re.compile(r"(?:[A-Z][a-z]?\d*)+\Z")


class R1521CurrentSnapshotError(RuntimeError):
    """The blocked R1521 handoff is malformed or no longer matches current main."""


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _contract_digest(handoff: dict) -> str:
    contract = dict(handoff)
    contract.pop("target_contract_sha256", None)
    return hashlib.sha256(
        json.dumps(contract, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def _canonical_json_file_digest(path: Path) -> str:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise R1521CurrentSnapshotError(
            f"cannot read R1521 evidence dependency {path}: {exc}"
        ) from exc
    _require(isinstance(value, dict), f"R1521 evidence dependency is not an object: {path}")
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise R1521CurrentSnapshotError(message)


def _is_digest(value) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _validate_handoff(handoff: dict) -> None:
    _require(isinstance(handoff, dict), "R1521 handoff must be an object")
    _require(
        set(handoff)
        == {
            "schema_version", "curation_id", "activation", "source_model", "source_sha256",
            "source_model_fingerprint", "evidence_dependencies", "legacy_metabolite_contract", "local_counterfactual",
            "global_nad_migration_gate", "reaction_contracts",
            "partial_known_formula_residual_contract",
            "frontier_closure_gate", "blockers", "target_contract_sha256",
        },
        "unknown or missing R1521 handoff key",
    )
    _require(type(handoff["schema_version"]) is int and handoff["schema_version"] == 2, "unsupported R1521 handoff schema")
    _require(handoff["curation_id"] == "iyali26_r1521_current_snapshot_handoff", "R1521 handoff identity drifted")
    _require(
        handoff["activation"]
        == {
            "state": "blocked",
            "ready_for_activation": False,
            "production_apply_forbidden": True,
            "approval_rule": handoff["activation"].get("approval_rule"),
        }
        and bool(handoff["activation"]["approval_rule"]),
        "R1521 handoff must remain explicitly blocked",
    )
    _require(handoff["source_model"] == "../iyali26_gem/model.xml", "unexpected current-main source path")
    _require(_is_digest(handoff["source_sha256"]), "invalid R1521 source SHA")
    _require(_is_digest(handoff["source_model_fingerprint"]), "invalid R1521 source fingerprint")
    _require(handoff["target_contract_sha256"] == _TRUSTED_CONTRACT_SHA256, "R1521 contract digest drifted")
    _require(_contract_digest(handoff) == _TRUSTED_CONTRACT_SHA256, "R1521 contract content drifted")
    _require(
        handoff["evidence_dependencies"]
        == {
            "er_vlcfa_3r_stereochemistry": {
                "path": "data/er_vlcfa_3r_stereochemistry.json",
                "canonical_contract_sha256": "1a527ee0bb99b1703f2cc8e18b4af147aa3bde7ca1dcb634d8220d4cb375ab80",
            },
            "r1521_rhea_frontier_source_audit": {
                "path": "data/r1521_rhea_frontier_source_audit.json",
                "canonical_contract_sha256": "ab72557b047c50e8c42af9732eadd44cdbf21cacdaeb0c0a7b5aa05e4cafa925",
            },
            "r1521_kegg_mnx_frontier_source_audit": {
                "path": "data/r1521_kegg_mnx_frontier_source_audit.json",
                "canonical_contract_sha256": "91cb0278234f2a349c4d61d843bf3d0fee440bdb2784bd2796376e43395db17c",
            },
            "r1521_unresolved_frontier_source_audit": {
                "path": "data/r1521_unresolved_frontier_source_audit.json",
                "canonical_contract_sha256": "1fc387d8f7270e0486058ba8307b279fceb87208525a79a36324f8cc29c959bf",
            },
        },
        "R1521 evidence dependency contract drifted",
    )

    metabolite = handoff["legacy_metabolite_contract"]
    _require(isinstance(metabolite, dict) and set(metabolite) == {"identity_name", "expected_copy_count", "copies"}, "invalid R1521 metabolite contract")
    _require(metabolite["identity_name"] == "(r)-3-hydroxyhexacosanoyl-coa", "R1521 metabolite identity drifted")
    _require(type(metabolite["expected_copy_count"]) is int and metabolite["expected_copy_count"] == 1, "R1521 copy count drifted")
    _require(isinstance(metabolite["copies"], list) and len(metabolite["copies"]) == 1, "R1521 copy inventory drifted")
    copy = metabolite["copies"][0]
    _require(
        isinstance(copy, dict)
        and set(copy) == {"id", "name", "legacy_tuple", "absent_annotation_keys"}
        and copy["id"] == "m1546[C_pe]"
        and copy["name"] == "(R)-3-hydroxyhexacosanoyl-CoA_"
        and copy["legacy_tuple"] == {"formula": None, "charge": 0}
        and set(copy["absent_annotation_keys"]) == _IDENTITY_KEYS
        and len(copy["absent_annotation_keys"]) == len(_IDENTITY_KEYS),
        "R1521 legacy metabolite contract drifted",
    )
    candidate = handoff["local_counterfactual"]
    _require(
        isinstance(candidate, dict)
        and set(candidate)
        == {
            "evidence", "supersedes", "rhea", "target_tuples", "gpr_evidence",
            "annotation_conflicts",
        },
        "invalid R1521 candidate transaction",
    )
    _require("MFE2" in candidate["evidence"] and "3R" in candidate["evidence"] and "78635" in candidate["evidence"] and "78638" in candidate["evidence"], "R1521 MFE2/3R evidence drifted")
    _require("er_vlcfa_3r_stereochemistry.json" in candidate["supersedes"] and "3S" in candidate["supersedes"], "R1521 supersession record drifted")
    _require(candidate["rhea"] == {"hydratase_master": "39211", "hydratase_direction": "39213", "dehydrogenase_master": "78635", "dehydrogenase_bidirectional": "78638", "dehydrogenase_written_direction": "78637"}, "R1521 Rhea directionality evidence drifted")
    expected_targets = {
        "m1546[C_pe]": ("C47H82N7O18P3S", -4),
        "m226[C_pe]": ("C47H80N7O17P3S", -4),
        "m116[C_pe]": ("C21H26N7O14P2", -1),
        "m119[C_pe]": ("C21H27N7O14P2", -2),
        "m185[C_pe]": ("C47H80N7O18P3S", -4),
        "m182[C_pe]": ("C21H32N7O16P3S", -4),
        "m183[C_pe]": ("C23H34N7O17P3S", -4),
        "m186[C_pe]": ("C45H78N7O17P3S", -4),
    }
    targets = candidate["target_tuples"]
    _require(isinstance(targets, list) and {entry.get("id") for entry in targets} == set(expected_targets), "R1521 candidate target inventory drifted")
    for entry in targets:
        _require(set(entry) in ({"id", "legacy_tuple", "formula", "charge"}, {"id", "name", "legacy_tuple", "formula", "charge", "annotation_target"}), "invalid R1521 candidate target")
        _require((entry["formula"], entry["charge"]) == expected_targets[entry["id"]], f"R1521 candidate target drifted for {entry['id']}")
        _require(entry["legacy_tuple"] == {"formula": entry["legacy_tuple"].get("formula"), "charge": entry["legacy_tuple"].get("charge")}, f"invalid R1521 legacy tuple for {entry['id']}")
    r1521_target = next(entry for entry in targets if entry["id"] == "m1546[C_pe]")
    _require(r1521_target["name"] == "(3R)-3-hydroxyhexacosanoyl-CoA" and r1521_target["annotation_target"] == {"chebi": "CHEBI:76378"}, "R1521 3R target identity drifted")
    gpr = candidate["gpr_evidence"]
    _require(
        isinstance(gpr, dict)
        and gpr.get("systematic_id") == "YALI1E18441g"
        and gpr.get("established_name") == "MFE2"
        and gpr.get("legacy_locus_tag") == "YALI0E15378g"
        and gpr.get("protein_accessions") == ["Q9P4D9", "F2Z6I5"]
        and gpr.get("model_gpr") == "YALI1E18441g"
        and gpr.get("reaction_ids") == ["R1504", "R1521"]
        and gpr.get("model_gpr_evidence_status", "").startswith(
            "model/GPR assignment only"
        )
        and isinstance(gpr.get("protein_function"), str)
        and gpr["protein_function"]
        and gpr.get("evidence_status", {}).get("exact_c26_substrate_turnover")
        == "unverified"
        and isinstance(gpr.get("sources"), list)
        and gpr["sources"],
        "R1521 MFE2 gene/GPR evidence drifted",
    )
    _require(
        [entry.get("reaction_id") for entry in candidate["annotation_conflicts"]]
        == ["R1504", "R1521"]
        and all(
            "no production annotation edit" in entry.get("decision", "").lower()
            for entry in candidate["annotation_conflicts"]
        ),
        "R1521 stale annotation evidence drifted",
    )
    global_gate = handoff["global_nad_migration_gate"]
    _require(
        isinstance(global_gate, dict)
        and set(global_gate)
        == {
            "baseline", "nad_plus", "nadh", "expected",
            "historical_partial_formula_baseline", "production_apply_forbidden",
        }
        and global_gate["baseline"]
        == "source_snapshot_only; no CoA/local transaction pre-applied"
        and global_gate["production_apply_forbidden"] is True,
        "invalid global NAD migration gate",
    )
    for label, identity_name, ids, legacy, target in (
        ("nad_plus", "nad_c21h27n7o14p2", ["m27[C_mi]", "m116[C_pe]", "m122[C_cy]", "m817[C_nu]", "m960[C_er]", "m1458[C_em]"], {"formula": "C21H27N7O14P2", "charge": 0}, {"formula": "C21H26N7O14P2", "charge": -1}),
        ("nadh", "nadh_c21h29n7o14p2", ["m30[C_mi]", "m119[C_pe]", "m123[C_cy]", "m957[C_er]", "m1457[C_em]"], {"formula": "C21H29N7O14P2", "charge": 0}, {"formula": "C21H27N7O14P2", "charge": -2}),
    ):
        entry = global_gate[label]
        _require(isinstance(entry, dict) and entry == {"identity_name": entry.get("identity_name"), "ids": entry.get("ids"), "legacy_tuple": entry.get("legacy_tuple"), "target_tuple": entry.get("target_tuple")} and entry["identity_name"] == identity_name and entry["ids"] == ids and entry["legacy_tuple"] == legacy and entry["target_tuple"] == target, f"R1521 global {label} migration drifted")
    expected = global_gate["expected"]
    _require(
        isinstance(expected, dict) and set(expected) == _BALANCE_EXPECTATION_KEYS,
        "invalid global NAD migration expectation",
    )
    _require(
        tuple(
            expected[key]
            for key in (
                "migrated_metabolite_count", "incident_reaction_count",
                "evaluable_reaction_count", "unassessable_reaction_count",
                "changed_residual_count", "newly_unbalanced_count",
                "repaired_count", "still_unbalanced_changed_count",
                "balanced_unchanged_count", "unbalanced_unchanged_count",
            )
        )
        == (11, 121, 89, 32, 83, 73, 2, 8, 6, 0)
        and len(expected["unassessable_reaction_ids"]) == 32
        and "R364" in expected["unassessable_reaction_ids"]
        and expected["repaired_reaction_ids"] == ["R1889", "R570"],
        "global NAD migration partition drifted",
    )
    historical = global_gate["historical_partial_formula_baseline"]
    _require(
        historical.get("reproduced_on_current_snapshot") is True
        and tuple(
            historical.get(key)
            for key in (
                "incident_reaction_count", "changed_residual_count",
                "newly_unbalanced_count", "repaired_count",
            )
        ) == (121, 114, 73, 2)
        and "not exact frontier-closure evidence" in historical.get("warning", ""),
        "historical NAD partial-formula baseline drifted",
    )
    _require(
        all(_is_digest(value) for key, value in expected.items() if key.endswith("sha256")),
        "global NAD migration digest drifted",
    )

    frontier = handoff["frontier_closure_gate"]
    _require(
        isinstance(frontier, dict)
        and set(frontier)
        == {
            "transaction_policy", "balance_ledger_contract", "source_audit_coverage",
            "local_counterfactual_full_model", "atomic_identity_transaction",
            "r364_symbolic_contract", "production_apply_forbidden",
        }
        and frontier["production_apply_forbidden"] is True
        and "never as a bulk repair" in frontier["transaction_policy"],
        "invalid R1521 frontier closure gate",
    )
    _require(
        frontier["balance_ledger_contract"].get("schema_version") == 1
        and "allow_nan=false"
        in frontier["balance_ledger_contract"].get("digest_algorithm", ""),
        "R1521 balance-ledger contract drifted",
    )
    coverage = frontier["source_audit_coverage"]
    _require(
        coverage.get("batch_reaction_counts")
        == {"rhea": 53, "kegg_mnx": 22, "unresolved": 21}
        and coverage.get("audited_regression_reaction_count") == 96
        and coverage.get("outside_global_incident_count") == 23
        and coverage.get("global_incident_not_source_audited_count") == 48
        and coverage.get("atomic_newly_unbalanced_total") == 146
        and coverage.get("atomic_newly_unbalanced_covered") == 95
        and coverage.get("atomic_newly_unbalanced_not_audited") == 51
        and coverage.get("all_atomic_frontier_identities_closed") is False
        and coverage.get("isolated_proton_edits_authorized") is False
        and _is_digest(coverage.get("audited_regression_reaction_ids_sha256"))
        and _is_digest(coverage.get("atomic_newly_unbalanced_not_audited_ids_sha256")),
        "R1521 frontier source-audit coverage drifted",
    )
    global_unaudited = coverage.get("global_incident_not_source_audited_reaction_ids")
    atomic_unaudited = coverage.get(
        "atomic_newly_unbalanced_not_audited_reaction_ids"
    )
    _require(
        isinstance(global_unaudited, list)
        and len(global_unaudited) == 48
        and global_unaudited == sorted(global_unaudited)
        and _ids_digest(global_unaudited)
        == coverage.get("global_incident_not_source_audited_ids_sha256")
        and isinstance(atomic_unaudited, list)
        and len(atomic_unaudited) == 51
        and atomic_unaudited == sorted(atomic_unaudited)
        and _ids_digest(atomic_unaudited)
        == coverage.get("atomic_newly_unbalanced_not_audited_ids_sha256")
        and {"R613", "R2076", "R1708", "R2004"} <= set(atomic_unaudited),
        "R1521 unaudited frontier ledger drifted",
    )
    local_expected = frontier["local_counterfactual_full_model"].get("expected")
    atomic = frontier["atomic_identity_transaction"]
    atomic_expected = atomic.get("expected")
    _require(
        isinstance(local_expected, dict)
        and set(local_expected) == _BALANCE_EXPECTATION_KEYS
        and local_expected["incident_reaction_count"] == 70
        and local_expected["newly_unbalanced_count"] == 30
        and local_expected["unassessable_reaction_count"] == 23
        and local_expected["repaired_count"] == 0
        and frontier["local_counterfactual_full_model"].get("zero_new_regressions_required") is True,
        "R1521 local full-model gate drifted",
    )
    _require(
        atomic.get("group_ids")
        == [
            "coenzyme_a", "acetyl_coa", "trans_hexacos_2_enoyl_coa",
            "er_vlcfa_3r_hydroxyhexacosanoyl_coa", "tetracosanoyl_coa",
            "3_oxohexacosanoyl_coa", "nad_plus", "nadh",
        ]
        and atomic.get("expected_group_count") == 8
        and atomic.get("curation_path") == "data/coa_protonation_curation.json"
        and isinstance(atomic_expected, dict)
        and set(atomic_expected) == _BALANCE_EXPECTATION_KEYS
        and atomic_expected["migrated_metabolite_count"] == 36
        and atomic_expected["incident_reaction_count"] == 303
        and atomic_expected["newly_unbalanced_count"] == 146
        and atomic_expected["unassessable_reaction_count"] == 75
        and "R364" in atomic_expected["unassessable_reaction_ids"]
        and atomic_expected["repaired_reaction_ids"] == ["R1889", "R570"]
        and atomic.get("zero_new_regressions_required") is True
        and all(
            _is_digest(value)
            for expectation in (local_expected, atomic_expected)
            for key, value in expectation.items()
            if key.endswith("sha256")
        ),
        "R1521 atomic identity transaction drifted",
    )
    r364 = frontier["r364_symbolic_contract"]
    r364_gene = r364.get("gene_evidence", {})
    _require(
        r364.get("reaction_id") == "R364"
        and r364.get("canonical_h_coefficient") == 1
        and r364.get("model_h_coefficient") == 1
        and r364.get("h_edit_authorized") is False
        and r364.get("reaction_identity_evidence_status") == "supported"
        and r364.get("overall_evidence_status") == "partially_supported"
        and r364.get("balance_assessment_status")
        == "unassessable_symbolic_moiety"
        and r364.get("symbolic_metabolite_identity_status") == "unresolved"
        and r364_gene.get("systematic_id") == "YALI1D26360g"
        and r364_gene.get("established_name") == "dihydrolipoyl dehydrogenase"
        and r364_gene.get("legacy_locus_tag") == "YALI0D20768g"
        and r364_gene.get("protein_accessions") == ["Q6C8C6"]
        and r364_gene.get("model_gpr") == "YALI1D26360g"
        and r364_gene.get("reaction_ids") == ["R364"]
        and r364_gene.get("model_gpr_evidence_status", "").startswith(
            "model/GPR assignment only"
        )
        and r364_gene.get("evidence_status", {}).get("exact_r364_substrate_context")
        == "unverified",
        "R364 symbolic identity/GPR contract drifted",
    )

    contracts = handoff["reaction_contracts"]
    _require(isinstance(contracts, list) and [entry.get("id") for entry in contracts] == ["R1504", "R1521", "R80"], "R1521 adjacent reaction inventory drifted")
    for contract in contracts:
        _require(set(contract) == {"id", "bounds", "reversible", "stoichiometry"}, "invalid R1521 reaction contract")
        _require(
            isinstance(contract["bounds"], list)
            and len(contract["bounds"]) == 2
            and all(type(value) in (int, float) and isfinite(value) for value in contract["bounds"])
            and type(contract["reversible"]) is bool
            and isinstance(contract["stoichiometry"], dict)
            and contract["stoichiometry"]
            and all(isinstance(key, str) and key and type(value) in (int, float) and value for key, value in contract["stoichiometry"].items()),
            f"invalid R1521 reaction contract for {contract.get('id')!r}",
        )
    balances = handoff["partial_known_formula_residual_contract"]
    _require(isinstance(balances, dict) and set(balances) == {"R1504", "R1521", "R80"}, "R1521 balance contract drifted")
    _require(all(isinstance(value, dict) and all(isinstance(key, str) and type(number) in (int, float) for key, number in value.items()) for value in balances.values()), "invalid R1521 partial-residual contract")
    _require(
        isinstance(handoff["blockers"], list)
        and len(handoff["blockers"]) == 5
        and all(isinstance(value, str) and value for value in handoff["blockers"])
        and any("R364" in value for value in handoff["blockers"])
        and any("51" in value for value in handoff["blockers"])
        and any(
            all(reaction_id in value for reaction_id in ("R613", "R2076", "R1708", "R2004"))
            for value in handoff["blockers"]
        ),
        "R1521 blockers must be explicit",
    )


def _validate_runtime_evidence(handoff: dict) -> None:
    """Require every locked evidence file for every handoff use."""

    for name, dependency in handoff["evidence_dependencies"].items():
        path = REPOSITORY / dependency["path"]
        declared = dependency["canonical_contract_sha256"]
        _require(
            _canonical_json_file_digest(path) == declared,
            f"R1521 evidence dependency drifted: {name}",
        )


def _validate_runtime_handoff(handoff: dict) -> None:
    _validate_handoff(handoff)
    _validate_runtime_evidence(handoff)


def load_r1521_current_snapshot_handoff(path: Path | None = None) -> dict:
    """Load the locked, blocked handoff; source validation happens at audit/apply time."""

    try:
        handoff = json.loads((path or HANDOFF_PATH).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise R1521CurrentSnapshotError(f"cannot read R1521 handoff: {exc}") from exc
    _validate_runtime_handoff(handoff)
    return handoff


def _identity_key(name: str | None) -> str:
    return (name or "").strip().lower().rstrip("_")


def validate_r1521_current_snapshot_model(model, handoff: dict | None = None) -> dict:
    """Check exact legacy tuples and adjacent reaction contracts without mutating ``model``."""

    handoff = load_r1521_current_snapshot_handoff() if handoff is None else handoff
    _validate_runtime_handoff(handoff)
    errors: list[str] = []
    contract = handoff["legacy_metabolite_contract"]
    observed = [
        met for met in model.metabolites
        if _identity_key(met.name) in _R1521_IDENTITY_NAMES
    ]
    expected_ids = {copy["id"] for copy in contract["copies"]}
    observed_ids = {met.id for met in observed}
    if observed_ids != expected_ids:
        errors.append(f"m1546[C_pe]: copy completeness drifted; observed={sorted(observed_ids)!r}, expected={sorted(expected_ids)!r}")
    if len(observed) != contract["expected_copy_count"]:
        errors.append(f"m1546[C_pe]: observed {len(observed)} copies, expected {contract['expected_copy_count']}")
    for copy in contract["copies"]:
        try:
            met = model.metabolites.get_by_id(copy["id"])
        except KeyError:
            errors.append(f"{copy['id']}: missing expected legacy copy")
            continue
        actual_tuple = {"formula": met.formula, "charge": met.charge}
        if actual_tuple != copy["legacy_tuple"]:
            errors.append(f"{copy['id']}: third tuple {actual_tuple!r}; expected exact legacy tuple {copy['legacy_tuple']!r}")
        if met.name != copy["name"]:
            errors.append(f"{copy['id']}: legacy identity drifted to {met.name!r}")
        present_keys = sorted(_IDENTITY_KEYS.intersection(met.annotation or {}))
        if present_keys:
            errors.append(f"{copy['id']}: legacy identity unexpectedly has chemical annotations {present_keys!r}")

    for target in handoff["local_counterfactual"]["target_tuples"]:
        try:
            metabolite = model.metabolites.get_by_id(target["id"])
        except KeyError:
            errors.append(f"{target['id']}: missing local counterfactual member")
            continue
        actual_tuple = {"formula": metabolite.formula, "charge": metabolite.charge}
        if actual_tuple != target["legacy_tuple"]:
            errors.append(f"{target['id']}: third tuple {actual_tuple!r}; expected exact local legacy tuple {target['legacy_tuple']!r}")

    for label in ("nad_plus", "nadh"):
        nad_contract = handoff["global_nad_migration_gate"][label]
        observed = {met.id for met in model.metabolites if _identity_key(met.name) == nad_contract["identity_name"]}
        expected = set(nad_contract["ids"])
        if observed != expected:
            errors.append(f"{label}: copy completeness drifted; observed={sorted(observed)!r}, expected={sorted(expected)!r}")
        for metabolite_id in expected:
            try:
                metabolite = model.metabolites.get_by_id(metabolite_id)
            except KeyError:
                errors.append(f"{label}: missing expected copy {metabolite_id}")
                continue
            actual_tuple = {"formula": metabolite.formula, "charge": metabolite.charge}
            if actual_tuple != nad_contract["legacy_tuple"]:
                errors.append(f"{label}: {metabolite_id} has third tuple {actual_tuple!r}; expected {nad_contract['legacy_tuple']!r}")

    for reaction_contract in handoff["reaction_contracts"]:
        reaction_id = reaction_contract["id"]
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            errors.append(f"{reaction_id}: missing adjacent reaction")
            continue
        expected = {key: float(value) for key, value in reaction_contract["stoichiometry"].items()}
        observed_stoichiometry = {met.id: float(value) for met, value in reaction.metabolites.items()}
        if observed_stoichiometry != expected:
            errors.append(f"{reaction_id}: stoichiometry drifted")
        if not (
            isclose(reaction.lower_bound, reaction_contract["bounds"][0], rel_tol=0.0, abs_tol=1e-12)
            and isclose(reaction.upper_bound, reaction_contract["bounds"][1], rel_tol=0.0, abs_tol=1e-12)
            and bool(reaction.reversibility) is reaction_contract["reversible"]
        ):
            errors.append(f"{reaction_id}: bounds or reversibility drifted")
    return {"valid": not errors, "errors": errors}


def _local_counterfactual_on_copy(model, handoff: dict):
    """Apply only the locked counterfactual to a copy for balance/FBA evidence."""

    candidate = model.copy()
    for target in handoff["local_counterfactual"]["target_tuples"]:
        metabolite = candidate.metabolites.get_by_id(target["id"])
        actual = {"formula": metabolite.formula, "charge": metabolite.charge}
        _require(actual == target["legacy_tuple"], f"refusing local third tuple for {metabolite.id}: {actual!r}")
        metabolite.formula, metabolite.charge = target["formula"], target["charge"]
        if "name" in target:
            metabolite.name = target["name"]
        if "annotation_target" in target:
            annotation = dict(metabolite.annotation or {})
            annotation.update(target["annotation_target"])
            metabolite.annotation = annotation
    return candidate


def _ids_digest(ids: list[str]) -> str:
    return hashlib.sha256(json.dumps(sorted(ids), separators=(",", ":")).encode()).hexdigest()


def _complete_mass_balance(reaction) -> dict | None:
    """Return exact element/charge balance only when every tuple is parseable."""

    for metabolite in reaction.metabolites:
        if (
            metabolite.formula is None
            or metabolite.charge is None
            or _COMPLETE_FORMULA.fullmatch(metabolite.formula) is None
        ):
            return None
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            elements = metabolite.elements
        if (
            caught
            or not elements
            or not set(elements) <= set(elements_and_molecular_weights)
        ):
            return None
    balance = _safe_mass_balance(reaction)
    if balance is None or any(not isfinite(float(value)) for value in balance.values()):
        return None
    return {key: float(value) for key, value in sorted(balance.items())}


def _balance_ledger(model, candidate, migrated_ids: list[str], expected: dict) -> dict:
    """Classify every incident reaction, keeping unparseable balances fail-closed."""

    incident = sorted(
        {
            reaction.id
            for metabolite_id in migrated_ids
            for reaction in model.metabolites.get_by_id(metabolite_id).reactions
        }
    )
    categories = {
        name: []
        for name in (
            "changed", "newly_unbalanced", "repaired", "still_unbalanced_changed",
            "balanced_unchanged", "unbalanced_unchanged", "unassessable",
        )
    }
    rows = []
    for reaction_id in incident:
        before = _complete_mass_balance(model.reactions.get_by_id(reaction_id))
        after = _complete_mass_balance(candidate.reactions.get_by_id(reaction_id))
        rows.append({"reaction_id": reaction_id, "before": before, "after": after})
        if before is None or after is None:
            categories["unassessable"].append(reaction_id)
        elif before != after:
            categories["changed"].append(reaction_id)
            if before == {}:
                categories["newly_unbalanced"].append(reaction_id)
            elif after == {}:
                categories["repaired"].append(reaction_id)
            else:
                categories["still_unbalanced_changed"].append(reaction_id)
        elif before == {}:
            categories["balanced_unchanged"].append(reaction_id)
        else:
            categories["unbalanced_unchanged"].append(reaction_id)
    actual = {
        "migrated_metabolite_count": len(migrated_ids),
        "incident_reaction_count": len(incident),
        "evaluable_reaction_count": len(incident) - len(categories["unassessable"]),
        "unassessable_reaction_count": len(categories["unassessable"]),
        "changed_residual_count": len(categories["changed"]),
        "newly_unbalanced_count": len(categories["newly_unbalanced"]),
        "repaired_count": len(categories["repaired"]),
        "still_unbalanced_changed_count": len(categories["still_unbalanced_changed"]),
        "balanced_unchanged_count": len(categories["balanced_unchanged"]),
        "unbalanced_unchanged_count": len(categories["unbalanced_unchanged"]),
        "migrated_metabolite_ids_sha256": _ids_digest(migrated_ids),
        "incident_reaction_ids_sha256": _ids_digest(incident),
        "changed_reaction_ids_sha256": _ids_digest(categories["changed"]),
        "newly_unbalanced_ids_sha256": _ids_digest(categories["newly_unbalanced"]),
        "repaired_ids_sha256": _ids_digest(categories["repaired"]),
        "still_unbalanced_changed_ids_sha256": _ids_digest(
            categories["still_unbalanced_changed"]
        ),
        "unassessable_ids_sha256": _ids_digest(categories["unassessable"]),
        "ledger_sha256": hashlib.sha256(
            json.dumps(
                rows,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
                allow_nan=False,
            ).encode()
        ).hexdigest(),
        "unassessable_reaction_ids": categories["unassessable"],
        "repaired_reaction_ids": categories["repaired"],
    }
    mismatched_fields = [
        key for key, value in expected.items() if actual[key] != value
    ]
    matches_expected = not mismatched_fields
    return {
        "actual": actual,
        "matches_expected": matches_expected,
        "mismatched_fields": mismatched_fields,
        "zero_new_regressions": matches_expected
        and actual["newly_unbalanced_count"] == 0,
        "frontier_closed": matches_expected
        and actual["newly_unbalanced_count"] == 0
        and actual["still_unbalanced_changed_count"] == 0
        and actual["unbalanced_unchanged_count"] == 0
        and actual["unassessable_reaction_count"] == 0,
    }


def _global_nad_migration_audit(model, handoff: dict) -> dict:
    """Quantify why all-copy NAD migration is blocked, entirely on a copy."""

    gate = handoff["global_nad_migration_gate"]
    candidate = model.copy()
    migrated_ids: list[str] = []
    for group in ("nad_plus", "nadh"):
        target = gate[group]["target_tuple"]
        for metabolite_id in gate[group]["ids"]:
            metabolite = candidate.metabolites.get_by_id(metabolite_id)
            actual = {"formula": metabolite.formula, "charge": metabolite.charge}
            _require(actual == gate[group]["legacy_tuple"], f"refusing global {group} third tuple for {metabolite.id}: {actual!r}")
            metabolite.formula, metabolite.charge = target["formula"], target["charge"]
            migrated_ids.append(metabolite_id)
    return _balance_ledger(model, candidate, migrated_ids, gate["expected"])


def _atomic_identity_frontier_audit(
    model, handoff: dict, source_model_path: Path
) -> dict:
    """Apply the eight existing identity groups on a copy and audit all incidents."""

    gate = handoff["frontier_closure_gate"]["atomic_identity_transaction"]
    curation = load_coa_protonation_curation(
        COA_CURATION_PATH, source_model_path=source_model_path
    )
    groups = {group["id"]: group for group in curation["groups"]}
    _require(
        set(gate["group_ids"]) <= set(groups),
        "R1521 atomic identity group is absent from CoA curation",
    )
    candidate = model.copy()
    migrated_ids: list[str] = []
    for group_id in gate["group_ids"]:
        group = groups[group_id]
        legacy_by_id = {
            metabolite_id: {"formula": row["formula"], "charge": row["charge"]}
            for row in group["legacy_tuples"]
            for metabolite_id in row["ids"]
        }
        _require(
            set(legacy_by_id) == set(group["expected_ids"]),
            f"R1521 atomic identity inventory drifted for {group_id}",
        )
        for metabolite_id in group["expected_ids"]:
            source_metabolite = model.metabolites.get_by_id(metabolite_id)
            actual = {
                "formula": source_metabolite.formula,
                "charge": source_metabolite.charge,
            }
            _require(
                actual == legacy_by_id[metabolite_id],
                f"refusing atomic third tuple for {metabolite_id}: {actual!r}",
            )
            metabolite = candidate.metabolites.get_by_id(metabolite_id)
            target = group["target_tuple"]
            metabolite.formula, metabolite.charge = target["formula"], target["charge"]
            migrated_ids.append(metabolite_id)
    _require(
        len(migrated_ids) == gate["expected"]["migrated_metabolite_count"],
        "R1521 atomic identity copy count drifted",
    )
    return _balance_ledger(model, candidate, migrated_ids, gate["expected"])


def _source_snapshot_report(model_path: Path, model, handoff: dict) -> dict:
    actual_sha = _sha256_file(model_path)
    actual_fingerprint = _model_snapshot_fingerprint(model)
    return {
        "declared_model_locator": handoff["source_model"],
        "physical_input_locator_recorded_in_execution_metadata": True,
        "declared_sha256": handoff["source_sha256"],
        "actual_sha256": actual_sha,
        "sha256_verified": actual_sha == handoff["source_sha256"],
        "declared_fingerprint": handoff["source_model_fingerprint"],
        "actual_fingerprint": actual_fingerprint,
        "fingerprint_verified": actual_fingerprint == handoff["source_model_fingerprint"],
    }


def _source_path(model_path: Path | None) -> Path:
    _require(model_path is not None, "explicit R1521 source model path is required")
    return Path(model_path).resolve()


def _validate_output_path(
    output_path: Path, source_path: Path, handoff: dict | None = None
) -> None:
    handoff = load_r1521_current_snapshot_handoff() if handoff is None else handoff
    _validate_runtime_handoff(handoff)
    protected = {
        source_path.resolve(),
        (REPOSITORY / handoff["source_model"]).resolve(),
        HANDOFF_PATH.resolve(),
        COA_CURATION_PATH.resolve(),
        *(
            (REPOSITORY / dependency["path"]).resolve()
            for dependency in handoff["evidence_dependencies"].values()
        ),
    }
    _require(
        output_path.resolve() not in protected,
        "R1521 audit output must not overwrite the model or an authoritative input",
    )


def audit_r1521_current_snapshot(model_path: Path | None = None, handoff: dict | None = None) -> dict:
    """Produce a current-snapshot, blocked report; no candidate model is written."""

    handoff = load_r1521_current_snapshot_handoff() if handoff is None else handoff
    _validate_runtime_handoff(handoff)
    source = _source_path(model_path)
    _require(source.is_file(), f"R1521 source model is unavailable: {source}")
    model = read_sbml_model(str(source))
    validation = validate_r1521_current_snapshot_model(model, handoff)
    snapshot = _source_snapshot_report(source, model, handoff)
    _require(validation["valid"], "; ".join(validation["errors"]))
    _require(snapshot["sha256_verified"] and snapshot["fingerprint_verified"], "R1521 source snapshot drifted; rebase handoff before any counterfactual")
    balances = {}
    for contract in handoff["reaction_contracts"]:
        reaction = model.reactions.get_by_id(contract["id"])
        actual = _safe_mass_balance(reaction)
        expected = handoff["partial_known_formula_residual_contract"][reaction.id]
        balances[reaction.id] = {
            "partial_residual": actual,
            "expected_partial_residual": expected,
            "matches_contract": actual == expected,
            "exact_balance_assessable": _complete_mass_balance(reaction) is not None,
        }
    objective_before = model.slim_optimize(error_value=None)
    candidate = _local_counterfactual_on_copy(model, handoff)
    candidate_balances = {
        contract["id"]: _complete_mass_balance(
            candidate.reactions.get_by_id(contract["id"])
        )
        for contract in handoff["reaction_contracts"]
    }
    local_gate = _balance_ledger(
        model,
        candidate,
        [entry["id"] for entry in handoff["local_counterfactual"]["target_tuples"]],
        handoff["frontier_closure_gate"]["local_counterfactual_full_model"]["expected"],
    )
    objective_after = candidate.slim_optimize(error_value=None)
    objective_unchanged = (
        objective_before is not None
        and objective_after is not None
        and isclose(float(objective_before), float(objective_after), rel_tol=0.0, abs_tol=1e-9)
    )
    global_nad_gate = _global_nad_migration_audit(model, handoff)
    atomic_gate = _atomic_identity_frontier_audit(model, handoff, source)
    for label, gate in (
        ("local counterfactual", local_gate),
        ("global NAD migration", global_nad_gate),
        ("atomic identity transaction", atomic_gate),
    ):
        _require(
            gate["matches_expected"],
            f"R1521 {label} drifted from the locked current snapshot: "
            f"{gate['mismatched_fields']!r}",
        )
    coverage = handoff["frontier_closure_gate"]["source_audit_coverage"]
    r364 = handoff["frontier_closure_gate"]["r364_symbolic_contract"]
    return {
        "schema_version": 2,
        "artifact_type": "r1521_current_snapshot_handoff_audit",
        "evidence_role": "local_cross_check_not_hpcc_acceptance",
        "activation": handoff["activation"],
        "handoff_contract_sha256": handoff["target_contract_sha256"],
        "evidence_dependencies": handoff["evidence_dependencies"],
        "source_snapshot": snapshot,
        "validation": validation,
        "gpr_evidence": handoff["local_counterfactual"]["gpr_evidence"],
        "annotation_conflicts": handoff["local_counterfactual"]["annotation_conflicts"],
        "source_partial_known_formula_residual": balances,
        "candidate_adjacent_exact_balance": candidate_balances,
        "candidate_balances_named_adjacent_reactions_exactly": all(
            value == {} for value in candidate_balances.values()
        ),
        "local_counterfactual_full_model_gate": local_gate,
        "global_nad_migration_gate": global_nad_gate,
        "atomic_identity_frontier_gate": atomic_gate,
        "frontier_source_audit_coverage": coverage,
        "r364_symbolic_contract": r364,
        "runtime_inputs": {
            "atomic_identity_curation": {
                "declared_locator": handoff["frontier_closure_gate"]
                ["atomic_identity_transaction"]["curation_path"],
                "sha256": _sha256_file(COA_CURATION_PATH),
            }
        },
        "frontier_closed": local_gate["frontier_closed"]
        and global_nad_gate["frontier_closed"]
        and atomic_gate["frontier_closed"]
        and coverage["all_atomic_frontier_identities_closed"]
        and r364["symbolic_metabolite_identity_status"] == "resolved",
        "objective_before": None if objective_before is None else round(float(objective_before), 12),
        "objective_after": None if objective_after is None else round(float(objective_after), 12),
        "objective_unchanged": objective_unchanged,
        "blockers": handoff["blockers"],
        "production_gate_passed": False,
        "ready_for_activation": False,
    }


def apply_r1521_current_snapshot_handoff(model, handoff: dict | None = None) -> None:
    """Fail before writes: a blocked handoff has no proposed target transaction."""

    handoff = load_r1521_current_snapshot_handoff() if handoff is None else handoff
    _validate_runtime_handoff(handoff)
    validation = validate_r1521_current_snapshot_model(model, handoff)
    _require(validation["valid"], "; ".join(validation["errors"]))
    _require(_model_snapshot_fingerprint(model) == handoff["source_model_fingerprint"], "R1521 source fingerprint drifted before blocked apply")
    raise CoAProtonationActivationBlocked(
        "R1521 current-snapshot handoff is blocked; no stereochemistry, tuple, or proton repair is authorized"
    )


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write(json.dumps(value, sort_keys=True, indent=2) + "\n")
        temporary = Path(handle.name)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    handoff = load_r1521_current_snapshot_handoff()
    source = _source_path(args.model)
    _validate_output_path(args.output, source, handoff)
    _write_json(args.output, audit_r1521_current_snapshot(source, handoff))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
