#!/usr/bin/env python3
"""Audit the current-main R1521 handoff without proposing a model patch."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from math import isclose, isfinite
from pathlib import Path

from cobra.io import read_sbml_model

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.gem_annotate.patches import (
    CoAProtonationActivationBlocked,
    _model_snapshot_fingerprint,
    _safe_mass_balance,
)


REPOSITORY = Path(__file__).resolve().parents[1]
HANDOFF_PATH = REPOSITORY / "data" / "r1521_current_snapshot_handoff.json"
ER_EVIDENCE_PATH = REPOSITORY / "data" / "er_vlcfa_3r_stereochemistry.json"
_TRUSTED_CONTRACT_SHA256 = "2f80ce48c930c948c0c0a0d38151b8a40b49f4af7babae586b4b16e04431b5b0"
_IDENTITY_KEYS = frozenset(
    {
        "chebi", "hmdb", "inchi", "inchikey", "lipidmaps", "lipidmapsm",
        "metacyc.compound", "metacycm", "metanetx.chemical", "seed.compound", "seedm",
    }
)
_R1521_IDENTITY_NAMES = frozenset(
    {"(r)-3-hydroxyhexacosanoyl-coa", "(3r)-3-hydroxyhexacosanoyl-coa"}
)


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
            "global_nad_migration_gate", "reaction_contracts", "absolute_balance_contract",
            "blockers", "target_contract_sha256",
        },
        "unknown or missing R1521 handoff key",
    )
    _require(type(handoff["schema_version"]) is int and handoff["schema_version"] == 1, "unsupported R1521 handoff schema")
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
            }
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
    _require(isinstance(candidate, dict) and set(candidate) == {"evidence", "supersedes", "rhea", "target_tuples"}, "invalid R1521 candidate transaction")
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
    global_gate = handoff["global_nad_migration_gate"]
    _require(isinstance(global_gate, dict) and set(global_gate) == {"baseline", "nad_plus", "nadh", "expected", "production_apply_forbidden"} and global_gate["baseline"] == "source_snapshot_only; no CoA/local transaction pre-applied" and global_gate["production_apply_forbidden"] is True, "invalid global NAD migration gate")
    for label, identity_name, ids, legacy, target in (
        ("nad_plus", "nad_c21h27n7o14p2", ["m27[C_mi]", "m116[C_pe]", "m122[C_cy]", "m817[C_nu]", "m960[C_er]", "m1458[C_em]"], {"formula": "C21H27N7O14P2", "charge": 0}, {"formula": "C21H26N7O14P2", "charge": -1}),
        ("nadh", "nadh_c21h29n7o14p2", ["m30[C_mi]", "m119[C_pe]", "m123[C_cy]", "m957[C_er]", "m1457[C_em]"], {"formula": "C21H29N7O14P2", "charge": 0}, {"formula": "C21H27N7O14P2", "charge": -2}),
    ):
        entry = global_gate[label]
        _require(isinstance(entry, dict) and entry == {"identity_name": entry.get("identity_name"), "ids": entry.get("ids"), "legacy_tuple": entry.get("legacy_tuple"), "target_tuple": entry.get("target_tuple")} and entry["identity_name"] == identity_name and entry["ids"] == ids and entry["legacy_tuple"] == legacy and entry["target_tuple"] == target, f"R1521 global {label} migration drifted")
    expected = global_gate["expected"]
    _require(isinstance(expected, dict) and set(expected) == {"incident_reaction_count", "changed_residual_count", "newly_unbalanced_count", "repaired_count", "incident_reaction_ids_sha256", "changed_reaction_ids_sha256", "newly_unbalanced_ids_sha256", "representative_newly_unbalanced"}, "invalid global NAD migration expectation")
    _require((expected["incident_reaction_count"], expected["changed_residual_count"], expected["newly_unbalanced_count"], expected["repaired_count"]) == (121, 114, 73, 2), "global NAD migration counts drifted")
    _require(all(_is_digest(expected[key]) for key in ("incident_reaction_ids_sha256", "changed_reaction_ids_sha256", "newly_unbalanced_ids_sha256")), "global NAD migration digest drifted")
    _require(expected["representative_newly_unbalanced"] == ["R1517", "R1518", "R1522", "R1524", "R535", "R564"], "global NAD representative reactions drifted")

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
    balances = handoff["absolute_balance_contract"]
    _require(isinstance(balances, dict) and set(balances) == {"R1504", "R1521", "R80"}, "R1521 balance contract drifted")
    _require(all(isinstance(value, dict) and all(isinstance(key, str) and type(number) in (int, float) for key, number in value.items()) for value in balances.values()), "invalid R1521 absolute-balance contract")
    _require(isinstance(handoff["blockers"], list) and len(handoff["blockers"]) == 3 and all(isinstance(value, str) and value for value in handoff["blockers"]) and "stale metadata" in handoff["blockers"][2], "R1521 blockers must be explicit")


def _validate_runtime_evidence(handoff: dict) -> None:
    """Require the current ER evidence file for every handoff use."""

    declared_er_digest = handoff["evidence_dependencies"]["er_vlcfa_3r_stereochemistry"]["canonical_contract_sha256"]
    _require(
        _canonical_json_file_digest(ER_EVIDENCE_PATH) == declared_er_digest,
        "R1521 ER stereochemistry evidence dependency drifted",
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


def _global_nad_migration_audit(model, handoff: dict) -> dict:
    """Quantify why all-copy NAD migration is blocked, entirely on a copy."""

    gate = handoff["global_nad_migration_gate"]
    before = {reaction.id: _safe_mass_balance(reaction) for reaction in model.reactions}
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
    after = {reaction.id: _safe_mass_balance(reaction) for reaction in candidate.reactions}
    incident = sorted({reaction.id for metabolite_id in migrated_ids for reaction in model.metabolites.get_by_id(metabolite_id).reactions})
    changed = sorted(reaction_id for reaction_id in before if before[reaction_id] != after[reaction_id])
    newly_unbalanced = sorted(reaction_id for reaction_id in before if before[reaction_id] == {} and after[reaction_id] not in ({}, None))
    repaired = sorted(reaction_id for reaction_id in before if before[reaction_id] not in ({}, None) and after[reaction_id] == {})
    actual = {
        "incident_reaction_count": len(incident),
        "changed_residual_count": len(changed),
        "newly_unbalanced_count": len(newly_unbalanced),
        "repaired_count": len(repaired),
        "incident_reaction_ids_sha256": _ids_digest(incident),
        "changed_reaction_ids_sha256": _ids_digest(changed),
        "newly_unbalanced_ids_sha256": _ids_digest(newly_unbalanced),
        "representative_newly_unbalanced_present": all(reaction_id in newly_unbalanced for reaction_id in gate["expected"]["representative_newly_unbalanced"]),
    }
    expected = gate["expected"]
    checks = {key: actual[key] == expected[key] for key in expected if key != "representative_newly_unbalanced"}
    checks["representative_newly_unbalanced"] = actual["representative_newly_unbalanced_present"]
    return {"actual": actual, "matches_expected": all(checks.values()), "checks": checks}


def _source_snapshot_report(model_path: Path, model, handoff: dict) -> dict:
    actual_sha = _sha256_file(model_path)
    actual_fingerprint = _model_snapshot_fingerprint(model)
    return {
        "model_path": handoff["source_model"],
        "declared_sha256": handoff["source_sha256"],
        "actual_sha256": actual_sha,
        "sha256_verified": actual_sha == handoff["source_sha256"],
        "declared_fingerprint": handoff["source_model_fingerprint"],
        "actual_fingerprint": actual_fingerprint,
        "fingerprint_verified": actual_fingerprint == handoff["source_model_fingerprint"],
    }


def _source_path(model_path: Path | None, handoff: dict) -> Path:
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
        ER_EVIDENCE_PATH.resolve(),
    }
    _require(
        output_path.resolve() not in protected,
        "R1521 audit output must not overwrite the model or an authoritative input",
    )


def audit_r1521_current_snapshot(model_path: Path | None = None, handoff: dict | None = None) -> dict:
    """Produce a current-snapshot, blocked report; no candidate model is written."""

    handoff = load_r1521_current_snapshot_handoff() if handoff is None else handoff
    _validate_runtime_handoff(handoff)
    source = _source_path(model_path, handoff)
    _require(source.is_file(), f"R1521 source model is unavailable: {source}")
    model = read_sbml_model(str(source))
    validation = validate_r1521_current_snapshot_model(model, handoff)
    snapshot = _source_snapshot_report(source, model, handoff)
    _require(validation["valid"], "; ".join(validation["errors"]))
    _require(snapshot["sha256_verified"] and snapshot["fingerprint_verified"], "R1521 source snapshot drifted; rebase handoff before any counterfactual")
    balances = {}
    for contract in handoff["reaction_contracts"]:
        reaction = model.reactions.get_by_id(contract["id"])
        actual = reaction.check_mass_balance()
        expected = handoff["absolute_balance_contract"][reaction.id]
        balances[reaction.id] = {
            "actual": actual,
            "expected": expected,
            "matches_exactly": actual == expected,
            "is_balanced": actual == {},
        }
    objective_before = model.slim_optimize(error_value=None)
    candidate = _local_counterfactual_on_copy(model, handoff)
    candidate_balances = {
        contract["id"]: candidate.reactions.get_by_id(contract["id"]).check_mass_balance()
        for contract in handoff["reaction_contracts"]
    }
    objective_after = candidate.slim_optimize(error_value=None)
    objective_unchanged = (
        objective_before is not None
        and objective_after is not None
        and isclose(float(objective_before), float(objective_after), rel_tol=0.0, abs_tol=1e-9)
    )
    global_nad_gate = _global_nad_migration_audit(model, handoff)
    _require(
        global_nad_gate["matches_expected"],
        "R1521 global NAD migration gate drifted from the locked current snapshot",
    )
    return {
        "schema_version": 1,
        "artifact_type": "r1521_current_snapshot_handoff_audit",
        "activation": handoff["activation"],
        "handoff_contract_sha256": handoff["target_contract_sha256"],
        "evidence_dependencies": handoff["evidence_dependencies"],
        "source_snapshot": snapshot,
        "validation": validation,
        "adjacent_absolute_balance": balances,
        "candidate_adjacent_absolute_balance": candidate_balances,
        "candidate_closes_adjacent_reactions": all(value == {} for value in candidate_balances.values()),
        "global_nad_migration_gate": global_nad_gate,
        "objective_before": None if objective_before is None else round(float(objective_before), 12),
        "objective_after": None if objective_after is None else round(float(objective_after), 12),
        "objective_unchanged": objective_unchanged,
        "blockers": handoff["blockers"],
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
    source = _source_path(args.model, handoff)
    _validate_output_path(args.output, source, handoff)
    _write_json(args.output, audit_r1521_current_snapshot(source, handoff))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
