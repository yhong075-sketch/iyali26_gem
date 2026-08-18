#!/usr/bin/env python3
"""Validate the controlled fatty acyl-CoA handoff and emit a read-only report."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

from cobra.io import read_sbml_model

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.gem_annotate.patches import (
    _model_snapshot_fingerprint,
    audit_coa_protonation_curation,
    load_coa_protonation_curation,
)
from scripts.r1521_current_snapshot_handoff import (
    ER_EVIDENCE_PATH,
    HANDOFF_PATH as R1521_HANDOFF_PATH,
    load_r1521_current_snapshot_handoff,
)


REPOSITORY = Path(__file__).resolve().parents[1]
HANDOFF_PATH = REPOSITORY / "data" / "fatty_acyl_coa_handoff.json"
COA_CURATION_PATH = REPOSITORY / "data" / "coa_protonation_curation.json"
LEDGER_PATH = REPOSITORY / "data" / "lipid_moiety_ledger_spec.json"
_GROUP_IDS = (
    "lauroyl_coa",
    "myristoyl_coa",
    "palmitoyl_coa",
    "palmitoleoyl_coa",
    "stearoyl_coa",
    "oleoyl_coa",
    "linoleoyl_coa",
)


class HandoffError(RuntimeError):
    """The handoff is incomplete, ambiguous, or inconsistent with its inputs."""


def _read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise HandoffError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise HandoffError(f"JSON root must be an object: {path}")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_handoff(path: Path = HANDOFF_PATH) -> dict:
    """Load and fail closed on the handoff's machine-readable invariants."""

    handoff = _read_json(path)
    if handoff.get("schema_version") != 1:
        raise HandoffError("unsupported handoff schema")
    activation = handoff.get("activation")
    if not isinstance(activation, dict) or activation != {
        "state": "blocked",
        "ready_for_activation": False,
        "production_apply_forbidden": True,
        "approval_rule": activation.get("approval_rule"),
    } or not activation["approval_rule"]:
        raise HandoffError("handoff activation state must be explicitly blocked")

    groups = handoff.get("concrete_acyl_coa_groups")
    if not isinstance(groups, list) or tuple(group.get("id") for group in groups) != _GROUP_IDS:
        raise HandoffError("handoff must contain the seven controlled acyl-CoA groups")
    copies = [copy for group in groups for copy in group.get("copies", [])]
    if len(copies) != 36 or len({copy.get("id") for copy in copies}) != 36:
        raise HandoffError("handoff must contain exactly 36 unique concrete copies")
    if not all(
        isinstance(copy.get("id"), str)
        and copy.get("compartment") in handoff["compartments"]
        for copy in copies
    ):
        raise HandoffError("concrete copy IDs or compartments are invalid")

    pools = handoff.get("generic_acyl_coa_pools")
    expected_pool_ids = ("xPOOL_AC_EM", "xPOOL_AC_LP", "xPOOL_AC_MM")
    if not isinstance(pools, list) or tuple(pool.get("reaction_id") for pool in pools) != expected_pool_ids:
        raise HandoffError("handoff must contain exactly the three generic acyl-CoA pools")
    for pool in pools:
        if (
            pool.get("formula") is not None
            or pool.get("charge") is not None
            or pool.get("live_model_charge_policy") != "zero_placeholder_only"
            or pool.get("status") != "generic_nonchemical_interface"
            or tuple(pool.get("member_group_ids", [])) != _GROUP_IDS
        ):
            raise HandoffError("generic acyl-CoA pools must not carry chemical tuples")
    if handoff.get("acceptance_checks", {}).get("ready_for_activation_expected") is not False:
        raise HandoffError("handoff must not permit activation")
    return handoff


def validate_handoff_inputs(
    handoff: dict, *, source_model_path: str | Path | None = None
) -> None:
    """Cross-check handoff entries against the two canonical pipeline JSON inputs."""

    curation = load_coa_protonation_curation(
        source_model_path=source_model_path
    )
    ledger = _read_json(LEDGER_PATH)
    r1521_handoff = load_r1521_current_snapshot_handoff()
    declared = handoff["authoritative_inputs"]["coa_protonation_curation"]
    if declared["source_model_sha256"] != curation["source_sha256"]:
        raise HandoffError("handoff source SHA differs from CoA curation")
    if declared["source_model_fingerprint"] != curation["source_model_fingerprint"]:
        raise HandoffError("handoff source fingerprint differs from CoA curation")
    current_review = handoff["authoritative_inputs"].get("r1521_current_snapshot_handoff")
    if not isinstance(current_review, dict) or current_review != {
        "path": "data/r1521_current_snapshot_handoff.json",
        "source_model": r1521_handoff["source_model"],
        "source_model_sha256": r1521_handoff["source_sha256"],
        "source_model_fingerprint": r1521_handoff["source_model_fingerprint"],
        "role": current_review.get("role") if isinstance(current_review, dict) else None,
        "activation_state": "blocked",
        "ready_for_activation": False,
    } or not current_review["role"]:
        raise HandoffError("R1521 current-snapshot handoff declaration drifted")
    curation_groups = {group["id"]: group for group in curation["groups"]}
    mappings = ledger["source_chemical_convention"]["acyl_coas"]
    for group in handoff["concrete_acyl_coa_groups"]:
        curation_group = curation_groups.get(group["id"])
        chain_id = group["id"].removesuffix("_coa")
        if curation_group is None or group["target_tuple"] != curation_group["target_tuple"]:
            raise HandoffError(f"target tuple drift for {group['id']}")
        legacy = curation_group["legacy_tuples"]
        if len(legacy) != 1 or group["legacy_tuple"] != {
            "formula": legacy[0]["formula"], "charge": legacy[0]["charge"]
        }:
            raise HandoffError(f"legacy tuple drift for {group['id']}")
        if [copy["id"] for copy in group["copies"]] != curation_group["expected_ids"]:
            raise HandoffError(f"copy inventory drift for {group['id']}")
        if group["canonical_mapping"] != mappings[chain_id]["canonical_mapping"]:
            raise HandoffError(f"canonical mapping drift for {group['id']}")


def audit_handoff(model_path: Path) -> dict:
    """Produce a current-snapshot, read-only acceptance report."""

    model_path = Path(model_path).resolve()
    handoff = load_handoff()
    validate_handoff_inputs(handoff, source_model_path=model_path)
    model = read_sbml_model(str(model_path))
    coa_audit = audit_coa_protonation_curation(
        model, source_model_path=model_path
    )
    groups = {group["id"]: group for group in handoff["concrete_acyl_coa_groups"]}
    pool_rows = []
    for pool in handoff["generic_acyl_coa_pools"]:
        try:
            reaction = model.reactions.get_by_id(pool["reaction_id"])
            metabolite = model.metabolites.get_by_id(pool["product_id"])
        except KeyError as exc:
            raise HandoffError(f"generic acyl-CoA pool is absent from model: {exc.args[0]}") from exc
        expected_members = {
            copy["id"]
            for group_id in pool["member_group_ids"]
            for copy in groups[group_id]["copies"]
            if copy["compartment"] == pool["compartment"]
        }
        actual_members = {
            candidate.id
            for candidate, coefficient in reaction.metabolites.items()
            if coefficient < 0
        }
        pool_rows.append(
            {
                "reaction_id": pool["reaction_id"],
                "product_id": pool["product_id"],
                "formula_is_absent": metabolite.formula is None,
                "handoff_tuple_is_null": pool["formula"] is None and pool["charge"] is None,
                "live_model_charge": metabolite.charge,
                "charge_is_nonchemical_placeholder": metabolite.charge in {None, 0},
                "chebi_is_absent": "chebi" not in (metabolite.annotation or {}),
                "no_complete_chemical_tuple": (
                    metabolite.formula is None
                    and metabolite.charge in {None, 0}
                    and "chebi" not in (metabolite.annotation or {})
                ),
                "product_is_the_only_positive_member": {
                    candidate.id
                    for candidate, coefficient in reaction.metabolites.items()
                    if coefficient > 0
                } == {pool["product_id"]},
                "expected_member_ids": sorted(expected_members),
                "actual_member_ids": sorted(actual_members),
                "inventory_linked": actual_members == expected_members,
            }
        )
    files = {
        "source_model": model_path,
        "handoff": HANDOFF_PATH,
        "coa_curation": COA_CURATION_PATH,
        "lipid_moiety_ledger": LEDGER_PATH,
        "r1521_handoff": R1521_HANDOFF_PATH,
    }
    return {
        "schema_version": 1,
        "artifact_type": "fatty_acyl_coa_handoff_acceptance_report",
        "model_path": handoff["authoritative_inputs"]["coa_protonation_curation"]["source_model"],
        "model_sha256": _sha256(model_path),
        "model_fingerprint": _model_snapshot_fingerprint(model),
        "declared_source_model_sha256": handoff["authoritative_inputs"]["coa_protonation_curation"]["source_model_sha256"],
        "declared_source_model_fingerprint": handoff["authoritative_inputs"]["coa_protonation_curation"]["source_model_fingerprint"],
        "current_model_matches_declared_source": coa_audit["source_snapshot"]["model_fingerprint_verified"],
        "file_sha256": {name: _sha256(path) for name, path in files.items()},
        "inventory": {"specific_groups": 7, "specific_copies": 36, "generic_pools": 3},
        "generic_pool_policy": pool_rows,
        "coa_curation_audit": coa_audit,
        "current_snapshot_component_reviews": {
            "R1521": handoff["authoritative_inputs"]["r1521_current_snapshot_handoff"]
        },
        "remaining_blockers": handoff["remaining_activation_blockers"],
        "ready_for_activation": False,
    }


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, sort_keys=True, indent=2) + "\n"
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write(payload)
        temporary = Path(handle.name)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _declared_input_path(path: str) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else REPOSITORY / candidate


def _validate_output_path(
    output_path: Path, model_path: Path, handoff: dict | None = None
) -> None:
    handoff = load_handoff() if handoff is None else handoff
    inputs = handoff["authoritative_inputs"]
    protected = {
        model_path.resolve(),
        HANDOFF_PATH.resolve(),
        ER_EVIDENCE_PATH.resolve(),
        *(
            _declared_input_path(entry[key]).resolve()
            for entry in inputs.values()
            for key in ("path", "source_model")
            if isinstance(entry, dict) and isinstance(entry.get(key), str)
        ),
    }
    if output_path.resolve() in protected:
        raise HandoffError(
            "handoff report output must not overwrite the model or an authoritative input"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=REPOSITORY / "model.xml")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    handoff = load_handoff()
    _validate_output_path(args.output, args.model, handoff)
    _write_json(args.output, audit_handoff(args.model))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
