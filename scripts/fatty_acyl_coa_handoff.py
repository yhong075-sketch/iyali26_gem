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
    _BALANCE_EXPECTATION_KEYS,
    _balance_ledger,
    _complete_mass_balance,
    _ids_digest,
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
_D_COMPONENT_GROUP_IDS = {
    "r1708": ("acetate", "succinate", "acetyl_coa", "succinyl_coa"),
    "r2004": ("malonyl_coa", "pyruvate", "acetyl_coa", "oxaloacetate"),
    "union": (
        "acetate", "succinate", "acetyl_coa", "succinyl_coa",
        "malonyl_coa", "pyruvate", "oxaloacetate",
    ),
}
_D_FOCAL_REACTION_IDS = ("R1708", "R2004", "R613", "R2076")
_D_REVIEW_SHA256 = "1e7874f2ace5d2288e3c6f062972038cdab7ca33d8df40b6f3f01f1059a4a08c"


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


def _validate_d_review_contract(handoff: dict) -> None:
    review = handoff.get("d_connected_component_review")
    if not isinstance(review, dict):
        raise HandoffError("D connected-component review is absent")
    if hashlib.sha256(
        json.dumps(review, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest() != _D_REVIEW_SHA256:
        raise HandoffError("D connected-component review contract drifted")
    declared = handoff["authoritative_inputs"]["coa_protonation_curation"]
    if (
        review.get("source_model_sha256") != declared["source_model_sha256"]
        or review.get("source_model_fingerprint") != declared["source_model_fingerprint"]
    ):
        raise HandoffError("D review source snapshot drifted")
    if review.get("status") != "current_snapshot_candidate_audited_frontier_open":
        raise HandoffError("D review status must remain fail closed")
    inventory = review.get("group_inventory")
    if (
        not isinstance(inventory, list)
        or tuple(row.get("id") for row in inventory) != _D_COMPONENT_GROUP_IDS["union"]
        or len({metabolite_id for row in inventory for metabolite_id in row.get("expected_ids", [])}) != 25
    ):
        raise HandoffError("D review must contain seven groups and 25 unique copies")
    components = review.get("components")
    if not isinstance(components, dict) or set(components) != set(_D_COMPONENT_GROUP_IDS):
        raise HandoffError("D review component set drifted")
    for name, group_ids in _D_COMPONENT_GROUP_IDS.items():
        component = components[name]
        expected = component.get("expected")
        if (
            tuple(component.get("group_ids", [])) != group_ids
            or not isinstance(expected, dict)
            or set(expected) != _BALANCE_EXPECTATION_KEYS
            or not isinstance(expected.get("ledger_sha256"), str)
            or len(expected["ledger_sha256"]) != 64
        ):
            raise HandoffError(f"D review {name} ledger contract drifted")
    union = components["union"]
    for field, count_field, digest_field in (
        ("newly_unbalanced_reaction_ids", "newly_unbalanced_count", "newly_unbalanced_ids_sha256"),
        ("still_unbalanced_changed_reaction_ids", "still_unbalanced_changed_count", "still_unbalanced_changed_ids_sha256"),
    ):
        reaction_ids = union.get(field)
        expected = union["expected"]
        if (
            not isinstance(reaction_ids, list)
            or reaction_ids != sorted(set(reaction_ids))
            or len(reaction_ids) != expected[count_field]
            or _ids_digest(reaction_ids) != expected[digest_field]
        ):
            raise HandoffError(f"D union {field} drifted")
    if review.get("focal_balance_contracts") != {
        "d_only_union_candidate": {
            "R1708": {}, "R2004": {},
            "R613": {"H": -5.0, "charge": -5.0}, "R2076": {},
        },
        "complete_coa_candidate": {
            "R1708": {}, "R2004": {},
            "R613": {"H": -1.0, "charge": -1.0},
            "R2076": {"H": 4.0, "charge": 4.0},
        },
    }:
        raise HandoffError("D focal balance contract drifted")
    reaction_reviews = review.get("reaction_reviews")
    if not isinstance(reaction_reviews, dict):
        raise HandoffError("D reaction reviews are absent")
    r1708 = reaction_reviews.get("R1708")
    r2004 = reaction_reviews.get("R2004")
    if not isinstance(r1708, dict) or not isinstance(r2004, dict):
        raise HandoffError("D focal reaction reviews are malformed")
    r89 = r1708.get("r89_adjacent_candidate")
    candidates = r2004.get("candidates")
    if (
        not isinstance(r89, dict)
        or not isinstance(candidates, list)
        or not candidates
        or not all(isinstance(candidate, dict) for candidate in candidates)
        or any(value is not False for value in (
            r1708.get("production_authorized"),
            r1708.get("production_decision_applied"),
            r89.get("production_authorized"),
            r89.get("production_decision_applied"),
            r2004.get("production_authorized"),
            r2004.get("production_decision_applied"),
        ))
        or any(candidate.get("production_authorized") is not False for candidate in candidates)
    ):
        raise HandoffError("D reaction candidates must remain production-blocked")
    if (
        r1708.get("candidate_reaction_name") != "acetate:succinate CoA-transferase"
        or r1708.get("candidate_ec") != "2.8.3.18"
        or r1708.get("rejected_legacy_ec") != "3.1.2.1"
        or r1708.get("preserve_database_reaction_ids")
        != {
            "rhea": {
                "master": "35711", "left_to_right": "35712",
                "right_to_left": "35713", "bidirectional": "35714",
            },
            "kegg": "R10343",
        }
        or r1708.get("ph_7_3_aggregate_each_side")
        != {"formula": "C27H38N7O21P3S", "charge": -6}
        or r1708.get("explicit_proton_coefficient") != 0
        or {
            key: r1708.get("gpr_evidence", {}).get(key)
            for key in (
                "systematic_identifier", "established_name",
                "legacy_systematic_identifier", "uniprot_accession",
                "protein_function", "model_gpr", "evidence_status",
            )
        }
        != {
            "systematic_identifier": "YALI1E36437g",
            "established_name": "ACH1",
            "legacy_systematic_identifier": "YALI0E30965g",
            "uniprot_accession": "Q6C3Z9",
            "protein_function": "acetate:succinate CoA-transferase",
            "model_gpr": "YALI1E36437g",
            "evidence_status": "experimentally_supported_function_direct_protein_localization_unverified",
        }
        or r2004.get("current_malonyl_coa_metabolite_id") != "m200[C_cy]"
        or r2004.get("forbidden_historical_metabolite_id") != "m1855[C_cy]"
        or r2004.get("historical_metabolite_status") != "merged_absent_do_not_restore"
        or r2004.get("preserve_compartment") != "C_cy"
        or r2004.get("preserve_reversible_bounds") != [-1000.0, 1000.0]
        or r2004.get("database_chemistry")
        != {
            "kegg_reaction": "R00353",
            "metanetx_reaction": "MNXR106487",
            "modelseed_reaction": "rxn00258",
            "evidence_status": "database_chemistry_supported_non_organism_specific",
        }
        or r2004.get("ph_7_3_aggregate_each_side")
        != {"formula": "C27H36N7O22P3S", "charge": -6}
        or r2004.get("explicit_proton_coefficient") != 0
        or r2004.get("preferred_candidate_decision") != "remove"
        or {candidate.get("decision"): candidate.get("status") for candidate in candidates}
        != {
            "retain_reversible_without_gpr": "allowed_for_sensitivity_only",
            "remove": "recommended_fail_closed_candidate",
            "replace_with_R00930": "forbidden_wrong_substrates",
        }
    ):
        raise HandoffError("D focal scientific decision contract drifted")
    r1708_provenance = r1708.get("provenance")
    identifier_mapping = (
        r1708_provenance.get("identifier_mapping")
        if isinstance(r1708_provenance, dict) else None
    )
    functional_evidence = (
        r1708_provenance.get("functional_evidence")
        if isinstance(r1708_provenance, dict) else None
    )
    proteomics_evidence = (
        r1708_provenance.get("proteomics_evidence")
        if isinstance(r1708_provenance, dict) else None
    )
    if (
        not isinstance(r1708_provenance, dict)
        or not isinstance(identifier_mapping, dict)
        or not isinstance(functional_evidence, dict)
        or not isinstance(proteomics_evidence, dict)
        or r1708_provenance.get("accessed_on") != "2026-08-18"
        or r1708_provenance.get("reaction_identity") != {
            "ec": "2.8.3.18",
            "rhea_master": "35711",
            "rhea_bidirectional": "35714",
            "kegg": "R10343",
            "evidence_status": "curated_reaction_identity",
        }
        or {
            key: identifier_mapping.get(key)
            for key in ("model_identifier", "s2_identifier", "legacy_identifier", "ncbi_gene_id")
        } != {
            "model_identifier": "YALI1E36437g",
            "s2_identifier": "YALI1_E36437g",
            "legacy_identifier": "YALI0E30965g",
            "ncbi_gene_id": "2911920",
        }
        or functional_evidence.get("evidence_status") != "experimentally_supported_condition_specific"
        or proteomics_evidence.get("evidence_status")
        != "condition_specific_protein_detection_not_function_or_organelle_localization"
    ):
        raise HandoffError("R1708 evidence provenance drifted")
    r2004_gpr = r2004.get("gpr_evidence")
    search_provenance = (
        r2004_gpr.get("search_provenance") if isinstance(r2004_gpr, dict) else None
    )
    zero_hit_searches = (
        search_provenance.get("zero_hit_searches")
        if isinstance(search_provenance, dict) else None
    )
    positive_controls = (
        search_provenance.get("positive_controls")
        if isinstance(search_provenance, dict) else None
    )
    if (
        not isinstance(search_provenance, dict)
        or not isinstance(zero_hit_searches, list)
        or not all(isinstance(row, dict) for row in zero_hit_searches)
        or not isinstance(positive_controls, list)
        or not all(isinstance(row, dict) for row in positive_controls)
        or any(r2004_gpr.get(key) is not None for key in (
            "systematic_identifier", "established_name", "protein_function", "model_gpr",
        ))
        or r2004_gpr.get("source_model_gene_reaction_rule") != ""
        or search_provenance.get("accessed_on") != "2026-08-18"
        or len(zero_hit_searches) != 6
        or any(row.get("hit_count") != 0 for row in zero_hit_searches)
        or len(positive_controls) != 5
        or any(row.get("hit_count", 0) <= 0 for row in positive_controls)
    ):
        raise HandoffError("R2004 unverified-GPR evidence provenance drifted")
    if any(review.get(key) is not value for key, value in {
        "frontier_closed": False,
        "production_gate_passed": False,
        "production_apply_performed": False,
        "ready_for_activation": False,
        "human_gate_required": True,
    }.items()):
        raise HandoffError("D review must remain fail closed")


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
    try:
        _validate_d_review_contract(handoff)
    except HandoffError:
        raise
    except (AttributeError, KeyError, TypeError) as exc:
        raise HandoffError(f"D review structure is malformed: {exc}") from exc
    return handoff


def validate_handoff_inputs(
    handoff: dict, *, source_model_path: str | Path | None = None
) -> dict:
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
    d_review = handoff["d_connected_component_review"]
    if [
        {
            "id": row["id"],
            "expected_ids": row["expected_ids"],
            "target_tuple": row["target_tuple"],
        }
        for row in d_review["group_inventory"]
    ] != [
        {
            "id": group_id,
            "expected_ids": curation_groups[group_id]["expected_ids"],
            "target_tuple": curation_groups[group_id]["target_tuple"],
        }
        for group_id in _D_COMPONENT_GROUP_IDS["union"]
    ]:
        raise HandoffError("D review copy inventory drifted from CoA curation")
    contracts = {row["reaction_id"]: row for row in curation["reaction_contracts"]}
    if not {"R1708", "R2004"} <= set(contracts) or any(
        contracts.get(reaction_id, {}).get("source_stoichiometry")
        != contracts.get(reaction_id, {}).get("target_stoichiometry")
        for reaction_id in ("R1708", "R2004")
    ):
        raise HandoffError("D focal reaction source-to-target contract drifted")
    identities = {
        row["reaction_id"]: row
        for row in curation["reaction_identity_corrections"]
    }
    identity = identities.get("R1708")
    if not isinstance(identity, dict) or identity.get("target_name") != "acetate:succinate CoA-transferase" or identity.get("annotation_target") != {"ec-code": "2.8.3.18"}:
        raise HandoffError("R1708 identity correction drifted")
    if "R2004" in identities:
        raise HandoffError("R2004 identity or GPR correction is not authorized")
    if {row["reaction_id"] for row in curation["reaction_corrections"]} & {"R89", "R1708", "R2004"}:
        raise HandoffError("D review does not authorize focal proton corrections")
    blockers = {
        row["reaction_id"]: row.get("residual_after_curated_tuples")
        for row in curation["activation_blockers"]
    }
    if any(blockers.get(reaction_id) != {} for reaction_id in ("R1708", "R2004")) or not {
        "R1708", "R2004"
    } <= set(curation["target_balance_reaction_ids"]):
        raise HandoffError("D focal reactions must remain explicit activation blockers")
    return curation


def _d_candidate(model, curation_groups: dict, group_ids: tuple[str, ...]):
    candidate = model.copy()
    migrated_ids = []
    for group_id in group_ids:
        group = curation_groups[group_id]
        legacy = {
            metabolite_id: (row["formula"], row["charge"])
            for row in group["legacy_tuples"]
            for metabolite_id in row["ids"]
        }
        for metabolite_id in group["expected_ids"]:
            source = model.metabolites.get_by_id(metabolite_id)
            if (source.formula, source.charge) != legacy[metabolite_id]:
                raise HandoffError(f"refusing D third tuple for {metabolite_id}")
            target = candidate.metabolites.get_by_id(metabolite_id)
            target.formula = group["target_tuple"]["formula"]
            target.charge = group["target_tuple"]["charge"]
            migrated_ids.append(metabolite_id)
    return candidate, migrated_ids


def _d_changed_category_ids(model, candidate, migrated_ids: list[str]) -> dict:
    incident = sorted({
        reaction.id
        for metabolite_id in migrated_ids
        for reaction in model.metabolites.get_by_id(metabolite_id).reactions
    })
    categories = {"newly_unbalanced_reaction_ids": [], "still_unbalanced_changed_reaction_ids": []}
    for reaction_id in incident:
        before = _complete_mass_balance(model.reactions.get_by_id(reaction_id))
        after = _complete_mass_balance(candidate.reactions.get_by_id(reaction_id))
        if before is None or after is None or before == after:
            continue
        if before == {}:
            categories["newly_unbalanced_reaction_ids"].append(reaction_id)
        elif after != {}:
            categories["still_unbalanced_changed_reaction_ids"].append(reaction_id)
    return categories


def _audit_d_components(model, curation: dict, coa_audit: dict, review: dict) -> dict:
    groups = {group["id"]: group for group in curation["groups"]}
    ledgers = {}
    union_candidate = None
    union_migrated_ids = None
    for name, group_ids in _D_COMPONENT_GROUP_IDS.items():
        candidate, migrated_ids = _d_candidate(model, groups, group_ids)
        expected = review["components"][name]["expected"]
        first = _balance_ledger(model, candidate, migrated_ids, expected)
        second_candidate, second_migrated_ids = _d_candidate(model, groups, group_ids)
        second = _balance_ledger(model, second_candidate, second_migrated_ids, expected)
        if migrated_ids != second_migrated_ids or first != second or not first["matches_expected"]:
            raise HandoffError(f"D {name} ledger drifted: {first['mismatched_fields']}")
        ledgers[name] = first
        if name == "union":
            union_candidate, union_migrated_ids = candidate, migrated_ids
    categories = _d_changed_category_ids(model, union_candidate, union_migrated_ids)
    declared_union = review["components"]["union"]
    if any(categories[field] != declared_union[field] for field in categories):
        raise HandoffError("D union frontier reaction IDs drifted")
    d_focal = {
        reaction_id: _complete_mass_balance(
            union_candidate.reactions.get_by_id(reaction_id)
        )
        for reaction_id in _D_FOCAL_REACTION_IDS
    }
    full_coa_unbalanced = coa_audit["target_reactions_unbalanced"]
    expected_focal = review["focal_balance_contracts"]
    complete_expected = expected_focal["complete_coa_candidate"]
    if any(
        (expected == {}) == (reaction_id in full_coa_unbalanced)
        for reaction_id, expected in complete_expected.items()
    ):
        raise HandoffError("D complete-CoA focal evaluation state drifted")
    full_coa_focal = {
        reaction_id: (
            {} if expected == {} else full_coa_unbalanced[reaction_id]
        )
        for reaction_id, expected in complete_expected.items()
    }
    if d_focal != expected_focal["d_only_union_candidate"] or full_coa_focal != expected_focal["complete_coa_candidate"]:
        raise HandoffError("D focal balance residual drifted")
    r1708 = model.reactions.get_by_id("R1708")
    r2004 = model.reactions.get_by_id("R2004")
    if (
        r1708.gene_reaction_rule != "YALI1E36437g"
        or r2004.gene_reaction_rule
        or r2004.metabolites.get(model.metabolites.get_by_id("m200[C_cy]")) != -1.0
        or model.metabolites.has_id("m1855[C_cy]")
    ):
        raise HandoffError("D focal runtime identity drifted")
    return {
        "status": review["status"],
        "deterministic_recomputation": True,
        "ledgers": ledgers,
        "union_frontier_reaction_ids": categories,
        "d_only_union_focal_balances": d_focal,
        "complete_coa_candidate_focal_balances": full_coa_focal,
        "reaction_reviews": review["reaction_reviews"],
        "frontier_closed": False,
        "production_gate_passed": False,
        "production_apply_performed": False,
        "ready_for_activation": False,
        "human_gate_required": True,
    }


def audit_handoff(model_path: Path) -> dict:
    """Produce a current-snapshot, read-only acceptance report."""

    model_path = Path(model_path).resolve()
    handoff = load_handoff()
    curation = validate_handoff_inputs(handoff, source_model_path=model_path)
    model = read_sbml_model(str(model_path))
    source_fingerprint = _model_snapshot_fingerprint(model)
    source_focal_identities = {
        reaction_id: (
            model.reactions.get_by_id(reaction_id).name,
            model.reactions.get_by_id(reaction_id).gene_reaction_rule,
            json.dumps(model.reactions.get_by_id(reaction_id).annotation or {}, sort_keys=True),
        )
        for reaction_id in _D_FOCAL_REACTION_IDS
    }
    coa_audit = audit_coa_protonation_curation(
        model, curation=curation, source_model_path=model_path
    )
    d_review = _audit_d_components(
        model, curation, coa_audit, handoff["d_connected_component_review"]
    )
    current_focal_identities = {
        reaction_id: (
            model.reactions.get_by_id(reaction_id).name,
            model.reactions.get_by_id(reaction_id).gene_reaction_rule,
            json.dumps(model.reactions.get_by_id(reaction_id).annotation or {}, sort_keys=True),
        )
        for reaction_id in _D_FOCAL_REACTION_IDS
    }
    if (
        _model_snapshot_fingerprint(model) != source_fingerprint
        or current_focal_identities != source_focal_identities
    ):
        raise HandoffError("D audit mutated the source model")
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
        "evidence_role": "local_cross_check_not_hpcc_acceptance",
        "declared_model_locator": handoff["authoritative_inputs"]["coa_protonation_curation"]["source_model"],
        "physical_input_locator_recorded_in_execution_metadata": True,
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
            "R1521": handoff["authoritative_inputs"]["r1521_current_snapshot_handoff"],
            "D_R1708_R2004": d_review,
        },
        "remaining_blockers": handoff["remaining_activation_blockers"],
        "activation_state": handoff["activation"]["state"],
        "production_gate_passed": False,
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
    r1521 = load_r1521_current_snapshot_handoff()
    inputs = handoff["authoritative_inputs"]
    protected = {
        model_path.resolve(),
        HANDOFF_PATH.resolve(),
        ER_EVIDENCE_PATH.resolve(),
        *(
            (REPOSITORY / dependency["path"]).resolve()
            for dependency in r1521["evidence_dependencies"].values()
        ),
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
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    handoff = load_handoff()
    _validate_output_path(args.output, args.model, handoff)
    _write_json(args.output, audit_handoff(args.model))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
