"""Fail-closed ER VLCFA C26 stereochemistry curation."""

from __future__ import annotations

import hashlib
import json
import logging
from copy import deepcopy
from math import isclose
from pathlib import Path
from typing import Any


logger = logging.getLogger(__name__)
_REPOSITORY = Path(__file__).resolve().parents[2]
_CURATION_PATH = _REPOSITORY / "data" / "er_vlcfa_3r_stereochemistry.json"
_TRUSTED_CURATION_SHA256 = "1a527ee0bb99b1703f2cc8e18b4af147aa3bde7ca1dcb634d8220d4cb375ab80"
_CHEMICAL_IDENTITY_KEYS = frozenset(
    {
        "chebi", "hmdb", "inchi", "inchikey", "lipidmaps", "lipidmapsm",
        "metacyc.compound", "metacycm", "metanetx.chemical", "seed.compound", "seedm",
    }
)


class ERVLCFAStereochemistryError(RuntimeError):
    """ER VLCFA stereochemistry curation is malformed or does not match the model."""


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _target_contract_digest(curation: dict[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(curation, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    ).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ERVLCFAStereochemistryError(message)


def _validate_curation(curation: dict[str, Any]) -> None:
    _require(isinstance(curation, dict), "ER VLCFA curation must be an object")
    _require(
        set(curation) == {"schema_version", "curation_id", "source_model", "source_sha256", "scope", "neutral_patch", "biochemical_ph_identity", "evidence"},
        "unknown or missing top-level ER VLCFA curation key",
    )
    _require(
        type(curation.get("schema_version")) is int and curation["schema_version"] == 2,
        "unsupported ER VLCFA curation schema",
    )
    _require(isinstance(curation.get("curation_id"), str) and curation["curation_id"], "invalid curation identifier")
    _require(curation.get("source_model") == "data/iyli21.xml", "unexpected ER VLCFA source model")
    _require(
        isinstance(curation.get("source_sha256"), str)
        and len(curation["source_sha256"]) == 64
        and all(character in "0123456789abcdef" for character in curation["source_sha256"]),
        "invalid ER VLCFA source SHA",
    )
    _require(
        _target_contract_digest(curation) == _TRUSTED_CURATION_SHA256,
        "ER VLCFA target contract drifted from its locked SHA",
    )
    neutral = curation.get("neutral_patch")
    _require(isinstance(neutral, dict), "missing ER VLCFA neutral patch")
    metabolite = neutral.get("metabolite")
    _require(isinstance(metabolite, dict), "missing ER VLCFA metabolite contract")
    _require(set(metabolite) == {"id", "legacy_name", "target_name", "allowed_source_tuples", "target_tuple", "annotation_target", "replace_annotation_keys"}, "unknown metabolite contract key")
    for key in ("id", "legacy_name", "target_name"):
        _require(isinstance(metabolite.get(key), str) and metabolite[key], f"invalid metabolite {key}")
    _require(
        metabolite["id"] == "m1446[C_em]"
        and metabolite["legacy_name"] == "(S)-3-hydroxyhexacosanoyl-CoA_C47H86N7O18P3S"
        and metabolite["target_name"] == "(3R)-3-hydroxyhexacosanoyl-CoA",
        "ER VLCFA metabolite identity contract drifted",
    )
    _require(isinstance(metabolite.get("allowed_source_tuples"), list) and metabolite["allowed_source_tuples"], "missing source tuples")
    _require(isinstance(metabolite.get("target_tuple"), dict), "missing target tuple")
    _validate_tuple_records(metabolite["allowed_source_tuples"], "metabolite source")
    _validate_target_tuple(metabolite["target_tuple"], "metabolite")
    _require(
        metabolite["target_tuple"] == {"formula": "C47H86N7O18P3S", "charge": 0},
        "ER VLCFA neutral target tuple drifted",
    )
    _require(metabolite.get("annotation_target") == {"chebi": "CHEBI:76465"}, "invalid neutral annotation target")
    _require(
        isinstance(metabolite.get("replace_annotation_keys"), list)
        and set(metabolite["replace_annotation_keys"]) == _CHEMICAL_IDENTITY_KEYS
        and len(metabolite["replace_annotation_keys"]) == len(_CHEMICAL_IDENTITY_KEYS),
        "invalid annotation replacement keys",
    )
    _require(set(neutral) == {"metabolite", "neutral_closure_metabolites", "reactions"}, "unknown neutral patch key")
    partners = neutral.get("neutral_closure_metabolites")
    _require(isinstance(partners, list) and len(partners) == 7, "missing neutral closure partners")
    _require({partner.get("id") for partner in partners} == {"m377[C_em]", "m1437[C_em]", "m1438[C_em]", "m1439[C_em]", "m1447[C_em]", "m1454[C_em]", "m185[C_pe]"}, "neutral closure partner set drifted")
    for partner in partners:
        _require(set(partner) == {"id", "allowed_source_tuples", "target_tuple"}, "unknown closure partner key")
        _require(isinstance(partner["id"], str) and partner["id"], "invalid closure partner ID")
        _require(isinstance(partner["allowed_source_tuples"], list) and partner["allowed_source_tuples"], "missing closure partner source tuples")
        _validate_tuple_records(partner["allowed_source_tuples"], "closure partner source")
        _validate_target_tuple(partner["target_tuple"], "closure partner")
    _require(isinstance(neutral.get("reactions"), list) and len(neutral["reactions"]) == 2, "expected two ER reactions")
    expected_ids = {"R1419", "R1426"}
    _require({entry.get("id") for entry in neutral["reactions"]} == expected_ids, "wrong ER reaction set")
    for reaction in neutral["reactions"]:
        _require(set(reaction) == {"id", "legacy_name", "target_name", "bounds", "reversible", "stoichiometry", "rhea"}, "unknown reaction contract key")
        _require(isinstance(reaction.get("legacy_name"), str) and isinstance(reaction.get("target_name"), str), "invalid reaction name contract")
        _require(
            isinstance(reaction.get("bounds"), list)
            and len(reaction["bounds"]) == 2
            and all(type(bound) in (int, float) for bound in reaction["bounds"]),
            "invalid reaction bounds",
        )
        _require(isinstance(reaction.get("reversible"), bool), "invalid reaction reversibility")
        _require(
            isinstance(reaction.get("stoichiometry"), dict)
            and reaction["stoichiometry"]
            and all(isinstance(met_id, str) and met_id and type(coefficient) in (int, float) for met_id, coefficient in reaction["stoichiometry"].items()),
            "missing reaction stoichiometry",
        )
        rhea = reaction.get("rhea")
        _require(isinstance(rhea, dict) and set(rhea) == {"master", "directional"}, "invalid directional Rhea contract")
        _require(all(isinstance(rhea[key], str) and rhea[key].isdigit() for key in rhea), "invalid Rhea identifier")
    scope = curation.get("scope")
    _require(
        isinstance(scope, dict)
        and set(scope) == {
            "metabolite_ids",
            "reaction_ids",
            "excluded_metabolite_ids",
            "excluded_identity_source_state",
        }
        and scope["metabolite_ids"] == ["m1446[C_em]"]
        and scope["reaction_ids"] == ["R1419", "R1426"]
        and scope["excluded_metabolite_ids"] == ["m1546[C_pe]"],
        "ER VLCFA scope drifted",
    )
    excluded_state = scope["excluded_identity_source_state"]
    _require(
        isinstance(excluded_state, dict)
        and set(excluded_state) == {"id", "name", "tuple", "absent_annotation_keys"}
        and excluded_state["id"] == "m1546[C_pe]"
        and excluded_state["name"] == "(R)-3-hydroxyhexacosanoyl-CoA_"
        and isinstance(excluded_state["tuple"], dict)
        and isinstance(excluded_state["absent_annotation_keys"], list)
        and set(excluded_state["absent_annotation_keys"]) == _CHEMICAL_IDENTITY_KEYS
        and len(excluded_state["absent_annotation_keys"]) == len(_CHEMICAL_IDENTITY_KEYS),
        "excluded peroxisomal identity source state drifted",
    )
    _validate_tuple_records([excluded_state["tuple"]], "excluded peroxisomal source")
    _require(
        excluded_state["tuple"] == {"formula": None, "charge": 0},
        "excluded peroxisomal source tuple drifted",
    )
    biochemical = curation.get("biochemical_ph_identity")
    _require(
        isinstance(biochemical, dict)
        and set(biochemical) == {"chebi", "biochemical_ph_tuple"}
        and biochemical["chebi"] == "CHEBI:76378"
        and isinstance(biochemical["biochemical_ph_tuple"], dict),
        "invalid biochemical pH identity",
    )
    _validate_target_tuple(biochemical["biochemical_ph_tuple"], "biochemical pH")
    _require(
        biochemical["biochemical_ph_tuple"] == {"formula": "C47H82N7O18P3S", "charge": -4},
        "ER VLCFA biochemical pH tuple drifted",
    )
    _require(isinstance(curation.get("evidence"), str) and curation["evidence"], "missing ER VLCFA evidence")


def _validate_tuple_records(records: list[dict[str, Any]], label: str) -> None:
    for record in records:
        _require(isinstance(record, dict) and set(record) == {"formula", "charge"}, f"invalid {label} tuple")
        _require(record["formula"] is None or isinstance(record["formula"], str), f"invalid {label} formula")
        _require(record["charge"] is None or type(record["charge"]) is int, f"invalid {label} charge")


def _validate_target_tuple(value: dict[str, Any], label: str) -> None:
    _require(
        set(value) == {"formula", "charge"}
        and isinstance(value["formula"], str)
        and value["formula"]
        and type(value["charge"]) is int,
        f"invalid {label} target tuple",
    )


def load_er_vlcfa_stereochemistry_curation(
    curation_path: str | Path | None = None,
) -> dict[str, Any]:
    """Load the single curated source and verify its source/target contracts."""

    path = Path(curation_path) if curation_path is not None else _CURATION_PATH
    try:
        curation = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ERVLCFAStereochemistryError(f"cannot read ER VLCFA curation: {exc}") from exc
    _validate_curation(curation)
    _validate_source_provenance(curation)
    return curation


def _validate_source_provenance(curation: dict[str, Any]) -> None:
    source = _REPOSITORY / curation["source_model"]
    _require(source.is_file(), f"ER VLCFA source model is unavailable: {source}")
    _require(_sha256_file(source) == curation["source_sha256"], "ER VLCFA source model SHA drifted")


def _same_number(left: float, right: float) -> bool:
    return isclose(float(left), float(right), rel_tol=0.0, abs_tol=1e-12)


def _validate_live_reaction(model, contract: dict[str, Any]) -> Any:
    try:
        reaction = model.reactions.get_by_id(contract["id"])
    except KeyError as exc:
        raise ERVLCFAStereochemistryError(f"missing ER reaction {contract['id']}") from exc
    bounds = contract["bounds"]
    _require(
        _same_number(reaction.lower_bound, bounds[0]) and _same_number(reaction.upper_bound, bounds[1]),
        f"{reaction.id} bounds drifted from curated direction",
    )
    _require(bool(reaction.reversibility) is contract["reversible"], f"{reaction.id} reversibility drifted")
    observed = {met.id: float(coefficient) for met, coefficient in reaction.metabolites.items()}
    expected = {met_id: float(coefficient) for met_id, coefficient in contract["stoichiometry"].items()}
    _require(observed == expected, f"{reaction.id} stoichiometry drifted")
    return reaction


def _allowed_tuple(value: tuple[str | None, int | None], allowed: list[dict[str, Any]]) -> bool:
    return any(value == (entry.get("formula"), entry.get("charge")) for entry in allowed)


def _audit_neutral_closure(model, neutral: dict[str, Any]) -> None:
    """Balance the curated neutral closure on a copy; never write partners live."""

    candidate = model.copy()
    met_contract = neutral["metabolite"]
    target = candidate.metabolites.get_by_id(met_contract["id"])
    target.formula = met_contract["target_tuple"]["formula"]
    target.charge = met_contract["target_tuple"]["charge"]
    for partner in neutral["neutral_closure_metabolites"]:
        copied = candidate.metabolites.get_by_id(partner["id"])
        copied.formula = partner["target_tuple"]["formula"]
        copied.charge = partner["target_tuple"]["charge"]
    for contract in neutral["reactions"]:
        reaction = candidate.reactions.get_by_id(contract["id"])
        _require(reaction.check_mass_balance() == {}, f"{reaction.id} neutral mass/charge postcondition failed")


def _validate_live_neutral_closure_tuples(model, neutral: dict[str, Any]) -> None:
    """Reject any uncurated live partner state before counterfactual balancing."""

    for partner in neutral["neutral_closure_metabolites"]:
        try:
            live_partner = model.metabolites.get_by_id(partner["id"])
        except KeyError as exc:
            raise ERVLCFAStereochemistryError(
                f"missing ER closure metabolite {partner['id']}"
            ) from exc
        tuple_value = (live_partner.formula, live_partner.charge)
        target = (partner["target_tuple"]["formula"], partner["target_tuple"]["charge"])
        _require(
            tuple_value == target or _allowed_tuple(tuple_value, partner["allowed_source_tuples"]),
            f"{live_partner.id} has an uncurated source tuple {tuple_value!r}",
        )


def _snapshot(objects: list[Any]) -> list[tuple[Any, str, str | None, int | None, dict]]:
    return [(obj, obj.name, getattr(obj, "formula", None), getattr(obj, "charge", None), deepcopy(obj.annotation or {})) for obj in objects]


def _restore(snapshot: list[tuple[Any, str, str | None, int | None, dict]]) -> None:
    for obj, name, formula, charge, annotation in snapshot:
        obj.name = name
        if hasattr(obj, "formula"):
            obj.formula, obj.charge = formula, charge
        obj.annotation = annotation


def _live_excluded_metabolites(model, curation: dict[str, Any]) -> list[Any]:
    """Resolve curation-excluded copies so their scope is not merely textual."""

    excluded = []
    state = curation["scope"]["excluded_identity_source_state"]
    for metabolite_id in curation["scope"]["excluded_metabolite_ids"]:
        try:
            excluded.append(model.metabolites.get_by_id(metabolite_id))
        except KeyError as exc:
            raise ERVLCFAStereochemistryError(
                f"missing curation-excluded metabolite {metabolite_id}"
            ) from exc
    for metabolite in excluded:
        _require(metabolite.id == state["id"], "excluded peroxisomal ID drifted")
        _require(metabolite.name == state["name"], "excluded peroxisomal identity drifted")
        _require(
            (metabolite.formula, metabolite.charge)
            == (state["tuple"]["formula"], state["tuple"]["charge"]),
            "excluded peroxisomal tuple drifted",
        )
        annotation = metabolite.annotation or {}
        _require(
            all(key not in annotation for key in state["absent_annotation_keys"]),
            "excluded peroxisomal chemical annotation drifted",
        )
    return excluded


def _require_unchanged(snapshot: list[tuple[Any, str, str | None, int | None, dict]], label: str) -> None:
    for object_, name, formula, charge, annotation in snapshot:
        _require(
            object_.name == name
            and getattr(object_, "formula", None) == formula
            and getattr(object_, "charge", None) == charge
            and (object_.annotation or {}) == annotation,
            f"{label} changed during ER VLCFA stereochemistry correction",
        )


def correct_er_vlcfa_3r_stereochemistry(model, curation: dict[str, Any] | None = None) -> dict[str, int]:
    """Apply the neutral identity correction after exact source-contract checks.

    This intentionally does not apply the pH 7.3 tuple; that transaction is
    separately blocked until R1419's 3-oxo partner copies are curated together.
    """

    if curation is None:
        curation = load_er_vlcfa_stereochemistry_curation()
    else:
        _validate_curation(curation)
        _validate_source_provenance(curation)
    neutral = curation["neutral_patch"]
    met_contract = neutral["metabolite"]
    excluded_metabolites = _live_excluded_metabolites(model, curation)
    try:
        metabolite = model.metabolites.get_by_id(met_contract["id"])
    except KeyError as exc:
        raise ERVLCFAStereochemistryError(f"missing ER metabolite {met_contract['id']}") from exc
    reactions = {entry["id"]: _validate_live_reaction(model, entry) for entry in neutral["reactions"]}
    _require(
        metabolite in reactions["R1419"].metabolites and metabolite in reactions["R1426"].metabolites,
        "m1446[C_em] is not shared by the curated R1419/R1426 pair",
    )
    live_tuple = (metabolite.formula, metabolite.charge)
    target_tuple = (met_contract["target_tuple"].get("formula"), met_contract["target_tuple"].get("charge"))
    _require(
        live_tuple == target_tuple or _allowed_tuple(live_tuple, met_contract["allowed_source_tuples"]),
        f"{metabolite.id} has an uncurated source tuple {live_tuple!r}",
    )
    _validate_live_neutral_closure_tuples(model, neutral)
    _audit_neutral_closure(model, neutral)

    changes = {"metabolites": 0, "reactions": 0, "annotations": 0}
    saved = _snapshot([metabolite, *reactions.values()])
    excluded_saved = _snapshot(excluded_metabolites)
    try:
        if metabolite.name == met_contract["legacy_name"]:
            metabolite.name = met_contract["target_name"]
            changes["metabolites"] += 1
        elif metabolite.name != met_contract["target_name"]:
            raise ERVLCFAStereochemistryError(f"{metabolite.id} identity drifted: {metabolite.name!r}")
        if live_tuple != target_tuple:
            metabolite.formula, metabolite.charge = target_tuple
            changes["metabolites"] += 1
        annotation = dict(metabolite.annotation or {})
        corrected = {key: value for key, value in annotation.items() if key not in met_contract["replace_annotation_keys"]}
        corrected.update(met_contract["annotation_target"])
        if corrected != annotation:
            metabolite.annotation = corrected
            changes["annotations"] += 1
        for contract in neutral["reactions"]:
            reaction = reactions[contract["id"]]
            if reaction.name == contract["target_name"]:
                pass
            elif reaction.name == contract["legacy_name"]:
                reaction.name = contract["target_name"]
                changes["reactions"] += 1
            else:
                raise ERVLCFAStereochemistryError(f"{reaction.id} identity drifted: {reaction.name!r}")
            annotation = dict(reaction.annotation or {})
            target = [contract["rhea"]["master"], contract["rhea"]["directional"]]
            corrected = dict(annotation)
            corrected.update({"rh": target, "rhea": target})
            if corrected != annotation:
                reaction.annotation = corrected
                changes["annotations"] += 1
        _verify_live_target_state(model, curation)
        _require_unchanged(excluded_saved, "excluded peroxisomal metabolite")
    except Exception:
        _restore([*saved, *excluded_saved])
        raise

    logger.info("ER VLCFA stereochemistry: m1446[C_em], R1419 and R1426 verified and corrected")
    return changes


def verify_er_vlcfa_3r_stereochemistry_target(model) -> None:
    """Read-only final gate: reject any later pipeline drift before SBML write."""

    curation = load_er_vlcfa_stereochemistry_curation()
    _verify_live_target_state(model, curation)


def _verify_live_target_state(model, curation: dict[str, Any]) -> None:
    """Verify the live target plus a non-mutating neutral chemistry closure."""

    _live_excluded_metabolites(model, curation)
    neutral = curation["neutral_patch"]
    met_contract = neutral["metabolite"]
    try:
        metabolite = model.metabolites.get_by_id(met_contract["id"])
    except KeyError as exc:
        raise ERVLCFAStereochemistryError(f"missing final ER metabolite {met_contract['id']}") from exc
    _require(metabolite.name == met_contract["target_name"], "final ER metabolite identity drift")
    _require(
        (metabolite.formula, metabolite.charge)
        == (met_contract["target_tuple"]["formula"], met_contract["target_tuple"]["charge"]),
        "final ER metabolite tuple drift",
    )
    _require(metabolite.annotation.get("chebi") == met_contract["annotation_target"]["chebi"], "final ER metabolite annotation drift")
    _require(
        all(
            key == "chebi" or key not in metabolite.annotation
            for key in met_contract["replace_annotation_keys"]
        ),
        "final ER metabolite chemical annotation drift",
    )
    for contract in neutral["reactions"]:
        reaction = _validate_live_reaction(model, contract)
        expected = [contract["rhea"]["master"], contract["rhea"]["directional"]]
        _require(reaction.name == contract["target_name"], f"final {reaction.id} identity drift")
        _require(reaction.annotation.get("rhea") == expected and reaction.annotation.get("rh") == expected, f"final {reaction.id} Rhea drift")
    _validate_live_neutral_closure_tuples(model, neutral)
    _audit_neutral_closure(model, neutral)
