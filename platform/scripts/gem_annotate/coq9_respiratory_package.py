"""User-authorized 2026-10-08 CoQ9/respiratory-chain package, applied on the E5 + qcycle stage.

Every item is data in ``coq9_respiratory_package.json``. Each step checks the exact current
definition of its targets (stoichiometry, bounds, Boolean GPR, species identity) against the
recorded ``before`` or ``after`` state and fails closed on anything else, so a rerun is a no-op.
Genes are never removed; metabolites are removed only when the curated file lists them and
they have become orphans.
"""

import json

from cobra import Reaction
from cobra.core.gene import GPR

from .coq9 import apply_coq9_functional_gpr, boolean_key, exact_residual
from .config import MODEL

PACKAGE_PATH = MODEL.curation_file("coq9_respiratory_package.json")
NOTE = "coq9_respiratory_package"


def _stoichiometry(reaction):
    return {m.id: float(c) for m, c in reaction.metabolites.items()}


def _matches(reaction, state):
    return (_stoichiometry(reaction) == state["stoichiometry"]
            and list(reaction.bounds) == state["bounds"]
            and boolean_key(reaction.gpr.body) == boolean_key(GPR.from_string(state["gpr"]).body))


def _check_species(model, species):
    for mid, identity in species.items():
        if mid not in model.metabolites:
            raise ValueError(f"Package species missing: {mid}")
        met = model.metabolites.get_by_id(mid)
        if {k: getattr(met, k) for k in ("formula", "charge", "compartment")} != identity:
            raise ValueError(f"Package species identity differs: {mid}")


def _check_genes(model, genes):
    for gid, identity in genes.items():
        if gid not in model.genes:
            raise ValueError(f"Package gene missing: {gid}")
        gene = model.genes.get_by_id(gid)
        for key, value in identity.items():
            if gene.annotation.get(key) not in (value, [value]):
                raise ValueError(f"Package gene identity differs: {gid} {key}")


def _state(model, rid, before, after):
    """Return 'before', 'after' or 'absent'; raise on any other definition."""
    if rid not in model.reactions:
        if after is not None:
            raise ValueError(f"Package target reaction missing: {rid}")
        return "absent"
    reaction = model.reactions.get_by_id(rid)
    if _matches(reaction, before):
        return "before"
    if after is not None and _matches(reaction, after):
        return "after"
    raise ValueError(f"{rid}: definition differs from the curated before/after state; no edits applied")


def _set(model, rid, after, note):
    reaction = model.reactions.get_by_id(rid)
    reaction.subtract_metabolites(dict(reaction.metabolites))
    reaction.add_metabolites({model.metabolites.get_by_id(mid): c
                              for mid, c in after["stoichiometry"].items()})
    reaction.bounds = tuple(after["bounds"])
    reaction.gene_reaction_rule = after["gpr"]
    reaction.notes = {**reaction.notes, NOTE: note}


def _check_balance(model, rid, stoichiometry):
    probe = Reaction(f"{rid}_balance_probe")
    probe.add_metabolites({model.metabolites.get_by_id(mid).copy(): c
                           for mid, c in stoichiometry.items()})
    if exact_residual(probe):
        raise ValueError(f"{rid}: curated target is not element/charge balanced")


def _preflight(model, item):
    """Validate every target of one item before any edit."""
    states = {}
    for rid, row in item.get("set", {}).items():
        _check_species(model, row["after_species"])
        if set(row["after_species"]) != set(row["after"]["stoichiometry"]):
            raise ValueError(f"{rid}: curated species list does not match the target equation")
        _check_balance(model, rid, row["after"]["stoichiometry"])
        states[rid] = _state(model, rid, row["before"], row["after"])
    for rid, before in item.get("remove", {}).items():
        states[rid] = _state(model, rid, before, None)
    _check_genes(model, item.get("genes", {}))
    for mid in item.get("remove_orphan_metabolites", []):
        if mid not in model.metabolites:
            continue
        # After this item, only its removed reactions or rewritten targets may have used it.
        releasing = set(item.get("remove", {})) | {
            rid for rid, row in item.get("set", {}).items() if mid not in row["after"]["stoichiometry"]}
        if not {r.id for r in model.metabolites.get_by_id(mid).reactions} <= releasing:
            raise ValueError(f"{mid} would not be orphaned by item {item['id']}; no edits applied")
    for rid, row in item.get("set_notes", {}).items():
        if rid not in model.reactions or model.reactions.get_by_id(rid).notes.get(row["key"]) not in (
                row["before"], row["after"]):
            raise ValueError(f"{rid} note {row['key']} differs from the curated before/after text")
    return states


def _apply_item(model, item):
    states = _preflight(model, item)
    for rid, row in item.get("set", {}).items():
        if states[rid] == "before":
            _set(model, rid, row["after"], item["note"])
    for rid, row in item.get("set_notes", {}).items():
        model.reactions.get_by_id(rid).notes[row["key"]] = row["after"]
    removed = [rid for rid, state in states.items() if rid in item.get("remove", {}) and state == "before"]
    model.remove_reactions([model.reactions.get_by_id(rid) for rid in removed], remove_orphans=False)
    orphans = []
    for mid in item.get("remove_orphan_metabolites", []):
        if mid not in model.metabolites:
            continue
        met = model.metabolites.get_by_id(mid)
        if met.reactions:
            raise ValueError(f"{mid} is not orphaned; refusing to remove it")
        orphans.append(met)
    model.remove_metabolites(orphans)
    changed = bool(removed or orphans or any(s == "before" for rid, s in states.items()
                                             if rid in item.get("set", {})))
    return {"item": item["id"], "status": "applied" if changed else "already_correct",
            "set": sorted(item.get("set", {})), "removed_reactions": removed,
            "removed_metabolites": [m.id for m in orphans], "evidence_tier": item["evidence_tier"]}


def _check_duplicate(model, item):
    """C1: the normalized reaction must equal its retained counterpart exactly."""
    dup = item["duplicate"]
    rid, keep = dup["reaction"], dup["duplicate_of"]
    if rid not in model.reactions:
        return
    reaction, counterpart = model.reactions.get_by_id(rid), model.reactions.get_by_id(keep)
    normalized = {}
    for mid, c in _stoichiometry(reaction).items():
        target = dup["species_map"].get(mid, mid)
        normalized[target] = normalized.get(target, 0.0) + c
    for mid, c in dup["add"].items():
        normalized[mid] = normalized.get(mid, 0.0) + c
    if (normalized != _stoichiometry(counterpart) or list(counterpart.bounds) != dup["normalized_bounds"]
            or not {g.id for g in reaction.genes} <= {g.id for g in counterpart.genes}):
        raise ValueError(f"{rid} is not an exact duplicate of {keep} after normalization")


def apply_coq9_respiratory_package(model, spec=None):
    spec = json.loads(PACKAGE_PATH.read_text()) if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("package_id") != "COQ9-RESP-PKG-20261008"
            or [item["id"] for item in spec["items"]] != ["A1", "A2", "B1", "B3", "C1", "C3", "C4"]):
        raise ValueError("Unexpected CoQ9/respiratory package scope")
    for item in spec["items"]:
        if "duplicate" in item:
            _check_duplicate(model, item)
        if item["id"] != "A2":
            _preflight(model, item)
    # A2 runs first: its own guard checks R695 before mutating, and no other item touches R695,
    # so a stale R695 stops the package before any edit.
    a2 = next(item for item in spec["items"] if item["id"] == "A2")
    records = {"A2": {**apply_coq9_functional_gpr(model), "item": "A2", "evidence_tier": a2["evidence_tier"]}}
    for item in spec["items"]:
        if item["id"] != "A2":
            records[item["id"]] = _apply_item(model, item)
    return {"package_id": spec["package_id"], "status": spec["status"],
            "items": [records[item["id"]] for item in spec["items"]]}
