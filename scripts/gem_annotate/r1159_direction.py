"""Conditional Golgi-to-cytosol passive proton leak; not a W29 measurement."""

import json

from .config import CURATION_DATA_DIR

CURATION_PATH = CURATION_DATA_DIR / "r1159_direction.json"


def apply_r1159_direction(model, spec=None):
    """Apply only to the reviewed proton-only reaction after field selection."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "R1159-GOLGI-LEAK-20260923"
            or spec.get("status") != "user_authorized_conditional"
            or spec.get("reaction_id") != "R1159"
            or spec.get("before_bounds") != [-1000, 1000]
            or spec.get("after_bounds") != [-1000, 0]):
        raise ValueError("Unexpected R1159 direction curation scope")
    reaction = model.reactions.get_by_id("R1159")
    if (reaction.name != spec["reaction_name"]
            or reaction.gene_reaction_rule != ""
            or list(reaction.bounds) not in (spec["before_bounds"], spec["after_bounds"])
            or {m.id: c for m, c in reaction.metabolites.items()} != spec["stoichiometry"]
            or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                for m in reaction.metabolites} != spec["species"]):
        raise ValueError("R1159 direction precondition differs; no edits applied")
    if any(k in reaction.notes and reaction.notes[k] != v for k, v in spec["notes"].items()):
        raise ValueError("Conflicting R1159 direction notes; no edits applied")
    notes = {**reaction.notes, **spec["notes"]}
    changed = list(reaction.bounds) != spec["after_bounds"] or notes != reaction.notes
    reaction.bounds = tuple(float(v) for v in spec["after_bounds"])
    reaction.notes = notes
    return {"item": reaction.id, "curation_id": spec["curation_id"],
            "status": "applied" if changed else "already_correct",
            "bounds": list(reaction.bounds), "evidence_status": spec["status"]}
