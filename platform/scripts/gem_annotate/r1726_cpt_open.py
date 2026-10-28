"""Opt-in opening of the cholinephosphotransferase R1726 (Kennedy CDP-choline step); bounds and notes only."""

import json

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r1726_cpt_open.json")


def _state(reaction):
    return ({m.id: c for m, c in reaction.metabolites.items()},
            {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")} for m in reaction.metabolites})


def apply_r1726_cpt_open(model, spec=None):
    """Set [0, 1000] on the inherited (0, 0) R1726 while R1725 stays closed; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    kept = spec.get("kept_closed", {})
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "R1726-CPT-OPEN-20261027"
            or spec.get("status") != "user_authorized"
            or spec.get("reaction_id") != "R1726"
            or spec.get("before_bounds") != [0, 0]
            or spec.get("after_bounds") != [0, 1000]
            or kept.get("reaction_id") != "R1725"
            or kept.get("bounds") != [0, 0]):
        raise ValueError("Unexpected R1726 opening curation scope")
    if "R1726" not in model.reactions or "R1725" not in model.reactions:
        raise ValueError("R1726 or R1725 missing; no edits applied")
    reaction = model.reactions.get_by_id("R1726")
    stoichiometry, species = _state(reaction)
    if (reaction.name != spec["reaction_name"]
            or reaction.gene_reaction_rule != spec["gpr"]
            or reaction.annotation.get("ec-code") not in (spec["ec_code"], [spec["ec_code"]])
            or list(reaction.bounds) not in (spec["before_bounds"], spec["after_bounds"])
            or stoichiometry != spec["stoichiometry"]
            or species != spec["species"]):
        raise ValueError("R1726 opening precondition differs; no edits applied")
    # The decision covers R1726 alone: an opened R1725 (sphingoid-lyase PE route) was never reviewed with it.
    closed = model.reactions.get_by_id("R1725")
    if (closed.name != kept["reaction_name"] or closed.gene_reaction_rule != kept["gpr"]
            or list(closed.bounds) != kept["bounds"]
            or _state(closed) != (kept["stoichiometry"], kept["species"])):
        raise ValueError("R1725 is not the closed reaction this curation reviewed; no edits applied")
    if any(k in reaction.notes and reaction.notes[k] != v for k, v in spec["notes"].items()):
        raise ValueError("Conflicting R1726 opening notes; no edits applied")
    # Closed again with these notes means something re-closed it after this curation; do not undo that.
    if list(reaction.bounds) == spec["before_bounds"] and any(k in reaction.notes for k in spec["notes"]):
        raise ValueError("R1726 is closed but carries this curation's notes; no edits applied")
    # Open without these notes came from somewhere else; the previous-bounds note would be false.
    if (list(reaction.bounds) == spec["after_bounds"]
            and any(reaction.notes.get(k) != v for k, v in spec["notes"].items())):
        raise ValueError("R1726 is already open without this curation; no edits applied")
    notes = {**reaction.notes, **spec["notes"]}
    changed = list(reaction.bounds) != spec["after_bounds"] or notes != reaction.notes
    reaction.bounds = tuple(float(v) for v in spec["after_bounds"])
    reaction.notes = notes
    return {"item": reaction.id, "curation_id": spec["curation_id"],
            "status": "applied" if changed else "already_correct",
            "bounds": list(reaction.bounds), "kept_closed": closed.id, "evidence_status": spec["status"]}
