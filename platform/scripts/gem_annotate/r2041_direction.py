"""Opt-in forward-only D-arabinokinase (R2041) after the metadata selection; bounds only."""

import json

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r2041_direction.json")


def apply_r2041_direction(model, spec=None):
    """Restore [0, 1000] on the selection's balanced 2 H+ R2041; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "R2041-FORWARD-20261019"
            or spec.get("status") != "user_authorized"
            or spec.get("reaction_id") != "R2041"
            or spec.get("before_bounds") != [-1000, 1000]
            or spec.get("after_bounds") != [0, 1000]):
        raise ValueError("Unexpected R2041 direction curation scope")
    if "R2041" not in model.reactions:
        raise ValueError("R2041 missing; no edits applied")
    reaction = model.reactions.get_by_id("R2041")
    # The species guard also stops this step after an ATP/ADP microspecies migration,
    # which would change the balanced proton count.
    if (reaction.name != spec["reaction_name"]
            or reaction.gene_reaction_rule != spec["gpr"]
            or reaction.annotation.get("ec-code") not in (spec["ec_code"], [spec["ec_code"]])
            or list(reaction.bounds) not in (spec["before_bounds"], spec["after_bounds"])
            or {m.id: c for m, c in reaction.metabolites.items()} != spec["stoichiometry"]
            or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                for m in reaction.metabolites} != spec["species"]):
        raise ValueError("R2041 direction precondition differs; no edits applied")
    if any(k in reaction.notes and reaction.notes[k] != v for k, v in spec["notes"].items()):
        raise ValueError("Conflicting R2041 direction notes; no edits applied")
    # Forward-only without these notes came from somewhere else; the previous-bounds note would be false.
    if (list(reaction.bounds) == spec["after_bounds"]
            and any(reaction.notes.get(k) != v for k, v in spec["notes"].items())):
        raise ValueError("R2041 is already forward-only without this curation; no edits applied")
    notes = {**reaction.notes, **spec["notes"]}
    changed = list(reaction.bounds) != spec["after_bounds"] or notes != reaction.notes
    reaction.bounds = tuple(float(v) for v in spec["after_bounds"])
    reaction.notes = notes
    return {"item": reaction.id, "curation_id": spec["curation_id"],
            "status": "applied" if changed else "already_correct",
            "bounds": list(reaction.bounds), "evidence_status": spec["status"]}
