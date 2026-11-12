"""Opt-in closure of R1172 (gene-less HMG-CoA transport between cytosol and mitochondrion); bounds and notes only.

R1172 has no GPR and no known carrier. In WT pFBA on e584151f and 2f352037 it carried mitochondrial HMG-CoA (made by
R83 and R412) to the cytosolic HMG-CoA reductase R410, so the cytosolic route R82 (ERG10-like) and R411 (ERG13-like)
carried none. The user reversed the 2026-09-09 retention for its flux on 2026-11-12: the reaction is kept, closed.
The decision and its evidence are data in ``r1172_hmgcoa_transport_close.json``. Every check runs before the edit.
"""

import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r1172_hmgcoa_transport_close.json")
CURATION_ID = "R1172-HMGCOA-TRANSPORT-CLOSE-20261112"
ROUTE = ("R82", "R411", "R410")


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _state(reaction):
    return ({m.id: c for m, c in reaction.metabolites.items()},
            {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")} for m in reaction.metabolites})


REQUIRED = ("reaction_name", "gpr", "stoichiometry", "species", "annotation")


def _check_spec(spec):
    route, notes, before = spec.get("required_route"), spec.get("notes"), spec.get("before_notes")
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != CURATION_ID
            or spec.get("status") != "user_authorized"
            or spec.get("reaction_id") != "R1172"
            or spec.get("gpr") != ""
            or spec.get("before_bounds") != [-1000, 1000]
            or spec.get("after_bounds") != [0, 0]
            or any(key not in spec for key in REQUIRED)
            or not isinstance(route, list)
            or tuple(row.get("reaction_id") for row in route) != ROUTE
            or any(key not in row for row in route for key in (*REQUIRED, "bounds"))
            or any(not row["gpr"] for row in route)
            or not isinstance(notes, dict) or not notes
            or notes.get("r1172_previous_bounds") != str(spec["before_bounds"])
            or not isinstance(before, dict) or set(notes) & set(before)):
        raise ValueError("Unexpected R1172 closure curation scope")


def apply_r1172_hmgcoa_transport_close(model, spec=None):
    """Set (0, 0) on R1172 while the cytosolic R82 / R411 / R410 route stays as reviewed; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    if any(rid not in model.reactions for rid in ("R1172", *ROUTE)):
        raise ValueError("R1172 or a reaction of the cytosolic route missing; no edits applied")
    reaction = model.reactions.get_by_id("R1172")
    if (reaction.name != spec["reaction_name"]
            or reaction.gene_reaction_rule != spec["gpr"]
            or list(reaction.bounds) not in (spec["before_bounds"], spec["after_bounds"])
            or _state(reaction) != (spec["stoichiometry"], spec["species"])
            or reaction.annotation != spec["annotation"]):
        raise ValueError("R1172 closure precondition differs; no edits applied")
    # The case for closing R1172 rests on the cytosolic route being there, as reviewed.
    for row in spec["required_route"]:
        kept = model.reactions.get_by_id(row["reaction_id"])
        if (kept.name != row["reaction_name"] or kept.gene_reaction_rule != row["gpr"]
                or list(kept.bounds) != row["bounds"]
                or _state(kept) != (row["stoichiometry"], row["species"])
                or kept.annotation != row["annotation"]):
            raise ValueError(f"{kept.id} is not the cytosolic route reaction this curation reviewed; no edits applied")
    notes, current = spec["notes"], _notes(reaction.notes)
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError("Conflicting R1172 closure notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != spec["before_notes"]:
        raise ValueError("R1172 notes outside this curation differ; no edits applied")
    closed = list(reaction.bounds) == spec["after_bounds"]
    if closed and any(current.get(k) != v for k, v in notes.items()):
        # Closed by something else: the previous-bounds note would be false.
        raise ValueError("R1172 is already closed without this curation; no edits applied")
    if not closed and any(k in current for k in notes):
        # Open again with these notes means something reopened it after this curation; do not undo that.
        raise ValueError("R1172 is open but carries this curation's notes; no edits applied")
    record = {"item": reaction.id, "curation_id": spec["curation_id"], "route": list(ROUTE),
              "evidence_status": spec["status"]}
    if closed:
        return {**record, "status": "already_correct", "bounds": list(reaction.bounds)}
    reaction.bounds = tuple(float(v) for v in spec["after_bounds"])  # history-tracked by cobra
    context = get_context(reaction)
    if context:
        # Notes are not history-tracked by cobra; restore them on context exit like the bounds.
        context(partial(setattr, reaction, "notes", reaction.notes))
    reaction.notes = {**reaction.notes, **notes}
    if list(reaction.bounds) != spec["after_bounds"]:
        raise RuntimeError("R1172 did not reach the curated bounds")
    return {**record, "status": "applied", "bounds_before": spec["before_bounds"], "bounds": list(reaction.bounds)}
