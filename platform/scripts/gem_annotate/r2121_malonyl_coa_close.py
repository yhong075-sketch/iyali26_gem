"""Opt-in closure of R2121 (malonate semialdehyde + CoA + NADP+ to malonyl-CoA + NADPH); bounds and notes only.

R2121 is KEGG R00740 (EC 1.2.1.75), the malonyl-CoA reductase of bacterial and archaeal CO2-fixation cycles. It came
from iYli21 on YALI1C02550g, whose KEGG ortholog group (K00140) catalyses R00705, kept in the model as R2119. In WT
pFBA on e584151f it made all cytosolic malonyl-CoA, so the acetyl-CoA carboxylase R88 carried none. The decision and
its evidence are data in ``r2121_malonyl_coa_close.json``. Every check runs before the edit.
"""

import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r2121_malonyl_coa_close.json")
CURATION_ID = "R2121-MALONYL-COA-CLOSE-20261112"


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _state(reaction):
    return ({m.id: c for m, c in reaction.metabolites.items()},
            {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")} for m in reaction.metabolites})


REQUIRED = ("reaction_name", "gpr", "stoichiometry", "species", "annotation")


def _check_spec(spec):
    kept, notes, before = spec.get("retained_reaction", {}), spec.get("notes"), spec.get("before_notes")
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != CURATION_ID
            or spec.get("status") != "user_authorized"
            or spec.get("reaction_id") != "R2121"
            or spec.get("before_bounds") != [-1000, 1000]
            or spec.get("after_bounds") != [0, 0]
            or any(key not in spec for key in REQUIRED)
            or kept.get("reaction_id") != "R2119"
            or kept.get("gpr") != spec.get("gpr")
            or kept.get("bounds") != [0, 1000]
            or any(key not in kept for key in REQUIRED)
            or not isinstance(notes, dict) or not notes
            or notes.get("r2121_previous_bounds") != str(spec["before_bounds"])
            or not isinstance(before, dict) or set(notes) & set(before)):
        raise ValueError("Unexpected R2121 closure curation scope")


def apply_r2121_malonyl_coa_close(model, spec=None):
    """Set (0, 0) on R2121 while R2119 keeps the gene's decarboxylating reaction; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    kept = spec["retained_reaction"]
    if "R2121" not in model.reactions or "R2119" not in model.reactions:
        raise ValueError("R2121 or R2119 missing; no edits applied")
    reaction = model.reactions.get_by_id("R2121")
    if (reaction.name != spec["reaction_name"]
            or reaction.gene_reaction_rule != spec["gpr"]
            or list(reaction.bounds) not in (spec["before_bounds"], spec["after_bounds"])
            or _state(reaction) != (spec["stoichiometry"], spec["species"])
            or reaction.annotation != spec["annotation"]):
        raise ValueError("R2121 closure precondition differs; no edits applied")
    # The case for closing R2121 rests on the gene keeping its real reaction, R2119 (R00705), as reviewed.
    retained = model.reactions.get_by_id("R2119")
    if (retained.name != kept["reaction_name"] or retained.gene_reaction_rule != kept["gpr"]
            or list(retained.bounds) != kept["bounds"]
            or _state(retained) != (kept["stoichiometry"], kept["species"])
            or retained.annotation != kept["annotation"]):
        raise ValueError("R2119 is not the retained reaction this curation reviewed; no edits applied")
    notes, current = spec["notes"], _notes(reaction.notes)
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError("Conflicting R2121 closure notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != spec["before_notes"]:
        raise ValueError("R2121 notes outside this curation differ; no edits applied")
    closed = list(reaction.bounds) == spec["after_bounds"]
    if closed and any(current.get(k) != v for k, v in notes.items()):
        # Closed by something else: the previous-bounds note would be false.
        raise ValueError("R2121 is already closed without this curation; no edits applied")
    if not closed and any(k in current for k in notes):
        # Open again with these notes means something reopened it after this curation; do not undo that.
        raise ValueError("R2121 is open but carries this curation's notes; no edits applied")
    record = {"item": reaction.id, "curation_id": spec["curation_id"], "retained": retained.id,
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
        raise RuntimeError("R2121 did not reach the curated bounds")
    return {**record, "status": "applied", "bounds_before": spec["before_bounds"], "bounds": list(reaction.bounds)}
