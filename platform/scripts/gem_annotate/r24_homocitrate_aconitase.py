"""Opt-in GPR correction of R24 (mitochondrial homocitrate dehydration) from LYS4 to the two aconitases.

R24 reads `YALI1E03387g`, the LYS4-like homoaconitase whose reviewed activity is the R395 step only; yeast and
fungal homoaconitases do not use homocitrate, which aconitase dehydrates. R24 gets `YALI1D11984g or YALI1E17935g`
(ACO1- and ACO2-type, either sufficient in S. cerevisiae); `YALI1E03387g` stays on R395. The decision and its
evidence are data in ``r24_homocitrate_aconitase.json``. The gene guards pin the state after
``apply_aconitase_aco2_split`` (`YALI1E17935g` with no reaction), so this step runs after it. Every check runs
before any edit.
"""

import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r24_homocitrate_aconitase.json")
TARGET = "R24"
BEFORE = "YALI1E03387g"
AFTER = "YALI1D11984g or YALI1E17935g"
GENES = {"YALI1E03387g", "YALI1D11984g", "YALI1E17935g"}


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _check_spec(spec):
    rows, genes = spec.get("reactions", {}), spec.get("genes", {})
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "R24-HOMOCITRATE-ACO-20261030"
            or spec.get("status") != "user_authorized_gpr_correction"
            or set(rows) != {TARGET}
            or rows[TARGET].get("before_gpr") != BEFORE or rows[TARGET].get("after_gpr") != AFTER
            or not isinstance(rows[TARGET].get("before_notes"), dict) or not rows[TARGET].get("notes")
            or set(genes) != GENES
            or genes["YALI1E03387g"].get("reactions_before") != ["R24", "R395"]
            or genes["YALI1E03387g"].get("reactions_after") != ["R395"]
            or genes["YALI1E17935g"].get("reactions_before") != []
            or genes["YALI1E17935g"].get("reactions_after") != [TARGET]
            or genes["YALI1D11984g"].get("reactions_after")
            != sorted([*genes["YALI1D11984g"].get("reactions_before", []), TARGET])):
        raise ValueError("Unexpected R24 homocitrate aconitase scope/status")


def apply_r24_homocitrate_aconitase(model, spec=None):
    """Assign R24 to the ACO1- or ACO2-type aconitase; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    row = spec["reactions"][TARGET]
    if TARGET not in model.reactions:
        raise ValueError(f"{TARGET} missing; no edits applied")
    reaction = model.reactions.get_by_id(TARGET)
    gpr = reaction.gene_reaction_rule
    if (reaction.name != row["reaction_name"]
            or list(reaction.bounds) != row["bounds"]
            or {m.id: c for m, c in reaction.metabolites.items()} != row["stoichiometry"]
            or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                for m in reaction.metabolites} != row["species"]
            or reaction.annotation != row["annotation"]
            or gpr not in (row["before_gpr"], row["after_gpr"])):
        raise ValueError(f"{TARGET} homocitrate aconitase precondition differs; no edits applied")
    notes, current = row["notes"], _notes(reaction.notes)
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError(f"Conflicting {TARGET} homocitrate aconitase notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != row["before_notes"]:
        raise ValueError(f"{TARGET} notes outside this curation differ; no edits applied")
    applied = gpr == row["after_gpr"]
    if applied and any(current.get(k) != v for k, v in notes.items()):
        raise ValueError(f"{TARGET} already carries this GPR without this curation; no edits applied")
    if not applied and any(k in current for k in notes):
        raise ValueError(f"{TARGET} carries this curation's notes without its GPR; no edits applied")
    for gene_id, gene_row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        expected = gene_row["reactions_after"] if applied else gene_row["reactions_before"]
        if gene.annotation != gene_row["annotation"] or sorted(r.id for r in gene.reactions) != expected:
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    record = {"item": [TARGET], "assignment_id": spec["assignment_id"], "evidence_status": spec["status"]}
    if applied:
        return {**record, "status": "already_correct", "gpr": AFTER}
    reaction.gene_reaction_rule = row["after_gpr"]  # history-tracked by cobra
    context = get_context(reaction)
    if context:
        # Notes are not history-tracked by cobra; restore them on context exit like the GPR.
        context(partial(setattr, reaction, "notes", reaction.notes))
    reaction.notes = {**reaction.notes, **notes}
    if any(sorted(r.id for r in model.genes.get_by_id(g).reactions) != spec["genes"][g]["reactions_after"]
           for g in GENES):
        raise RuntimeError("R24 homocitrate aconitase did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": BEFORE, "gpr": AFTER}
