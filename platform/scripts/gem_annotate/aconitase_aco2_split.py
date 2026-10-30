"""Opt-in GPR split of the aconitase reactions from the ACO2-type paralog.

R1556 (cytosol), R216 and R238 (mitochondrion) read `YALI1D11984g or YALI1E17935g`, as inherited from iYli21; R239
(cytosol) already reads `YALI1D11984g`. The three get the ACO1-type `YALI1D11984g` alone; the ACO2-type
`YALI1E17935g` (homocitrate-specific or inactive by homology to S. cerevisiae Aco2p) is kept in the model with no
reaction. The decision and its evidence are data in ``aconitase_aco2_split.json``. Every check on all three reactions
and both genes runs before any edit; a partly applied state is rejected.
"""

import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("aconitase_aco2_split.json")
BEFORE = "YALI1D11984g or YALI1E17935g"
AFTER = "YALI1D11984g"
TARGETS = ("R1556", "R216", "R238")
GENES = {"YALI1D11984g", "YALI1E17935g"}
UNLINKED = "YALI1E17935g"


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _check_spec(spec):
    rows, genes = spec.get("reactions", {}), spec.get("genes", {})
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "ACONITASE-ACO2-SPLIT-20261029"
            or spec.get("status") != "user_authorized_gpr_correction"
            or set(rows) != set(TARGETS)
            or any(rows[r].get("before_gpr") != BEFORE or rows[r].get("after_gpr") != AFTER for r in TARGETS)
            or any(not isinstance(rows[r].get("before_notes"), dict) or not rows[r].get("notes") for r in TARGETS)
            or set(genes) != GENES
            or spec.get("unlinked_gene", {}).get("id") != UNLINKED
            or genes[UNLINKED].get("reactions_before") != sorted(TARGETS)
            or genes[UNLINKED].get("reactions_after") != []
            or genes[AFTER].get("reactions_after") != genes[AFTER].get("reactions_before")):
        raise ValueError("Unexpected aconitase ACO2 split scope/status")


def apply_aconitase_aco2_split(model, spec=None):
    """Assign R1556, R216 and R238 to the ACO1-type aconitase alone; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    done = {}
    for rid, row in spec["reactions"].items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        gpr = reaction.gene_reaction_rule
        if (reaction.name != row["reaction_name"]
                or list(reaction.bounds) != row["bounds"]
                or {m.id: c for m, c in reaction.metabolites.items()} != row["stoichiometry"]
                or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                    for m in reaction.metabolites} != row["species"]
                or reaction.annotation != row["annotation"]
                or gpr not in (row["before_gpr"], row["after_gpr"])):
            raise ValueError(f"{rid} aconitase split precondition differs; no edits applied")
        notes, current = row["notes"], _notes(reaction.notes)
        if any(k in current and current[k] != v for k, v in notes.items()):
            raise ValueError(f"Conflicting {rid} aconitase split notes; no edits applied")
        if {k: v for k, v in current.items() if k not in notes} != row["before_notes"]:
            raise ValueError(f"{rid} notes outside this curation differ; no edits applied")
        is_done = gpr == row["after_gpr"]
        if is_done and any(current.get(k) != v for k, v in notes.items()):
            raise ValueError(f"{rid} already carries this GPR without this curation; no edits applied")
        if not is_done and any(k in current for k in notes):
            raise ValueError(f"{rid} carries this curation's notes without its GPR; no edits applied")
        done[rid] = is_done
    if len(set(done.values())) != 1:
        raise ValueError(f"Aconitase split partly applied ({done}); no edits applied")
    applied = all(done.values())
    for gene_id, row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        expected = row["reactions_after"] if applied else row["reactions_before"]
        if gene.annotation != row["annotation"] or sorted(r.id for r in gene.reactions) != expected:
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    record = {"item": sorted(TARGETS), "assignment_id": spec["assignment_id"], "unlinked_gene": UNLINKED,
              "evidence_status": spec["status"]}
    if applied:
        return {**record, "status": "already_correct", "gpr": AFTER}
    for rid, row in spec["reactions"].items():
        reaction = model.reactions.get_by_id(rid)
        reaction.gene_reaction_rule = row["after_gpr"]  # history-tracked by cobra
        context = get_context(reaction)
        if context:
            # Notes are not history-tracked by cobra; restore them on context exit like the GPR.
            context(partial(setattr, reaction, "notes", reaction.notes))
        reaction.notes = {**reaction.notes, **row["notes"]}
    if any(sorted(r.id for r in model.genes.get_by_id(g).reactions) != spec["genes"][g]["reactions_after"]
           for g in GENES):
        raise RuntimeError("Aconitase split did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": BEFORE, "gpr": AFTER}
