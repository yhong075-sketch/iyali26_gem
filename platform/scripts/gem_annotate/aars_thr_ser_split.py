"""Opt-in GPR split of the cytosolic threonyl- and seryl-tRNA synthetases from their mitochondrial paralogs.

R759 (cytosol) and R760 (mitochondrion) both read `YALI1A09121g or YALI1F14882g`, and R723 (cytosol) reads
`YALI1B07957g or YALI1F03908g`, as inherited from iYli21. Each cytosolic reaction gets its cytosolic enzyme, R760
its mitochondrial one; YALI1B07957g (DIA4-like) is kept in the model with no reaction. The decision and its
evidence are data in ``aars_thr_ser_split.json``. Every check on all three reactions and four genes runs before
any edit; a partly applied state is rejected.
"""

import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("aars_thr_ser_split.json")
AFTER = {"R759": "YALI1F14882g", "R760": "YALI1A09121g", "R723": "YALI1F03908g"}
BEFORE = {"R759": "YALI1A09121g or YALI1F14882g", "R760": "YALI1A09121g or YALI1F14882g",
          "R723": "YALI1B07957g or YALI1F03908g"}
GENES = {"YALI1A09121g", "YALI1F14882g", "YALI1B07957g", "YALI1F03908g"}
UNLINKED = "YALI1B07957g"


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _check_spec(spec):
    rows = spec.get("reactions", {})
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "AARS-THR-SER-SPLIT-20261029"
            or spec.get("status") != "user_authorized_gpr_correction"
            or set(rows) != set(AFTER)
            or any(rows[r].get("before_gpr") != BEFORE[r] or rows[r].get("after_gpr") != AFTER[r] for r in AFTER)
            or any(not isinstance(rows[r].get("before_notes"), dict) or not rows[r].get("notes") for r in AFTER)
            or set(spec.get("genes", {})) != GENES
            or spec.get("unlinked_gene", {}).get("id") != UNLINKED
            or spec["genes"][UNLINKED].get("reactions_after") != []):
        raise ValueError("Unexpected aaRS Thr/Ser split scope/status")


def apply_aars_thr_ser_split(model, spec=None):
    """Assign R759/R723 to the cytosolic enzymes and R760 to the mitochondrial one; fail closed otherwise."""
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
            raise ValueError(f"{rid} aaRS split precondition differs; no edits applied")
        notes, current = row["notes"], _notes(reaction.notes)
        if any(k in current and current[k] != v for k, v in notes.items()):
            raise ValueError(f"Conflicting {rid} aaRS split notes; no edits applied")
        if {k: v for k, v in current.items() if k not in notes} != row["before_notes"]:
            raise ValueError(f"{rid} notes outside this curation differ; no edits applied")
        is_done = gpr == row["after_gpr"]
        if is_done and any(current.get(k) != v for k, v in notes.items()):
            raise ValueError(f"{rid} already carries this GPR without this curation; no edits applied")
        if not is_done and any(k in current for k in notes):
            raise ValueError(f"{rid} carries this curation's notes without its GPR; no edits applied")
        done[rid] = is_done
    if len(set(done.values())) != 1:
        raise ValueError(f"aaRS split partly applied ({done}); no edits applied")
    applied = all(done.values())
    for gene_id, row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        expected = row["reactions_after"] if applied else row["reactions_before"]
        if gene.annotation != row["annotation"] or sorted(r.id for r in gene.reactions) != expected:
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    record = {"item": sorted(AFTER), "assignment_id": spec["assignment_id"], "unlinked_gene": UNLINKED,
              "evidence_status": spec["status"]}
    if applied:
        return {**record, "status": "already_correct", "gpr": dict(AFTER)}
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
        raise RuntimeError("aaRS split did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": dict(BEFORE), "gpr": dict(AFTER)}
