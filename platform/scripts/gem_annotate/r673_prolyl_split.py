"""Opt-in GPR split of the cytosolic prolyl-tRNA synthetase from its mitochondrial-type paralog.

R673 (cytosol) reads `YALI1D19110g or YALI1E05929g`, as inherited from iYli21. It gets the YHR020W-like
`YALI1E05929g` (TargetP / DeepLoc: no mitochondrial signal); the AIM10-like `YALI1D19110g` (predicted mitochondrial)
is kept in the model with no reaction, since the model has no mitochondrial prolyl-tRNA reaction. The decision and its
evidence are data in ``r673_prolyl_split.json``. Every check runs before the edit.
"""

import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r673_prolyl_split.json")
BEFORE = "YALI1D19110g or YALI1E05929g"
AFTER = "YALI1E05929g"
UNLINKED = "YALI1D19110g"
GENES = {"YALI1D19110g", "YALI1E05929g"}


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _check_spec(spec):
    row, genes = spec.get("reaction", {}), spec.get("genes", {})
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "R673-PROLYL-SPLIT-20261029"
            or spec.get("status") != "user_authorized_gpr_correction"
            or row.get("id") != "R673"
            or row.get("before_gpr") != BEFORE or row.get("after_gpr") != AFTER
            or not isinstance(row.get("before_notes"), dict) or not row.get("notes")
            or set(genes) != GENES
            or spec.get("unlinked_gene", {}).get("id") != UNLINKED
            or any(genes[g].get("reactions_before") != ["R673"] for g in genes)
            or any(genes[g].get("reactions_after") != ([] if g == UNLINKED else ["R673"]) for g in genes)):
        raise ValueError("Unexpected R673 prolyl split scope/status")


def apply_r673_prolyl_split(model, spec=None):
    """Assign R673 to the cytosolic prolyl-tRNA synthetase and unlink the AIM10-like paralog; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    row = spec["reaction"]
    if "R673" not in model.reactions:
        raise ValueError("R673 missing; no edits applied")
    reaction = model.reactions.R673
    gpr = reaction.gene_reaction_rule
    if (reaction.name != row["reaction_name"]
            or list(reaction.bounds) != row["bounds"]
            or {m.id: c for m, c in reaction.metabolites.items()} != row["stoichiometry"]
            or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                for m in reaction.metabolites} != row["species"]
            or reaction.annotation != row["annotation"]
            or gpr not in (row["before_gpr"], row["after_gpr"])):
        raise ValueError("R673 prolyl split precondition differs; no edits applied")
    notes, current = row["notes"], _notes(reaction.notes)
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError("Conflicting R673 prolyl split notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != row["before_notes"]:
        raise ValueError("R673 notes outside this curation differ; no edits applied")
    applied = gpr == row["after_gpr"]
    if applied and any(current.get(k) != v for k, v in notes.items()):
        raise ValueError("R673 already carries this GPR without this curation; no edits applied")
    if not applied and any(k in current for k in notes):
        raise ValueError("R673 carries this curation's notes without its GPR; no edits applied")
    for gene_id, gene_row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        expected = gene_row["reactions_after"] if applied else gene_row["reactions_before"]
        if gene.annotation != gene_row["annotation"] or sorted(r.id for r in gene.reactions) != expected:
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    record = {"item": ["R673"], "assignment_id": spec["assignment_id"], "unlinked_gene": UNLINKED,
              "evidence_status": spec["status"]}
    if applied:
        return {**record, "status": "already_correct", "gpr": AFTER}
    reaction.gene_reaction_rule = row["after_gpr"]  # history-tracked by cobra
    context = get_context(reaction)
    if context:
        # Notes are not history-tracked by cobra; restore them on context exit like the GPR.
        context(partial(setattr, reaction, "notes", reaction.notes))
    reaction.notes = {**reaction.notes, **notes}
    if (reaction.gene_reaction_rule != AFTER
            or any(sorted(r.id for r in model.genes.get_by_id(g).reactions) != spec["genes"][g]["reactions_after"]
                   for g in GENES)):
        raise RuntimeError("R673 prolyl split did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": BEFORE, "gpr": AFTER}
