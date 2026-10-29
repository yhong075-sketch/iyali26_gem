"""Opt-in R171 ATP synthase GPR: the single-gene OR branches and subunits e and j leave the rule.

R171 reads `(14-subunit AND) or YALI1B28025g or YALI1B15733g`, as inherited from iYli21; the OR branches are subunits
g and k, with no evidence that either synthesises ATP alone. The rule becomes the 12-gene AND of ATP6/8/9 and the nine nuclear core
subunits; g, k, e (`YALI1E38131g`) and j (`YALI0D17490g`) are kept in the model with no reaction. The decision and its
evidence are data in ``r171_atp_synthase_gpr.json``. Every check runs before the edit.
"""

import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r171_atp_synthase_gpr.json")
CORE = ("YALI1M00260g", "YALI1M00267g", "YALI1M00318g", "YALI1F04433g", "YALI1F04110g", "YALI1D27826g",
        "YALI1D15734g", "YALI1D14639g", "YALI1B08885g", "YALI1F07384g", "YALI1B05364g", "YALI1F26967g")
UNLINKED = ("YALI1B28025g", "YALI1B15733g", "YALI1E38131g", "YALI0D17490g")
BEFORE = ("(YALI1M00260g and YALI1M00267g and YALI1M00318g and YALI1F04433g and YALI1F04110g and YALI1D27826g and "
          "YALI1D15734g and YALI1D14639g and YALI1B08885g and YALI1F07384g and YALI1B05364g and YALI0D17490g and "
          "YALI1F26967g and YALI1E38131g) or YALI1B28025g or YALI1B15733g")
AFTER = " and ".join(CORE)


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _check_spec(spec):
    row, genes = spec.get("reaction", {}), spec.get("genes", {})
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "R171-ATP-SYNTHASE-GPR-20261029"
            or spec.get("status") != "user_authorized_gpr_correction"
            or row.get("id") != "R171"
            or row.get("before_gpr") != BEFORE or row.get("after_gpr") != AFTER
            or not isinstance(row.get("before_notes"), dict) or not row.get("notes")
            or set(genes) != set(CORE) | set(UNLINKED)
            or set(spec.get("unlinked_genes", {})) != set(UNLINKED)
            or any(genes[g].get("reactions_before") != ["R171"] for g in genes)
            or any(genes[g].get("reactions_after") != ([] if g in UNLINKED else ["R171"]) for g in genes)):
        raise ValueError("Unexpected R171 ATP synthase GPR scope/status")


def apply_r171_atp_synthase_gpr(model, spec=None):
    """Set R171 to the 12-gene core AND and unlink g, k, e and j; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    row = spec["reaction"]
    if "R171" not in model.reactions:
        raise ValueError("R171 missing; no edits applied")
    reaction = model.reactions.R171
    gpr = reaction.gene_reaction_rule
    if (reaction.name != row["reaction_name"]
            or list(reaction.bounds) != row["bounds"]
            or {m.id: c for m, c in reaction.metabolites.items()} != row["stoichiometry"]
            or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                for m in reaction.metabolites} != row["species"]
            or reaction.annotation != row["annotation"]
            or gpr not in (row["before_gpr"], row["after_gpr"])):
        raise ValueError("R171 ATP synthase GPR precondition differs; no edits applied")
    notes, current = row["notes"], _notes(reaction.notes)
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError("Conflicting R171 ATP synthase GPR notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != row["before_notes"]:
        raise ValueError("R171 notes outside this curation differ; no edits applied")
    applied = gpr == row["after_gpr"]
    if applied and any(current.get(k) != v for k, v in notes.items()):
        raise ValueError("R171 already carries this GPR without this curation; no edits applied")
    if not applied and any(k in current for k in notes):
        raise ValueError("R171 carries this curation's notes without its GPR; no edits applied")
    for gene_id, gene_row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        expected = gene_row["reactions_after"] if applied else gene_row["reactions_before"]
        if gene.annotation != gene_row["annotation"] or sorted(r.id for r in gene.reactions) != expected:
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    record = {"item": ["R171"], "assignment_id": spec["assignment_id"], "unlinked_genes": list(UNLINKED),
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
                   for g in spec["genes"])):
        raise RuntimeError("R171 ATP synthase GPR did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": BEFORE, "gpr": AFTER}
