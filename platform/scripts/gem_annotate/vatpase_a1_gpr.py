"""Opt-in V-ATPase rule for R794 (Golgi) and R795 (vacuole): one core AND either a-subunit, plus V0 e.

Both reactions read `(12-gene core with a1) or (a2 and c'')`, as inherited from iYli21; the second branch has no V1
head and no c/c' ring, so no single knockout closes either reaction. The rule becomes the 12 V1/V0 core subunits (c''
moved in) and V0 e (`YALI1F25663g`, a new model gene) AND (a1 or a2). The decision and its evidence are data in
``vatpase_a1_gpr.json``. Every check runs before the edit.
"""

import copy
import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("vatpase_a1_gpr.json")
REACTIONS = ("R794", "R795")
CORE = ("YALI1A11258g", "YALI1E14125g", "YALI1D00581g", "YALI1E10492g", "YALI1A09766g", "YALI1F20965g",
        "YALI0E16192g", "YALI1F31854g", "YALI1F21690g", "YALI1E37063g", "YALI1E32332g", "YALI1F13017g")
NEW_GENE = "YALI1F25663g"
A_SUBUNITS = ("YALI1E12482g", "YALI1F38820g")
BEFORE = ("(YALI1A11258g and YALI1E14125g and YALI1D00581g and YALI1E10492g and YALI1A09766g and YALI1F20965g and "
          "YALI0E16192g and YALI1F31854g and YALI1F21690g and YALI1E37063g and YALI1E32332g and YALI1E12482g) or "
          "(YALI1F38820g and YALI1F13017g)")
AFTER = " and ".join(CORE + (NEW_GENE,)) + f" and ({' or '.join(A_SUBUNITS)})"


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _check_spec(spec):
    rows, genes, new = spec.get("reactions", {}), spec.get("genes", {}), spec.get("new_gene", {})
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "VATPASE-A1-GPR-20261103"
            or spec.get("status") != "user_authorized_gpr_correction"
            or set(rows) != set(REACTIONS)
            or any(row.get("before_gpr") != BEFORE or row.get("after_gpr") != AFTER
                   or not isinstance(row.get("before_notes"), dict)
                   or not row.get("notes") or row["notes"] != spec.get("notes") for row in rows.values())
            or set(genes) != set(CORE) | set(A_SUBUNITS)
            or any(not set(REACTIONS) <= set(genes[g].get("reactions_before", ()))
                   or genes[g].get("reactions_after") != genes[g]["reactions_before"] for g in genes)
            or new.get("id") != NEW_GENE or new.get("reactions_after") != list(REACTIONS)
            or not new.get("annotation")):
        raise ValueError("Unexpected V-ATPase A1 GPR scope/status")


def _check_reaction(model, row, rid):
    """Pin one reaction's before- or after-state; return True when it already carries the curated rule."""
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
        raise ValueError(f"{rid} V-ATPase A1 GPR precondition differs; no edits applied")
    notes, current = row["notes"], _notes(reaction.notes)
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError(f"Conflicting {rid} V-ATPase A1 GPR notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != row["before_notes"]:
        raise ValueError(f"{rid} notes outside this curation differ; no edits applied")
    applied = gpr == row["after_gpr"]
    if applied and any(current.get(k) != v for k, v in notes.items()):
        raise ValueError(f"{rid} already carries this GPR without this curation; no edits applied")
    if not applied and any(k in current for k in notes):
        raise ValueError(f"{rid} carries this curation's notes without its GPR; no edits applied")
    return applied


def apply_vatpase_a1_gpr(model, spec=None):
    """Set R794 and R795 to the core AND (a1 or a2) rule with V0 e; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    states = {rid: _check_reaction(model, spec["reactions"][rid], rid) for rid in REACTIONS}
    if len(set(states.values())) != 1:
        raise ValueError("Only one of R794 and R795 carries the curated rule; no edits applied")
    applied = states["R794"]
    for gene_id, row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        if gene.annotation != row["annotation"] or sorted(r.id for r in gene.reactions) != row["reactions_before"]:
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    new = spec["new_gene"]
    if NEW_GENE in model.genes:
        gene = model.genes.get_by_id(NEW_GENE)
        # No verified gene name exists; the SBML writer stores the ID-derived default.
        if (not applied or gene.annotation != new["annotation"] or gene.name not in ("", f"G_{NEW_GENE}")
                or sorted(r.id for r in gene.reactions) != new["reactions_after"]):
            raise ValueError(f"{NEW_GENE} already exists with other identity or reactions; no edits applied")
    elif applied:
        raise ValueError(f"{NEW_GENE} missing although R794/R795 carry the curated GPR; no edits applied")
    record = {"item": list(REACTIONS), "assignment_id": spec["assignment_id"], "new_gene": NEW_GENE,
              "evidence_status": spec["status"]}
    if applied:
        return {**record, "status": "already_correct", "gpr": AFTER, "added_genes": []}
    # One new gene only, so cobra's set-ordered gene creation cannot depend on PYTHONHASHSEED.
    for rid in REACTIONS:
        reaction = model.reactions.get_by_id(rid)
        reaction.gene_reaction_rule = AFTER  # history-tracked by cobra
        context = get_context(reaction)
        if context:
            # Notes are not history-tracked by cobra; restore them on context exit like the GPR.
            context(partial(setattr, reaction, "notes", reaction.notes))
        reaction.notes = {**reaction.notes, **spec["notes"]}
    model.genes.get_by_id(NEW_GENE).annotation = copy.deepcopy(new["annotation"])
    if (any(model.reactions.get_by_id(rid).gene_reaction_rule != AFTER for rid in REACTIONS)
            or sorted(r.id for r in model.genes.get_by_id(NEW_GENE).reactions) != new["reactions_after"]
            or any(sorted(r.id for r in model.genes.get_by_id(g).reactions) != row["reactions_after"]
                   for g, row in spec["genes"].items())):
        raise RuntimeError("V-ATPase A1 GPR did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": BEFORE, "gpr": AFTER, "added_genes": [NEW_GENE]}
