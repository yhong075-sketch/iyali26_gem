"""User-authorized 2026-10-19 removal of the A1 peptidase gene from R153 (argininosuccinate synthase).

The decision and its evidence are data in ``r153_gpr_revision.json``. The step runs after the
R153/R2176 merge. Before any edit it checks the exact merged R153 definition, the recorded notes, and
the gene's identity and R153-only membership, and fails closed on anything else. The gene is kept with
no reactions, so a rerun is a no-op. Builds start from the source model and merge before this step, so
``merge_r153_r2176`` never sees the revised R153.
"""

import json

from .config import MODEL

SPEC_PATH = MODEL.curation_file("r153_gpr_revision.json")


def _definition(reaction):
    return ({m.id: float(c) for m, c in reaction.metabolites.items()},
            [float(b) for b in reaction.bounds],
            {m.id: {"formula": m.formula, "charge": m.charge, "compartment": m.compartment}
             for m in reaction.metabolites})


def apply_r153_gpr_revision(model, spec=None):
    spec = json.loads(SPEC_PATH.read_text()) if spec is None else spec
    before, gene_row = spec.get("before", {}), spec.get("gene", {})
    if (spec.get("schema_version") != 1 or spec.get("change_id") != "R153-GPR-UNASSIGNED-20261019"
            or spec.get("reaction") != "R153" or gene_row.get("id") != "YALI1D17462g"
            or before.get("gpr") != "YALI1D17462g" or spec.get("after", {}).get("gpr") != ""):
        raise ValueError("Unexpected R153 GPR revision scope")
    if "R153" not in model.reactions or "R2176" in model.reactions:
        raise ValueError("R153 GPR revision expects the merged R153 without R2176; no edits applied")
    reaction = model.reactions.get_by_id("R153")
    if _definition(reaction) != (before["stoichiometry"], before["bounds"], before["species"]):
        raise ValueError("R153: definition differs from the curated record; no edits applied")
    if gene_row["id"] not in model.genes:
        raise ValueError("R153 GPR revision gene missing: YALI1D17462g")
    gene = model.genes.get_by_id(gene_row["id"])
    if any(gene.annotation.get(k) != v for k, v in gene_row["annotation"].items()):
        raise ValueError("YALI1D17462g: identity differs from the curated record; no edits applied")
    if any(reaction.notes.get(k) != v for k, v in before["required_notes"].items()):
        raise ValueError("R153: merge notes differ from the curated record; no edits applied")
    replaced, added = spec["notes_replaced"], spec["notes_added"]
    after_notes = {**{k: row["after"] for k, row in replaced.items()}, **added}
    record = {"change_id": spec["change_id"], "status_label": spec["status"], "reaction": "R153",
              "unlinked_gene": gene.id}
    if (reaction.gene_reaction_rule == spec["after"]["gpr"] and not gene.reactions
            and all(reaction.notes.get(k) == v for k, v in after_notes.items())):
        return {**record, "status": "already_correct", "gpr": reaction.gene_reaction_rule}
    if reaction.gene_reaction_rule != before["gpr"] or {r.id for r in gene.reactions} != {"R153"}:
        raise ValueError("R153: GPR or gene membership differs from the curated before-state; no edits applied")
    if any(reaction.notes.get(k) != row["before"] for k, row in replaced.items()):
        raise ValueError("R153: notes differ from the curated before-state; no edits applied")
    if any(k in reaction.notes for k in added):
        raise ValueError("R153: revision notes already present; no edits applied")
    reaction.gene_reaction_rule = spec["after"]["gpr"]
    reaction.notes = {**reaction.notes, **after_notes}
    if gene.id not in model.genes or gene.reactions:
        raise RuntimeError("R153 GPR revision did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": before["gpr"], "gpr": reaction.gene_reaction_rule}
