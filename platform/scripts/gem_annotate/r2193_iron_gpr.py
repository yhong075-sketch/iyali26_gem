"""Opt-in GPR correction for R2193 (iron(2+) import): the two FET3/FTR1 pairs replace a bHLH regulator.

The decision and its evidence are data in ``r2193_iron_gpr.json``. The previous gene, YALI1F07747g, is
kept in the model with no reaction (as YALI1D17462g after the R153 revision). The four new genes get
their NCBI Gene, RefSeq, UniProt and KEGG cross-references. Every check runs before any edit.
"""

import copy
import json
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r2193_iron_gpr.json")
REACTION = "R2193"
OLD_GENE = "YALI1F07747g"
NEW_GENES = {"YALI1D09277g", "YALI1D09308g", "YALI1D08684g", "YALI1D08564g"}
AFTER_GPR = "(YALI1D09277g and YALI1D09308g) or (YALI1D08684g and YALI1D08564g)"


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def apply_r2193_iron_gpr(model, spec=None):
    """Replace R2193's GPR with the reviewed FET3/FTR1 rule; fail closed on any unexpected state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "R2193-IRON-GPR-20261022"
            or spec.get("status") != "user_authorized_gpr_correction"
            or spec.get("reaction_id") != REACTION
            or spec.get("before_gpr") != OLD_GENE
            or spec.get("after_gpr") != AFTER_GPR
            or spec.get("old_gene", {}).get("id") != OLD_GENE
            or set(spec.get("new_genes", {})) != NEW_GENES
            or not isinstance(spec.get("before_notes"), dict)):
        raise ValueError("Unexpected R2193 iron GPR scope/status")
    if REACTION not in model.reactions:
        raise ValueError(f"{REACTION} missing; no edits applied")
    reaction = model.reactions.get_by_id(REACTION)
    gpr = reaction.gene_reaction_rule
    if (reaction.name != spec["reaction_name"]
            or list(reaction.bounds) != spec["bounds"]
            or {m.id: c for m, c in reaction.metabolites.items()} != spec["stoichiometry"]
            or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                for m in reaction.metabolites} != spec["species"]
            or reaction.annotation != spec["annotation"]
            or gpr not in (spec["before_gpr"], spec["after_gpr"])):
        raise ValueError(f"{REACTION} iron GPR precondition differs; no edits applied")
    notes, current = spec["notes"], _notes(reaction.notes)
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError(f"Conflicting {REACTION} iron GPR notes; no edits applied")
    if {k: v for k, v in current.items() if k not in notes} != spec["before_notes"]:
        raise ValueError(f"{REACTION} notes outside this curation differ; no edits applied")
    done = gpr == spec["after_gpr"]
    if done and any(current.get(k) != v for k, v in notes.items()):
        # The after-state without these notes came from somewhere else; its evidence would be unrecorded.
        raise ValueError(f"{REACTION} already carries this GPR without this curation; no edits applied")
    if not done and any(k in current for k in notes):
        raise ValueError(f"{REACTION} carries this curation's notes without its GPR; no edits applied")
    old = spec["old_gene"]
    if OLD_GENE not in model.genes:
        raise ValueError(f"{OLD_GENE} missing; no edits applied")
    gene = model.genes.get_by_id(OLD_GENE)
    expected_reactions = [] if done else old["reactions_before"]
    if gene.annotation != old["annotation"] or sorted(r.id for r in gene.reactions) != expected_reactions:
        raise ValueError(f"{OLD_GENE} identity or reaction context differs; no edits applied")
    for gene_id, row in spec["new_genes"].items():
        if gene_id in model.genes:
            new = model.genes.get_by_id(gene_id)
            # No verified gene name exists; the SBML writer stores the ID-derived default.
            if (not done or new.annotation != row["annotation"] or new.name not in ("", f"G_{gene_id}")
                    or {r.id for r in new.reactions} != {REACTION}):
                raise ValueError(f"{gene_id} already exists with other identity or reactions; no edits applied")
        elif done:
            raise ValueError(f"{gene_id} missing although {REACTION} carries the curated GPR; no edits applied")
    record = {"item": REACTION, "assignment_id": spec["assignment_id"], "unlinked_gene": OLD_GENE,
              "evidence_status": spec["status"]}
    if done:
        return {**record, "status": "already_correct", "gpr": gpr, "added_genes": []}
    added = sorted(NEW_GENES)
    # cobra adds genes that first appear in a GPR in set order, which depends on PYTHONHASHSEED.
    # Growing the rule one new gene at a time creates them in sorted order through cobra's own
    # (context-aware) path, so the written model does not depend on the seed.
    for count in range(1, len(added) + 1):
        reaction.gene_reaction_rule = " or ".join(added[:count])
    reaction.gene_reaction_rule = spec["after_gpr"]
    for gene_id, row in spec["new_genes"].items():
        model.genes.get_by_id(gene_id).annotation = copy.deepcopy(row["annotation"])
    context = get_context(reaction)
    if context:
        # Notes are not history-tracked by cobra; restore them on context exit like the GPR.
        context(partial(setattr, reaction, "notes", reaction.notes))
    reaction.notes = {**reaction.notes, **notes}
    if gene.reactions or any({r.id for r in model.genes.get_by_id(g).reactions} != {REACTION} for g in NEW_GENES):
        raise RuntimeError(f"{REACTION} iron GPR did not reach the curated state")
    return {**record, "status": "applied", "gpr_before": spec["before_gpr"], "gpr": reaction.gene_reaction_rule,
            "added_genes": added}
