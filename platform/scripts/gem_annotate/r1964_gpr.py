"""Opt-in provisional GPR and D-xylulose cross-references for the orphan R1964 (arabinitol dehydrogenase)."""

import copy
import json

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r1964_gpr_assignment.json")


def apply_r1964_gpr_assignment(model, spec=None):
    """Assign the reviewed GPR, cross-references and notes to R1964; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("assignment_id") != "R1964-ARADH-GPR-20261019"
            or spec.get("status") != "user_authorized_provisional_assignment"
            or spec.get("reaction_id") != "R1964"
            or spec.get("before_gpr") != ""
            or spec.get("after_gpr") != "YALI1F03452g or YALI1E06616g"
            or set(spec.get("existing_genes", {})) != {"YALI1F03452g"}
            or set(spec.get("new_genes", {})) != {"YALI1E06616g"}
            or not isinstance(spec.get("before_notes"), dict)):
        raise ValueError("Unexpected R1964 GPR assignment scope/status")
    if "R1964" not in model.reactions:
        raise ValueError("R1964 missing; no edits applied")
    reaction = model.reactions.get_by_id("R1964")
    state = (reaction.gene_reaction_rule, reaction.annotation)
    after = (spec["after_gpr"], spec["after_annotation"])
    if (reaction.name != spec["reaction_name"]
            or list(reaction.bounds) != spec["bounds"]
            or {m.id: c for m, c in reaction.metabolites.items()} != spec["stoichiometry"]
            or {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                for m in reaction.metabolites} != spec["species"]
            or state not in ((spec["before_gpr"], spec["before_annotation"]), after)):
        raise ValueError("R1964 GPR assignment precondition differs; no edits applied")
    if any(k in reaction.notes and reaction.notes[k] != v for k, v in spec["notes"].items()):
        raise ValueError("Conflicting R1964 GPR assignment notes; no edits applied")
    if {k: v for k, v in reaction.notes.items() if k not in spec["notes"]} != spec["before_notes"]:
        raise ValueError("R1964 notes outside this curation differ; no edits applied")
    # The after-state without these notes came from somewhere else; its evidence would be unrecorded.
    if state == after and any(reaction.notes.get(k) != v for k, v in spec["notes"].items()):
        raise ValueError("R1964 already carries this GPR without this curation; no edits applied")
    for gene_id, row in spec["existing_genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        if (gene.annotation != row["annotation"]
                or sorted(r.id for r in gene.reactions - {reaction}) != row["other_reactions"]):
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    for gene_id, row in spec["new_genes"].items():
        if gene_id in model.genes:
            gene = model.genes.get_by_id(gene_id)
            # No verified gene name exists; the SBML writer stores the ID-derived default.
            if (gene.annotation != row["annotation"] or gene.name not in ("", f"G_{gene_id}")
                    or gene.reactions - {reaction}):
                raise ValueError(f"{gene_id} already exists with other identity or reactions; no edits applied")
    notes = {**reaction.notes, **spec["notes"]}
    added = sorted(g for g in spec["new_genes"] if g not in model.genes)
    changed = state != after or notes != reaction.notes or bool(added)
    reaction.gene_reaction_rule = spec["after_gpr"]
    for gene_id, row in spec["new_genes"].items():
        model.genes.get_by_id(gene_id).annotation = copy.deepcopy(row["annotation"])
    reaction.annotation = copy.deepcopy(spec["after_annotation"])
    reaction.notes = notes
    return {"item": reaction.id, "assignment_id": spec["assignment_id"],
            "status": "applied" if changed else "already_correct",
            "gpr": reaction.gene_reaction_rule, "added_genes": added, "evidence_status": spec["status"]}
