"""Opt-in co-products for the GPI-anchor chain: DAG on R373 (EtNP transfer from PE), Dol-P on R375 (third mannose)."""

import json

from .config import MODEL

CURATION_PATH = MODEL.curation_file("gpi_coproducts.json")
REACTION_IDS = {"R373", "R375"}
LINK_IDS = {"R_GPIGNT", "R_PEt_em_er", "R_DAGt_er_em", "R_ACt_er_cy", "R_DM_GPI"}


def _notes(spec, row):
    return {"gpi_coproducts_curation_id": spec["curation_id"],
            "gpi_coproducts_status": spec["status"],
            "gpi_coproducts_authorization": spec["authorization"],
            "gpi_coproducts_evidence": row["evidence"],
            "gpi_coproducts_evidence_tier": row["evidence_tier"],
            "gpi_coproducts_rationale": row["rationale"],
            "gpi_coproducts_previous_equation": row["previous_equation"]}


def apply_gpi_coproducts(model, spec=None):
    """Add the reviewed co-products to R373 and R375; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "GPI-COPRODUCTS-20261021"
            or spec.get("status") != "user_authorized_stoichiometry_correction"
            or set(spec.get("reactions", {})) != REACTION_IDS
            or set(spec.get("requires_reactions", {})) != LINK_IDS):
        raise ValueError("Unexpected GPI co-products scope/status")
    missing = sorted(rid for rid in LINK_IDS if rid not in model.reactions)
    if missing:
        raise ValueError(f"GPI co-products require the GPI chain links ({missing} missing); no edits applied")
    for rid, row in spec["requires_reactions"].items():
        link = model.reactions.get_by_id(rid)
        if (list(link.bounds) != row["bounds"]
                or {m.id: c for m, c in link.metabolites.items()} != row["stoichiometry"]):
            raise ValueError(f"{rid} differs from the GPI chain links; no edits applied")
    rows = spec["reactions"]
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
        if ("reactions_before" in row
                and sorted(r.id for r in metabolite.reactions if r.id not in REACTION_IDS) != row["reactions_before"]):
            raise ValueError(f"{mid} reaction context differs; no edits applied")
    done = set()
    for rid, row in rows.items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
        if (reaction.name != row["name"] or list(reaction.bounds) != row["bounds"]
                or reaction.gene_reaction_rule != row["gpr"] or reaction.annotation != row["annotation"]
                or stoichiometry not in (row["before_stoichiometry"], row["after_stoichiometry"])):
            raise ValueError(f"{rid} precondition differs; no edits applied")
        notes = _notes(spec, row)
        if any(key in reaction.notes and reaction.notes[key] != value for key, value in notes.items()):
            raise ValueError(f"Conflicting {rid} co-product notes; no edits applied")
        if {k: v for k, v in reaction.notes.items() if k not in notes} != row["notes_before"]:
            raise ValueError(f"{rid} notes outside this curation differ; no edits applied")
        if stoichiometry == row["after_stoichiometry"]:
            # The after-state without these notes came from somewhere else; its evidence would be unrecorded.
            if reaction.notes != {**row["notes_before"], **notes}:
                raise ValueError(f"{rid} already carries the co-product without this curation; no edits applied")
            done.add(rid)
        elif any(key in reaction.notes for key in notes):
            raise ValueError(f"{rid} carries this curation's notes without the co-product; no edits applied")
    if done and done != REACTION_IDS:
        raise ValueError("Only one GPI co-product is present; no edits applied")
    if not done:
        for rid, row in rows.items():
            reaction = model.reactions.get_by_id(rid)
            before, after = row["before_stoichiometry"], row["after_stoichiometry"]
            delta = {mid: after.get(mid, 0.0) - before.get(mid, 0.0) for mid in sorted(set(before) | set(after))}
            reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in delta.items() if c})
            reaction.notes = {**reaction.notes, **_notes(spec, row)}
    return {"item": "GPI co-products", "curation_id": spec["curation_id"],
            "status": "already_correct" if done else "applied",
            "changed_reactions": [] if done else sorted(rows), "evidence_status": spec["status"]}
