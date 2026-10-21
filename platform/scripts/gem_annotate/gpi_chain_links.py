"""Opt-in links that connect the blocked GPI-anchor chain (R371-R376): step 1, three ER links and an open exit."""

import copy
import json

from cobra import Reaction

from .config import MODEL

CURATION_PATH = MODEL.curation_file("gpi_chain_links.json")
REACTION_IDS = {"R_GPIGNT", "R_PEt_em_er", "R_DAGt_er_em", "R_ACt_er_cy", "R_DM_GPI"}


def _notes(spec, row):
    return {"gpi_chain_links_curation_id": spec["curation_id"],
            "gpi_chain_links_status": spec["status"],
            "gpi_chain_links_authorization": spec["authorization"],
            "gpi_chain_links_evidence": row["evidence"],
            "gpi_chain_links_evidence_tier": row["evidence_tier"],
            "gpi_chain_links_rationale": row["rationale"]}


def _definition(reaction):
    return (reaction.name, list(reaction.bounds), reaction.gene_reaction_rule,
            {m.id: c for m, c in reaction.metabolites.items()}, reaction.annotation, reaction.notes)


def apply_gpi_chain_links(model, spec=None):
    """Add the five reviewed GPI chain links; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "GPI-CHAIN-LINKS-20261021"
            or spec.get("status") != "user_authorized_connectivity_fix"
            or set(spec.get("reactions", {})) != REACTION_IDS):
        raise ValueError("Unexpected GPI chain links scope/status")
    rows = spec["reactions"]
    present = {rid for rid in REACTION_IDS if rid in model.reactions}
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
        # A chain species that gained a producer or consumer elsewhere would make these links unreviewed.
        if ("reactions_before" in row
                and sorted(r.id for r in metabolite.reactions if r.id not in present) != row["reactions_before"]):
            raise ValueError(f"{mid} reaction context differs; no edits applied")
    for rid, row in spec["chain_reactions"].items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        if (list(reaction.bounds) != row["bounds"]
                or {m.id: c for m, c in reaction.metabolites.items()} != row["stoichiometry"]):
            raise ValueError(f"{rid} differs from the recorded GPI chain; no edits applied")
    groups = {}
    for name in sorted({row["group"] for row in rows.values()}):
        found = [g for g in model.groups if g.name == name]
        if len(found) != 1:
            raise ValueError(f"Group {name!r} is missing or not unique; no edits applied")
        groups[name] = found[0]
    for rid in sorted(present):
        reaction, row = model.reactions.get_by_id(rid), rows[rid]
        after = (row["name"], row["bounds"], row["gpr"], row["stoichiometry"], row["annotation"], _notes(spec, row))
        if (_definition(reaction) != after
                or [g.name for g in model.groups if reaction in g.members] != [row["group"]]):
            raise ValueError(f"{rid} already exists with another definition; no edits applied")
    if present and present != REACTION_IDS:
        raise ValueError("Only some GPI chain links are present; no edits applied")
    if not present:
        for rid, row in rows.items():
            reaction = Reaction(rid, name=row["name"], subsystem=row["group"],
                                lower_bound=row["bounds"][0], upper_bound=row["bounds"][1])
            model.add_reactions([reaction])
            reaction.add_metabolites({model.metabolites.get_by_id(m): c for m, c in row["stoichiometry"].items()})
            reaction.gene_reaction_rule = row["gpr"]
            reaction.annotation = copy.deepcopy(row["annotation"])
            reaction.notes = _notes(spec, row)
            groups[row["group"]].add_members([reaction])
    return {"item": "GPI chain links", "curation_id": spec["curation_id"],
            "status": "already_correct" if present else "applied",
            "added_reactions": [] if present else list(rows), "evidence_status": spec["status"]}
