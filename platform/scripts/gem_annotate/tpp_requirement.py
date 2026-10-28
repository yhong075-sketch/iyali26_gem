"""Opt-in TPP requirement: TPP (m975) diluted into biomass_C and the gene-less TMP kinase R1342 closed."""

import json
from xml.sax.saxutils import unescape

from .config import MODEL

CURATION_PATH = MODEL.curation_file("tpp_requirement.json")
TARGET_IDS = {"biomass_C", "R1342"}
# Every reaction on cytosolic TPP besides the targets; R753 and R933 lead to dead ends m976 and m1116.
REQUIRED_IDS = {"R752", "R753", "R933"}
TPP = "m975[C_cy]"
SPECIES_IDS = {TPP, "m976[C_cy]", "m1116[C_mi]"}
XML_ENTITIES = {"&apos;": "'", "&quot;": '"'}


def _unescaped(notes):
    # biomass_C notes hold -> and " in a build, but read back from SBML as -&gt; and &quot;.
    return {key: unescape(value, XML_ENTITIES) for key, value in notes.items()}


def _notes(spec, row):
    return {"tpp_requirement_curation_id": spec["curation_id"],
            "tpp_requirement_status": spec["status"],
            "tpp_requirement_authorization": spec["authorization"],
            "tpp_requirement_evidence": row["evidence"],
            "tpp_requirement_evidence_tier": row["evidence_tier"],
            "tpp_requirement_rationale": row["rationale"],
            "tpp_requirement_previous_state": row["previous_state"]}


def _check_spec(spec):
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "TPP-REQUIREMENT-20261022"
            or spec.get("status") != "user_authorized_biomass_term_and_closure"
            or set(spec.get("targets", {})) != TARGET_IDS
            or set(spec.get("requires_reactions", {})) != REQUIRED_IDS
            or set(spec.get("species", {})) != SPECIES_IDS):
        raise ValueError("Unexpected TPP requirement scope/status")
    # The reviewed change is exactly one negative TPP term and R1342 closed; nothing else may ride along.
    biomass, closure = spec["targets"]["biomass_C"], spec["targets"]["R1342"]
    coefficient = biomass["after_stoichiometry"].get(TPP, 0.0)
    if (TPP in biomass["before_stoichiometry"] or not coefficient < 0
            or biomass["after_stoichiometry"] != {**biomass["before_stoichiometry"], TPP: coefficient}
            or biomass["bounds_after"] != biomass["bounds_before"]
            or closure["after_stoichiometry"] != closure["before_stoichiometry"]
            or closure["bounds_after"] != [0.0, 0.0]):
        raise ValueError("TPP requirement spec values differ from the reviewed change")


def apply_tpp_requirement(model, spec=None):
    """Add the reviewed TPP biomass term and close R1342; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    if {r.id for r in model.reactions if r.objective_coefficient} != {"biomass_C"}:
        raise ValueError("biomass_C is not the sole objective; no edits applied")
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
        # Another TPP producer would be a route the knockout analysis did not review.
        if sorted(r.id for r in metabolite.reactions if r.id != "biomass_C") != row["reactions_before"]:
            raise ValueError(f"{mid} reaction context differs; no edits applied")
    for rid, row in spec["requires_reactions"].items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        if (reaction.gene_reaction_rule != row["gpr"] or list(reaction.bounds) != row["bounds"]
                or {m.id: c for m, c in reaction.metabolites.items()} != row["stoichiometry"]):
            raise ValueError(f"{rid} differs from the reviewed TPP context; no edits applied")
    rows = spec["targets"]
    done = set()
    for rid, row in rows.items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
        before = (stoichiometry == row["before_stoichiometry"] and list(reaction.bounds) == row["bounds_before"])
        after = (stoichiometry == row["after_stoichiometry"] and list(reaction.bounds) == row["bounds_after"])
        if (reaction.name != row["name"] or reaction.gene_reaction_rule != row["gpr"]
                or reaction.annotation != row["annotation"] or not (before or after)):
            raise ValueError(f"{rid} precondition differs; no edits applied")
        notes, current = _notes(spec, row), _unescaped(reaction.notes)
        if any(key in current and current[key] != value for key, value in notes.items()):
            raise ValueError(f"Conflicting {rid} TPP requirement notes; no edits applied")
        if {k: v for k, v in current.items() if k not in notes} != row["notes_before"]:
            raise ValueError(f"{rid} notes outside this curation differ; no edits applied")
        if after:
            # The after-state without these notes came from somewhere else; its evidence would be unrecorded.
            if current != {**row["notes_before"], **notes}:
                raise ValueError(f"{rid} already carries the TPP change without this curation; no edits applied")
            done.add(rid)
        elif any(key in current for key in notes):
            raise ValueError(f"{rid} carries this curation's notes without the change; no edits applied")
    if done and done != TARGET_IDS:
        raise ValueError("Only one TPP requirement change is present; no edits applied")
    if not done:
        tpp = model.metabolites.get_by_id(TPP)
        model.reactions.R1342.bounds = tuple(rows["R1342"]["bounds_after"])
        model.reactions.biomass_C.add_metabolites({tpp: rows["biomass_C"]["after_stoichiometry"][TPP]})
        for rid, row in rows.items():
            reaction = model.reactions.get_by_id(rid)
            reaction.notes = {**reaction.notes, **_notes(spec, row)}
    return {"item": "TPP requirement", "curation_id": spec["curation_id"],
            "status": "already_correct" if done else "applied",
            "changed_reactions": [] if done else sorted(rows), "evidence_status": spec["status"]}
