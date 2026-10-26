"""Opt-in heme a and iron(2+) terms in biomass_C, at the Yeast8 consensus values (2026-10-22)."""

import json
import math

from .config import MODEL

CURATION_PATH = MODEL.curation_file("heme_iron_biomass.json")
BIOMASS = "biomass_C"
TERM_IDS = {"m1635[C_cy]", "m948[C_cy]"}
# Iron import (R1189, R2193, R863), heme b -> heme o -> heme a (R302, R384, R383) and heme a export (R1709).
SUPPLY_IDS = {"R1189", "R2193", "R863", "R302", "R384", "R383", "R1709"}
# Reviewed starting states of biomass_C: the base, and the base after the TPP requirement (2026-10-22).
BEFORE_STATE_IDS = ["base", "after_TPP-REQUIREMENT-20261022"]


def _notes(spec):
    terms = sorted(spec["terms"].items())
    return {"heme_iron_biomass_curation_id": spec["curation_id"],
            "heme_iron_biomass_status": spec["status"],
            "heme_iron_biomass_authorization": spec["authorization"],
            "heme_iron_biomass_terms": "; ".join(f"{mid} {row['coefficient_mmol_per_gDW']:g} mmol/gDW"
                                                 for mid, row in terms),
            "heme_iron_biomass_evidence": " | ".join(row["source"] for _, row in terms),
            "heme_iron_biomass_evidence_tier": spec["evidence_tier"],
            "heme_iron_biomass_caveat": spec["organism_caveat"]}


def _same(value, expected):
    return math.isclose(value, expected, rel_tol=1e-12, abs_tol=0.0)


def apply_heme_iron_biomass(model, spec=None):
    """Add heme a and iron(2+) to biomass_C; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "HEME-IRON-BIOMASS-20261022"
            or spec.get("status") != "user_authorized_biomass_term"
            or spec.get("biomass_reaction") != BIOMASS
            or set(spec.get("terms", {})) != TERM_IDS
            or set(spec.get("requires_reactions", {})) != SUPPLY_IDS
            or spec.get("requires_objective") != BIOMASS
            or [s.get("id") for s in spec.get("accepted_before_states", [])] != BEFORE_STATE_IDS
            or spec["accepted_before_states"][0]["added_stoichiometry"]
            or spec["accepted_before_states"][0]["added_note_keys"]
            or any(TERM_IDS & set(s["added_stoichiometry"]) for s in spec["accepted_before_states"])):
        raise ValueError("Unexpected heme/iron biomass scope/status")
    if BIOMASS not in model.reactions:
        raise ValueError(f"{BIOMASS} missing; no edits applied")
    if {r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient} != {BIOMASS: 1.0}:
        raise ValueError(f"The objective is not {BIOMASS} alone; no edits applied")
    for rid, row in spec["requires_reactions"].items():
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
        supply = model.reactions.get_by_id(rid)
        if (list(supply.bounds) != row["bounds"]
                or {m.id: c for m, c in supply.metabolites.items()} != row["stoichiometry"]):
            raise ValueError(f"{rid} differs from the reviewed supply route; no edits applied")
    terms = spec["terms"]
    for mid, row in terms.items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
        if sorted(r.id for r in metabolite.reactions if r.id != BIOMASS) != row["reactions_before"]:
            raise ValueError(f"{mid} reaction context differs; no edits applied")
    biomass = model.reactions.get_by_id(BIOMASS)
    stoichiometry = {m.id: c for m, c in biomass.metabolites.items()}
    others = {mid: c for mid, c in stoichiometry.items() if mid not in TERM_IDS}
    # Each accepted starting state is the reviewed base plus another curation's reviewed additions.
    matched = None
    for state in spec["accepted_before_states"]:
        before = {**spec["biomass_stoichiometry_before"], **state["added_stoichiometry"]}
        if (list(biomass.bounds) == spec["biomass_bounds_before"] and set(others) == set(before)
                and all(_same(others[mid], c) for mid, c in before.items())):
            matched = state
            break
    if matched is None:
        raise ValueError(f"{BIOMASS} composition differs from the reviewed baseline; no edits applied")
    notes = _notes(spec)
    if any(key in biomass.notes and biomass.notes[key] != value for key, value in notes.items()):
        raise ValueError(f"Conflicting {BIOMASS} heme/iron notes; no edits applied")
    if (sorted(key for key in biomass.notes if key not in notes)
            != sorted(spec["biomass_note_keys_before"] + matched["added_note_keys"])):
        raise ValueError(f"{BIOMASS} notes outside this curation differ; no edits applied")
    present = {mid for mid in TERM_IDS if mid in stoichiometry}
    for mid in sorted(present):
        if not _same(stoichiometry[mid], -terms[mid]["coefficient_mmol_per_gDW"]):
            raise ValueError(f"Conflicting {mid} biomass coefficient; no edits applied")
    if present and present != TERM_IDS:
        raise ValueError("Only one heme/iron biomass term is present; no edits applied")
    if present and not all(key in biomass.notes for key in notes):
        # Terms without these notes came from somewhere else; their evidence would be unrecorded.
        raise ValueError(f"{BIOMASS} already carries the heme/iron terms without this curation; no edits applied")
    if not present and any(key in biomass.notes for key in notes):
        raise ValueError(f"{BIOMASS} carries this curation's notes without the terms; no edits applied")
    if not present:
        biomass.add_metabolites({model.metabolites.get_by_id(mid): -row["coefficient_mmol_per_gDW"]
                                 for mid, row in sorted(terms.items())})
        biomass.notes = {**biomass.notes, **notes}
    return {"item": "heme a + iron biomass terms", "curation_id": spec["curation_id"],
            "status": "already_correct" if present else "applied", "before_state": matched["id"],
            "changed_terms": [] if present else sorted(TERM_IDS), "evidence_status": spec["status"]}
