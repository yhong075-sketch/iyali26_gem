"""Opt-in THF, FAD, riboflavin, NAD, NADP and CoA terms in biomass_C at the yeast-GEM 9.1.1 values (2026-10-27).

The Ponce-de-Leon 2013 census of a1237821 found these six cofactors synthesised but never drained
(coenzyme pseudo-gaps). The remedy, as for TDP and heme a, is a growth-dilution term in biomass_C.
The decision and its evidence are data in ``cofactor_biomass_terms.json``.

The step fails closed unless biomass_C is exactly the reviewed a1237821 composition (so it must run after
every other step that edits biomass_C) and each species has its reviewed identity and reaction context.
A rerun, also on a saved and reloaded model, is a no-op.
"""

import json
import math

from .config import MODEL

CURATION_PATH = MODEL.curation_file("cofactor_biomass_terms.json")
CURATION_ID = "COFACTOR-BIOMASS-TERMS-20261027"
BIOMASS = "biomass_C"
TERM_IDS = {"m479[C_cy]", "m528[C_cy]", "m926[C_cy]", "m122[C_cy]", "m41[C_cy]", "m69[C_cy]"}
# The one reviewed starting state of biomass_C: candidate a1237821 (2026-10-23).
BEFORE_STATE_ID = "a1237821_biomass_C"


def _notes(spec):
    terms = sorted(spec["terms"].items())
    return {"cofactor_biomass_terms_curation_id": spec["curation_id"],
            "cofactor_biomass_terms_status": spec["status"],
            "cofactor_biomass_terms_authorization": spec["authorization"],
            "cofactor_biomass_terms_terms": "; ".join(f"{mid} {row['coefficient_mmol_per_gDW']:g} mmol/gDW" for mid, row in terms),
            "cofactor_biomass_terms_evidence": " | ".join(row["source"] for _, row in terms),
            "cofactor_biomass_terms_evidence_tier": spec["evidence_tier"],
            "cofactor_biomass_terms_caveat": spec["organism_caveat"]}


def _same(value, expected):
    return math.isclose(value, expected, rel_tol=1e-12, abs_tol=0.0)


def apply_cofactor_biomass_terms(model, spec=None):
    """Add the six cofactor terms to biomass_C; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    before = spec.get("accepted_before_state") or {}
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != CURATION_ID
            or spec.get("status") != "user_authorized_biomass_terms"
            or spec.get("biomass_reaction") != BIOMASS
            or spec.get("requires_objective") != BIOMASS
            or set(spec.get("terms", {})) != TERM_IDS
            or not all(row.get("coefficient_mmol_per_gDW", 0) > 0 for row in spec["terms"].values())
            or before.get("id") != BEFORE_STATE_ID
            or not isinstance(before.get("biomass_stoichiometry"), dict)
            or TERM_IDS & set(before["biomass_stoichiometry"])
            or not isinstance(before.get("biomass_bounds"), list)
            or not isinstance(before.get("biomass_note_keys"), list)):
        raise ValueError("Unexpected cofactor biomass terms scope/status")
    if BIOMASS not in model.reactions:
        raise ValueError(f"{BIOMASS} missing; no edits applied")
    if {r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient} != {BIOMASS: 1.0}:
        raise ValueError(f"The objective is not {BIOMASS} alone; no edits applied")
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
    expected = before["biomass_stoichiometry"]
    if (list(biomass.bounds) != before["biomass_bounds"] or set(others) != set(expected)
            or not all(_same(others[mid], c) for mid, c in expected.items())):
        raise ValueError(f"{BIOMASS} composition differs from the reviewed baseline; no edits applied")
    notes = _notes(spec)
    if any(key in biomass.notes and biomass.notes[key] != value for key, value in notes.items()):
        raise ValueError(f"Conflicting {BIOMASS} cofactor notes; no edits applied")
    if sorted(key for key in biomass.notes if key not in notes) != sorted(before["biomass_note_keys"]):
        raise ValueError(f"{BIOMASS} notes outside this curation differ; no edits applied")
    present = {mid for mid in TERM_IDS if mid in stoichiometry}
    for mid in sorted(present):
        if not _same(stoichiometry[mid], -terms[mid]["coefficient_mmol_per_gDW"]):
            raise ValueError(f"Conflicting {mid} biomass coefficient; no edits applied")
    if present and present != TERM_IDS:
        raise ValueError("Only some of the cofactor biomass terms are present; no edits applied")
    if present and not all(key in biomass.notes for key in notes):
        # Terms without these notes came from somewhere else; their evidence would be unrecorded.
        raise ValueError(f"{BIOMASS} already carries the cofactor terms without this curation; no edits applied")
    if not present and any(key in biomass.notes for key in notes):
        raise ValueError(f"{BIOMASS} carries this curation's notes without the terms; no edits applied")
    if not present:
        biomass.add_metabolites({model.metabolites.get_by_id(mid): -row["coefficient_mmol_per_gDW"]
                                 for mid, row in sorted(terms.items())})
        biomass.notes = {**biomass.notes, **notes}
    return {"item": "cofactor biomass terms", "curation_id": spec["curation_id"],
            "status": "already_correct" if present else "applied", "before_state": before["id"],
            "changed_terms": [] if present else sorted(TERM_IDS), "evidence_status": spec["status"]}
