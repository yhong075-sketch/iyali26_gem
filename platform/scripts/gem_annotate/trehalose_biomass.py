"""Opt-in removal of the trehalose term from biomass_C (2026-10-30).

biomass_C drains trehalose m294[C_cy] at 0.00686 mmol/gDW, a value inherited from the source model with no
Yarrowia lipolytica measurement behind it (the iYali4 carbohydrate share, scaled). Flores et al. 2011 measured
below 1 nmol/mg dry weight in the wild type and found that a tps1 mutant without detectable trehalose grows on
glucose, so the term is removed. The decision and its evidence are data in ``trehalose_biomass.json``; biomass_C
is not re-normalised and no other reaction changes.

The step fails closed unless biomass_C is the maximised objective and equals the reviewed e1e7c868 state (every
other coefficient, the bounds and the note keys; note values are not pinned), so it runs after every other step that
edits biomass_C, and unless the term and this curation's notes are either both absent or both present. A rerun, also
on a saved and reloaded model, is a no-op.
"""

import json
import math
from functools import partial
from html import unescape

from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("trehalose_biomass.json")
CURATION_ID = "TREHALOSE-BIOMASS-20261030"
BIOMASS, TERM, COEFFICIENT = "biomass_C", "m294[C_cy]", -0.00686


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _same(value, expected):
    return math.isclose(value, expected, rel_tol=1e-12, abs_tol=0.0)


def _check_spec(spec):
    term, before = spec.get("biomass_term") or {}, spec.get("accepted_before_state") or {}
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != CURATION_ID
            or spec.get("status") != "user_authorized_biomass_term_removal"
            or spec.get("requires_objective") != BIOMASS
            or term.get("biomass_reaction") != BIOMASS
            or term.get("metabolite_id") != TERM
            or term.get("coefficient_before") != COEFFICIENT
            or BIOMASS not in term.get("reactions_before", [])
            or term.get("reactions_after") != sorted(r for r in term["reactions_before"] if r != BIOMASS)
            or not isinstance(term.get("notes"), dict) or not term["notes"]
            or not isinstance(before.get("biomass_stoichiometry_without_term"), dict)
            or TERM in before["biomass_stoichiometry_without_term"]
            or not isinstance(before.get("biomass_bounds"), list)
            or not isinstance(before.get("biomass_note_keys"), list)):
        raise ValueError("Unexpected trehalose biomass scope/status")


def apply_trehalose_biomass(model, spec=None):
    """Remove trehalose m294[C_cy] from biomass_C; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    term, before = spec["biomass_term"], spec["accepted_before_state"]
    if BIOMASS not in model.reactions or TERM not in model.metabolites:
        raise ValueError(f"{BIOMASS} or {TERM} missing; no edits applied")
    if ({r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient} != {BIOMASS: 1.0}
            or model.objective.direction != "max"):
        raise ValueError(f"The objective is not {BIOMASS} alone, maximised; no edits applied")
    metabolite = model.metabolites.get_by_id(TERM)
    if (metabolite.name, metabolite.formula, metabolite.charge, metabolite.compartment) != (
            term["name"], term["formula"], term["charge"], term["compartment"]):
        raise ValueError(f"{TERM} identity differs; no edits applied")
    biomass = model.reactions.get_by_id(BIOMASS)
    stoichiometry = {m.id: c for m, c in biomass.metabolites.items()}
    present = TERM in stoichiometry
    coefficient = stoichiometry.pop(TERM, None)
    expected = before["biomass_stoichiometry_without_term"]
    if (list(biomass.bounds) != before["biomass_bounds"] or set(stoichiometry) != set(expected)
            or any(not _same(stoichiometry[k], v) for k, v in expected.items())):
        raise ValueError(f"{BIOMASS} composition differs from the reviewed baseline; no edits applied")
    if present and not _same(coefficient, term["coefficient_before"]):
        raise ValueError(f"Conflicting {TERM} biomass coefficient; no edits applied")
    reactions = sorted(r.id for r in metabolite.reactions)
    if reactions != (term["reactions_before"] if present else term["reactions_after"]):
        raise ValueError(f"{TERM} reaction context differs; no edits applied")
    notes, current = term["notes"], _notes(biomass.notes)
    if sorted(k for k in current if k not in notes) != sorted(before["biomass_note_keys"]):
        raise ValueError(f"{BIOMASS} notes outside this curation differ; no edits applied")
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError(f"Conflicting {BIOMASS} trehalose notes; no edits applied")
    noted = [k in current for k in notes]
    if not present and not all(noted):
        raise ValueError(f"{BIOMASS} lacks {TERM} without this curation; no edits applied")
    if present and any(noted):
        raise ValueError(f"{BIOMASS} carries this curation's notes with the {TERM} term; no edits applied")
    record = {"item": f"{BIOMASS} {TERM}", "curation_id": spec["curation_id"], "before_state": before["id"],
              "evidence_status": spec["status"]}
    if not present:
        return {**record, "status": "already_correct"}
    biomass.add_metabolites({metabolite: 0.0}, combine=False)  # coefficient 0 removes the entry; context-aware
    context = get_context(biomass)
    if context:
        # Notes are not history-tracked by cobra; restore them on context exit like the stoichiometry.
        context(partial(setattr, biomass, "notes", biomass.notes))
    biomass.notes = {**biomass.notes, **notes}
    after = {m.id: c for m, c in biomass.metabolites.items()}
    if TERM in after or after != stoichiometry:
        raise RuntimeError("Trehalose biomass removal did not reach the curated state")
    return {**record, "status": "applied", "coefficient_before": coefficient}
