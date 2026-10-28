"""Opt-in join of the split dolichol route (gap-census UM22 + UM167, 2026-10-28).

The model makes dehydrodolichyl diphosphate (R2180, Rer2/Srt1) and has dolichol kinase (R287, Sec59), but
nothing turns the one into dolichol: in between sits R2178, a gene-less, un-annotated placeholder that burns
O2 to CO2. The decision, data in ``dolichol_route.json``, is three edits made together:

- R2178 is closed (bounds 0, 0) with the reason in its notes; it is kept, not deleted.
- R_DOLSYN, one lumped gene-less step dehydrodolichyl-PP -> dolichol with the net stoichiometry of the yeast
  dephosphorylation and Tda5 / Dfg10 / Env9 detour, is added to the N-Glycan biosynthesis group.
- Dolichyl phosphate m503[C_er] is drained by biomass_C, so the carrier is made and not only cycled.

The step fails closed unless R2178, R287, the route species and biomass_C are exactly the reviewed ce9e1112
state (so it runs after every other step that edits biomass_C), and the three edits are either all absent or
all present with this curation's definition. A rerun, also on a saved and reloaded model, is a no-op.
"""

import copy
import json
import math

from cobra import Reaction

from .config import MODEL

CURATION_PATH = MODEL.curation_file("dolichol_route.json")
CURATION_ID = "DOLICHOL-ROUTE-20261028"
CLOSED, KINASE, NEW, BIOMASS, TERM = "R2178", "R287", "R_DOLSYN", "biomass_C", "m503[C_er]"
# The one reviewed starting state of biomass_C: candidate ce9e1112 (2026-10-28).
BEFORE_STATE_ID = "ce9e1112_biomass_C"


def _stoichiometry(reaction):
    return {m.id: c for m, c in reaction.metabolites.items()}


def _groups(model, reaction):
    return sorted(g.name for g in model.groups if reaction in g.members)


def _same(value, expected):
    return math.isclose(value, expected, rel_tol=1e-12, abs_tol=0.0)


def _scope_ok(spec):
    closed, kept, new, term = (spec.get(k) or {} for k in ("closed_reaction", "kept_reaction", "new_reaction", "biomass_term"))
    before = term.get("accepted_before_state") or {}
    return (spec.get("schema_version") == 1
            and spec.get("curation_id") == CURATION_ID
            and spec.get("status") == "user_authorized"
            and spec.get("requires_objective") == BIOMASS
            and closed.get("reaction_id") == CLOSED
            and closed.get("bounds") == [0.0, 1000.0]
            and closed.get("after_bounds") == [0.0, 0.0]
            and kept.get("reaction_id") == KINASE
            and new.get("reaction_id") == NEW
            and new.get("bounds") == [0.0, 1000.0]
            and new.get("gpr") == ""
            and term.get("biomass_reaction") == BIOMASS
            and term.get("metabolite_id") == TERM
            and isinstance(term.get("coefficient_mmol_per_gDW"), (int, float))
            and term["coefficient_mmol_per_gDW"] > 0
            and before.get("id") == BEFORE_STATE_ID
            and isinstance(before.get("biomass_stoichiometry"), dict)
            and TERM not in before["biomass_stoichiometry"]
            and isinstance(before.get("biomass_bounds"), list)
            and isinstance(before.get("biomass_note_keys"), list)
            and all(isinstance(row.get("notes"), dict) and row["notes"] for row in (closed, new, term))
            # Every species the step writes must have a reviewed identity row, so no lookup can fail after an edit.
            and isinstance(spec.get("species"), dict)
            and isinstance(new.get("stoichiometry"), dict)
            and set(new["stoichiometry"]) | {TERM} <= set(spec["species"]))


def apply_dolichol_route(model, spec=None):
    """Close R2178, add R_DOLSYN and the Dol-P biomass term; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if not _scope_ok(spec):
        raise ValueError("Unexpected dolichol route curation scope/status")
    closed, kept, new, term = (spec[k] for k in ("closed_reaction", "kept_reaction", "new_reaction", "biomass_term"))
    before = term["accepted_before_state"]
    for rid in (CLOSED, KINASE, BIOMASS):
        if rid not in model.reactions:
            raise ValueError(f"{rid} missing; no edits applied")
    if {r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient} != {BIOMASS: 1.0}:
        raise ValueError(f"The objective is not {BIOMASS} alone; no edits applied")
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
        # A route species that gained a producer or consumer elsewhere would make this decision unreviewed.
        if ("reactions_before" in row
                and sorted(r.id for r in metabolite.reactions if r.id not in (NEW, BIOMASS)) != row["reactions_before"]):
            raise ValueError(f"{mid} reaction context differs; no edits applied")
    reaction = model.reactions.get_by_id(CLOSED)
    if (reaction.name != closed["reaction_name"] or reaction.gene_reaction_rule != closed["gpr"]
            or reaction.annotation != closed["annotation"] or _groups(model, reaction) != closed["groups"]
            or _stoichiometry(reaction) != closed["stoichiometry"]
            or list(reaction.bounds) not in (closed["bounds"], closed["after_bounds"])):
        raise ValueError(f"{CLOSED} differs from the reviewed reaction; no edits applied")
    kinase = model.reactions.get_by_id(KINASE)
    if (kinase.name != kept["reaction_name"] or kinase.gene_reaction_rule != kept["gpr"]
            or list(kinase.bounds) != kept["bounds"] or kinase.annotation != kept["annotation"]
            or _groups(model, kinase) != kept["groups"] or _stoichiometry(kinase) != kept["stoichiometry"]):
        raise ValueError(f"{KINASE} is not the reviewed dolichol kinase; no edits applied")
    found = [g for g in model.groups if g.name == new["group"]]
    if len(found) != 1:
        raise ValueError(f"Group {new['group']!r} is missing or not unique; no edits applied")
    group = found[0]
    biomass = model.reactions.get_by_id(BIOMASS)
    stoichiometry = _stoichiometry(biomass)
    others = {mid: c for mid, c in stoichiometry.items() if mid != TERM}
    expected = before["biomass_stoichiometry"]
    if (list(biomass.bounds) != before["biomass_bounds"] or set(others) != set(expected)
            or not all(_same(others[mid], c) for mid, c in expected.items())):
        raise ValueError(f"{BIOMASS} composition differs from the reviewed baseline; no edits applied")
    if sorted(key for key in biomass.notes if key not in term["notes"]) != sorted(before["biomass_note_keys"]):
        raise ValueError(f"{BIOMASS} notes outside this curation differ; no edits applied")
    for target, notes, label in ((reaction, closed["notes"], CLOSED), (biomass, term["notes"], BIOMASS)):
        if any(key in target.notes and target.notes[key] != value for key, value in notes.items()):
            raise ValueError(f"Conflicting {label} dolichol route notes; no edits applied")
    # Each edit is either absent or present exactly as this curation writes it.
    closed_done = list(reaction.bounds) == closed["after_bounds"]
    closed_noted = any(key in reaction.notes for key in closed["notes"])
    if closed_done and not all(key in reaction.notes for key in closed["notes"]):
        raise ValueError(f"{CLOSED} is already closed without this curation; no edits applied")
    if not closed_done and closed_noted:
        raise ValueError(f"{CLOSED} is open but carries this curation's notes; no edits applied")
    term_present = TERM in stoichiometry
    if term_present and not _same(stoichiometry[TERM], -term["coefficient_mmol_per_gDW"]):
        raise ValueError(f"Conflicting {TERM} biomass coefficient; no edits applied")
    if term_present and not all(key in biomass.notes for key in term["notes"]):
        raise ValueError(f"{BIOMASS} already carries {TERM} without this curation; no edits applied")
    if not term_present and any(key in biomass.notes for key in term["notes"]):
        raise ValueError(f"{BIOMASS} carries this curation's notes without the {TERM} term; no edits applied")
    new_present = NEW in model.reactions
    if new_present:
        added = model.reactions.get_by_id(NEW)
        if ((added.name, list(added.bounds), added.gene_reaction_rule, _stoichiometry(added), added.annotation,
             added.notes, _groups(model, added))
                != (new["reaction_name"], new["bounds"], new["gpr"], new["stoichiometry"], new["annotation"],
                    new["notes"], [new["group"]])):
            raise ValueError(f"{NEW} already exists with another definition; no edits applied")
    present = {closed_done, term_present, new_present}
    if present == {True, False}:
        raise ValueError("Only part of the dolichol route curation is present; no edits applied")
    applied = present == {False}
    if applied:
        reaction.bounds = tuple(closed["after_bounds"])
        reaction.notes = {**reaction.notes, **closed["notes"]}
        added = Reaction(NEW, name=new["reaction_name"], subsystem=new["group"],
                         lower_bound=new["bounds"][0], upper_bound=new["bounds"][1])
        model.add_reactions([added])
        added.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in sorted(new["stoichiometry"].items())})
        added.gene_reaction_rule = new["gpr"]
        added.annotation = copy.deepcopy(new["annotation"])
        added.notes = dict(new["notes"])
        group.add_members([added])
        biomass.add_metabolites({model.metabolites.get_by_id(TERM): -term["coefficient_mmol_per_gDW"]})
        biomass.notes = {**biomass.notes, **term["notes"]}
    return {"item": "dolichol route", "curation_id": spec["curation_id"],
            "status": "applied" if applied else "already_correct",
            "closed": CLOSED, "added_reaction": NEW, "biomass_term": f"{TERM} {term['coefficient_mmol_per_gDW']:g}",
            "before_state": before["id"], "evidence_status": spec["status"]}
