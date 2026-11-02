"""Opt-in ceramide pool term in biomass_C (2026-11-01).

biomass_C has no sphingolipid term, so the sphingoid backbone (serine palmitoyltransferase R722, 3-ketosphinganine
reductase R37, ceramide synthesis and hydroxylation) carries no flux in WT. Y. lipolytica contains about 5.4e-4 mmol
sphingolipid per gDW (Rupcic 1996, total long-chain bases; cross-checked against Ianutsevich 2019). One
pseudo-metabolite sphingolipid_ceramide_pool[C_er] is filled by gene-less pseudo-reactions xPOOL_CER_* from each of
the 10 ER ceramides, so dihydro, phyto and alpha-hydroxy forms all count (a Ylsur2 deletion is viable), and
biomass_C consumes the pool. The decision and its evidence are data in ``ceramide_pool_biomass.json``; biomass_C is
not re-normalised and no other reaction changes.

The step fails closed unless biomass_C is the maximised objective and equals the reviewed c1f70ac4 state (every
coefficient, the bounds and the note keys; note values are not pinned), so it runs after every other step that edits
biomass_C, and unless the 10 ceramides keep their reviewed identity. The pool, its 10 reactions, the term and this
curation's notes are either all absent or all present as curated; a rerun, also on a saved and reloaded model, is a
no-op.
"""

import json
import math
from functools import partial
from html import unescape

from cobra import Metabolite, Reaction
from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("ceramide_pool_biomass.json")
CURATION_ID = "CERAMIDE-POOL-BIOMASS-20261101"
BIOMASS, POOL, COEFFICIENT = "biomass_C", "sphingolipid_ceramide_pool[C_er]", -5.4e-4
CERAMIDES = ("m386", "m388", "m390", "m391", "m392", "m393", "m397", "m398", "m399", "m400")
REACTION_PREFIX = "xPOOL_CER_"


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _same(value, expected):
    return math.isclose(value, expected, rel_tol=1e-12, abs_tol=0.0)


def _check_spec(spec):
    pool, before = spec.get("pool") or {}, spec.get("accepted_before_state") or {}
    feeds = pool.get("feeds") or []
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != CURATION_ID
            or spec.get("status") != "user_authorized_biomass_term_addition"
            or spec.get("requires_objective") != BIOMASS
            or pool.get("metabolite_id") != POOL
            or pool.get("compartment") != "C_er"
            or pool.get("biomass_reaction") != BIOMASS
            or pool.get("coefficient") != COEFFICIENT
            or [f.get("ceramide", {}).get("id") for f in feeds] != [f"{c}[C_er]" for c in CERAMIDES]
            or [f.get("reaction_id") for f in feeds] != [f"{REACTION_PREFIX}{c}" for c in CERAMIDES]
            or any(f.get("bounds") != [0.0, 1000.0] for f in feeds)
            or not isinstance(spec.get("biomass_notes"), dict) or not spec["biomass_notes"]
            or not isinstance(before.get("biomass_stoichiometry"), dict)
            or POOL in before["biomass_stoichiometry"]
            or not isinstance(before.get("biomass_bounds"), list)
            or not isinstance(before.get("biomass_note_keys"), list)):
        raise ValueError("Unexpected ceramide pool biomass scope/status")


def _pool_reaction_state(reaction):
    return (reaction.name, list(reaction.bounds), reaction.gene_reaction_rule,
            {m.id: c for m, c in reaction.metabolites.items()})


def apply_ceramide_pool_biomass(model, spec=None):
    """Add the ER ceramide pool and its biomass_C term; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    pool_spec, before, notes = spec["pool"], spec["accepted_before_state"], spec["biomass_notes"]
    if BIOMASS not in model.reactions:
        raise ValueError(f"{BIOMASS} missing; no edits applied")
    if ({r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient} != {BIOMASS: 1.0}
            or model.objective.direction != "max"):
        raise ValueError(f"The objective is not {BIOMASS} alone, maximised; no edits applied")
    for feed in pool_spec["feeds"]:
        ceramide = feed["ceramide"]
        if ceramide["id"] not in model.metabolites:
            raise ValueError(f"{ceramide['id']} missing; no edits applied")
        met = model.metabolites.get_by_id(ceramide["id"])
        if (met.name, met.formula, met.charge, met.compartment) != (
                ceramide["name"], ceramide["formula"], ceramide["charge"], ceramide["compartment"]):
            raise ValueError(f"{ceramide['id']} identity differs; no edits applied")

    biomass = model.reactions.get_by_id(BIOMASS)
    stoichiometry = {m.id: c for m, c in biomass.metabolites.items()}
    present = POOL in stoichiometry
    coefficient = stoichiometry.pop(POOL, None)
    expected = before["biomass_stoichiometry"]
    if (list(biomass.bounds) != before["biomass_bounds"] or set(stoichiometry) != set(expected)
            or any(not _same(stoichiometry[k], v) for k, v in expected.items())):
        raise ValueError(f"{BIOMASS} composition differs from the reviewed baseline; no edits applied")
    if present and not _same(coefficient, pool_spec["coefficient"]):
        raise ValueError(f"Conflicting {POOL} biomass coefficient; no edits applied")
    current = _notes(biomass.notes)
    if sorted(k for k in current if k not in notes) != sorted(before["biomass_note_keys"]):
        raise ValueError(f"{BIOMASS} notes outside this curation differ; no edits applied")
    if any(k in current and current[k] != v for k, v in notes.items()):
        raise ValueError(f"Conflicting {BIOMASS} ceramide pool notes; no edits applied")

    reaction_ids = [f["reaction_id"] for f in pool_spec["feeds"]]
    existing = sorted(r.id for r in model.reactions if r.id.startswith(REACTION_PREFIX))
    pool_present = POOL in model.metabolites
    noted = [k in current for k in notes]
    parts = [present, pool_present, bool(existing), *noted]
    record = {"item": f"{BIOMASS} {POOL}", "curation_id": spec["curation_id"], "before_state": before["id"],
              "evidence_status": spec["status"]}
    if all(parts):
        pool = model.metabolites.get_by_id(POOL)
        wanted = {f["reaction_id"]: (f["name"], f["bounds"], "", {f["ceramide"]["id"]: -1.0, POOL: 1.0})
                  for f in pool_spec["feeds"]}
        if (existing != sorted(reaction_ids)
                or any(_pool_reaction_state(model.reactions.get_by_id(r)) != wanted[r] for r in reaction_ids)
                or (pool.name, pool.formula, pool.charge, pool.compartment) != (
                    pool_spec["name"], None, pool_spec["charge"], pool_spec["compartment"])
                or sorted(r.id for r in pool.reactions) != sorted([*reaction_ids, BIOMASS])):
            raise ValueError("The ceramide pool differs from the curated state; no edits applied")
        return {**record, "status": "already_correct"}
    if any(parts):
        raise ValueError("Partial ceramide pool state (pool, reactions, term or notes); no edits applied")

    pool = Metabolite(POOL, name=pool_spec["name"], compartment=pool_spec["compartment"], charge=pool_spec["charge"])
    reactions = []
    for feed in pool_spec["feeds"]:
        reaction = Reaction(feed["reaction_id"], name=feed["name"], lower_bound=feed["bounds"][0],
                            upper_bound=feed["bounds"][1])
        reaction.add_metabolites({model.metabolites.get_by_id(feed["ceramide"]["id"]): -1.0, pool: 1.0})
        reactions.append(reaction)
    model.add_reactions(reactions)  # context-aware: removed again on context exit
    biomass.add_metabolites({model.metabolites.get_by_id(POOL): pool_spec["coefficient"]})
    context = get_context(biomass)
    if context:
        # Notes are not history-tracked by cobra; restore them on context exit like the stoichiometry.
        context(partial(setattr, biomass, "notes", biomass.notes))
    biomass.notes = {**biomass.notes, **notes}
    after = {m.id: c for m, c in biomass.metabolites.items()}
    if set(after) != set(stoichiometry) | {POOL} or not _same(after[POOL], pool_spec["coefficient"]):
        raise RuntimeError("Ceramide pool biomass term did not reach the curated state")
    return {**record, "status": "applied", "reactions_added": reaction_ids, "coefficient": pool_spec["coefficient"]}
