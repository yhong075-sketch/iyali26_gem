"""Opt-in methylcitrate cycle curation (2026-11-16): matrix 2-oxobutanoate dehydrogenase added, R2075 closed.

On cbf2ef95 the methylcitrate cycle (R25, R95, R552) is blocked: matrix propionyl-CoA m78 has no producer that
leaves free CoA in the matrix. The decision, data in ``methylcitrate_cycle_20261116.json``, makes two edits:

- R_OBDHm (new): 2-oxobutanoate + CoA + NAD -> propionyl-CoA + CO2 + NADH + 4 H+ in C_mi, bounds (0, 1000), gene
  rule copied from R693 (pyruvate dehydrogenase). The rule is INFERRED (user: "just infer"); the data and the
  reaction notes say so;
- R2075 (gene-less methylmalonyl-CoA mutase): closed to (0, 0) and kept, with notes. Without the closure the new
  reaction opens a gene-less side route through R2113, R2112, R2074, R2076, R2075 and R2077.

Every check runs before any edit. Both targets must be in their reviewed cbf2ef95 state or both in the curated state
(a rerun, also on a saved and reloaded model, is a no-op). The species, R693, the low-degree species' reaction lists,
the five PDH genes' reaction links and the absence of any other reaction with the new chemistry are pinned.
"""

import copy
import json
from functools import partial
from html import unescape

from cobra import Reaction
from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("methylcitrate_cycle_20261116.json")
CURATION_ID = "METHYLCITRATE-CYCLE-20261116"
NEW_ID = "R_OBDHm"
CLOSED_ID = "R2075"
GPR_SOURCE = "R693"
PROTON = "m28[C_mi]"
SPECIES = ("m42[C_mi]", "m73[C_mi]", "m27[C_mi]", "m78[C_mi]", "m44[C_mi]", "m30[C_mi]", "m28[C_mi]",
           "m2047[C_cy]", "m2046[C_cy]")
CONTEXT = ("m42[C_mi]", "m78[C_mi]", "m2047[C_cy]", "m2046[C_cy]")
GENES = ("YALI1E31947g", "YALI1F27556g", "YALI1D26360g", "YALI1B13179g", "YALI1D30618g")
NEW_STOICHIOMETRY = {"m42[C_mi]": -1.0, "m73[C_mi]": -1.0, "m27[C_mi]": -1.0,
                     "m78[C_mi]": 1.0, "m44[C_mi]": 1.0, "m30[C_mi]": 1.0, "m28[C_mi]": 4.0}


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _stoichiometry(reaction):
    return {m.id: float(c) for m, c in reaction.metabolites.items()}


def _groups(model, reaction):
    return sorted(g.name for g in model.groups if reaction in g.members)


def _check_spec(spec):
    new, closed, source = spec.get("new_reaction", {}), spec.get("closed_reaction", {}), spec.get("gpr_source", {})
    genes = spec.get("genes", {})
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != CURATION_ID
            or spec.get("status") != "user_authorized"
            or new.get("reaction_id") != NEW_ID
            or new.get("bounds") != [0.0, 1000.0]
            or new.get("stoichiometry") != NEW_STOICHIOMETRY
            or new.get("gpr_status") != "inferred"
            or new.get("gpr_source_reaction") != GPR_SOURCE
            or source.get("reaction_id") != GPR_SOURCE
            or not new.get("gpr") or new.get("gpr") != source.get("gpr")
            or not new.get("name") or not new.get("subsystem") or not new.get("group")
            or not isinstance(new.get("annotation"), dict) or not new["annotation"]
            or not isinstance(new.get("notes"), dict) or not new["notes"]
            or "INFERRED" not in new["notes"].get("methylcitrate_cycle_gpr_status", "")
            or closed.get("reaction_id") != CLOSED_ID
            or closed.get("before_bounds") != [0.0, 1000.0]
            or closed.get("after_bounds") != [0.0, 0.0]
            or any(k not in closed for k in ("name", "gpr", "stoichiometry", "annotation", "before_notes"))
            or not isinstance(closed.get("notes"), dict) or not closed["notes"]
            or closed["notes"].get("methylcitrate_cycle_previous_bounds") != str(closed["before_bounds"])
            or set(closed["notes"]) & set(closed["before_notes"])
            or set(spec.get("species", {})) != set(SPECIES)
            or set(spec.get("species_reactions_before", {})) != set(CONTEXT)
            or set(genes) != set(GENES)
            or any(sorted(g.get("reactions_before", []) + [NEW_ID]) != g.get("reactions_after") for g in genes.values())
            # XML-special characters are escaped or dropped by the SBML writer, so the notes would not roundtrip.
            or any(not isinstance(text, str) or not text or any(ch in text for ch in "<>&'\"")
                   for text in [*new["notes"].values(), *closed["notes"].values()])):
        raise ValueError("Unexpected methylcitrate cycle curation scope/status")


def _closed_state(model, closed):
    """True when R2075 is in its curated state, False in its reviewed starting state."""
    if CLOSED_ID not in model.reactions:
        raise ValueError(f"{CLOSED_ID} missing; no edits applied")
    reaction = model.reactions.get_by_id(CLOSED_ID)
    if (reaction.name != closed["name"] or reaction.gene_reaction_rule != closed["gpr"]
            or _stoichiometry(reaction) != closed["stoichiometry"]
            or dict(reaction.annotation) != closed["annotation"]):
        raise ValueError(f"{CLOSED_ID} methylcitrate cycle precondition differs; no edits applied")
    bounds, notes = list(reaction.bounds), _notes(reaction.notes)
    if bounds == closed["before_bounds"] and notes == closed["before_notes"]:
        return False
    if bounds == closed["after_bounds"] and notes == {**closed["before_notes"], **closed["notes"]}:
        return True
    if bounds == closed["after_bounds"] and notes == closed["before_notes"]:
        # Closed by something else: the previous-bounds note would be false.
        raise ValueError(f"{CLOSED_ID} is already closed without this curation; no edits applied")
    raise ValueError(f"{CLOSED_ID} bounds or notes differ from the reviewed states; no edits applied")


def _new_state(model, new):
    """True when R_OBDHm exists with the curated definition, False when absent."""
    if NEW_ID not in model.reactions:
        return False
    added = model.reactions.get_by_id(NEW_ID)
    if ((added.name, added.subsystem, list(added.bounds), added.gene_reaction_rule, _stoichiometry(added),
         dict(added.annotation), _notes(added.notes), _groups(model, added))
            != (new["name"], new["subsystem"], new["bounds"], new["gpr"], new["stoichiometry"], new["annotation"],
                new["notes"], [new["group"]])):
        raise ValueError(f"{NEW_ID} already exists with another definition; no edits applied")
    return True


def _same_chemistry(reaction):
    """True when a reaction has the new chemistry, either direction, with any proton coefficient."""
    core = {k: v for k, v in _stoichiometry(reaction).items() if k != PROTON}
    target = {k: v for k, v in NEW_STOICHIOMETRY.items() if k != PROTON}
    return core == target or core == {k: -v for k, v in target.items()}


def apply_methylcitrate_cycle(model, spec=None):
    """Add R_OBDHm (rule as R693, inferred) and close R2075; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    new, closed, source = spec["new_reaction"], spec["closed_reaction"], spec["gpr_source"]

    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
    if GPR_SOURCE not in model.reactions:
        raise ValueError(f"{GPR_SOURCE} missing; no edits applied")
    pdh = model.reactions.get_by_id(GPR_SOURCE)
    if ((pdh.name, list(pdh.bounds), pdh.gene_reaction_rule, _stoichiometry(pdh))
            != (source["name"], source["bounds"], source["gpr"], source["stoichiometry"])):
        raise ValueError(f"{GPR_SOURCE} is not the pyruvate dehydrogenase this curation copied; no edits applied")

    states = {CLOSED_ID: _closed_state(model, closed), NEW_ID: _new_state(model, new)}
    if len(set(states.values())) != 1:
        raise ValueError(f"Only part of the methylcitrate cycle curation is present {states}; no edits applied")
    applied = states[NEW_ID]

    duplicates = sorted(r.id for r in model.reactions if r.id != NEW_ID and _same_chemistry(r))
    if duplicates:
        raise ValueError(f"Reactions with the new chemistry already exist {duplicates}; no edits applied")
    for mid, reactions in spec["species_reactions_before"].items():
        expected = sorted(reactions + ([NEW_ID] if applied and mid in new["stoichiometry"] else []))
        if sorted(r.id for r in model.metabolites.get_by_id(mid).reactions) != expected:
            raise ValueError(f"{mid} reaction context differs; no edits applied")
    for gene_id, row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        expected = row["reactions_after"] if applied else row["reactions_before"]
        if sorted(r.id for r in model.genes.get_by_id(gene_id).reactions) != expected:
            raise ValueError(f"{gene_id} reaction context differs; no edits applied")
    found = [g for g in model.groups if g.name == new["group"]]
    if len(found) != 1:
        raise ValueError(f"Group {new['group']!r} is missing or not unique; no edits applied")

    record = {"item": [NEW_ID, CLOSED_ID], "curation_id": spec["curation_id"], "evidence_status": spec["status"],
              "gpr_status": new["gpr_status"]}
    if applied:
        return {**record, "status": "already_correct"}

    reaction = model.reactions.get_by_id(CLOSED_ID)
    context = get_context(reaction)
    if context:
        # Notes are not history-tracked by cobra; restore them on context exit like the bounds.
        context(partial(setattr, reaction, "notes", reaction.notes))
    reaction.bounds = tuple(closed["after_bounds"])  # history-tracked by cobra
    reaction.notes = {**reaction.notes, **closed["notes"]}

    added = Reaction(NEW_ID, name=new["name"], subsystem=new["subsystem"],
                     lower_bound=new["bounds"][0], upper_bound=new["bounds"][1])
    model.add_reactions([added])  # history-tracked by cobra
    added.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in sorted(new["stoichiometry"].items())})
    added.gene_reaction_rule = new["gpr"]
    added.annotation = copy.deepcopy(new["annotation"])
    added.notes = dict(new["notes"])
    group = found[0]
    group.add_members([added])
    context = get_context(model)
    if context:
        context(partial(group.remove_members, [added]))

    if (added.check_mass_balance() != {}
            or not _new_state(model, new)
            or not _closed_state(model, closed)
            or any(sorted(r.id for r in model.genes.get_by_id(g).reactions) != spec["genes"][g]["reactions_after"]
                   for g in GENES)):
        raise RuntimeError("Methylcitrate cycle curation did not reach the curated state")
    return {**record, "status": "applied", "added_reactions": [NEW_ID], "closed_reactions": [CLOSED_ID],
            "gpr": new["gpr"]}
