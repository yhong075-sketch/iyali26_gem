"""Opt-in 5-formyl-THF (folinic acid) module curation (2026-11-12).

Cytosolic 5-formyl-THF m159 had no producer, so the FAU1-type cyclo-ligase R64 could never carry flux, and the
mitochondrial 5-formyl-THF source R748 carried the glycine-cleavage T protein gene (KEGG files R02300 under EC
2.1.2.10). The decision, data in ``folinic_acid_module.json``, makes three edits:

- R748 (mitochondrial methenyl-THF hydrolysis): rule to the SHMT genes, as on R359; name, EC and notes follow;
- R990 (gene-less lumped R64 + R543, "5-formyl-THF isomerase"): closed to (0, 0) and kept, with notes;
- R_THFAT (new): the same hydrolysis in the cytosol, methenyl-THF m160 + H2O -> m159 + H+, rule as R358.

Every check runs before any edit. The three targets must all be in their reviewed e584151f state or all in the
curated state (a rerun, also on a saved and reloaded model, is a no-op); the species, the folate group and the
three genes' reaction links are pinned in both states.
"""

import copy
import json
from functools import partial
from html import unescape

from cobra import Reaction
from cobra.util.context import get_context

from .config import MODEL

CURATION_PATH = MODEL.curation_file("folinic_acid_module.json")
CURATION_ID = "FOLINIC-ACID-MODULE-20261112"
TARGETS = ("R748", "R990")
NEW_ID = "R_THFAT"
GENES = ("YALI1D28569g", "YALI1E19557g", "YALI1F04082g")
SPECIES = ("m159[C_cy]", "m160[C_cy]", "m32[C_cy]", "m10[C_cy]")


def _notes(values):
    # SBML read escapes note text; compare the unescaped form.
    return {k: unescape(v) if isinstance(v, str) else v for k, v in values.items()}


def _groups(model, reaction):
    return sorted(g.name for g in model.groups if reaction in g.members)


def _definition(model, reaction):
    return {"name": reaction.name, "bounds": list(reaction.bounds), "gpr": reaction.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in reaction.metabolites.items()},
            "species": {m.id: {k: getattr(m, k) for k in ("formula", "charge", "compartment")}
                        for m in reaction.metabolites},
            "annotation": dict(reaction.annotation), "notes": _notes(reaction.notes),
            "groups": _groups(model, reaction)}


def _check_spec(spec):
    rows, new, genes = spec.get("reactions", {}), spec.get("new_reaction", {}), spec.get("genes", {})
    keys = {"name", "bounds", "gpr", "stoichiometry", "species", "annotation", "notes", "groups"}
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != CURATION_ID
            or spec.get("status") != "user_authorized"
            or set(rows) != set(TARGETS)
            or any(set(rows[r].get("before", {})) != keys or set(rows[r].get("after", {})) != keys for r in rows)
            # R748 keeps its chemistry and bounds; R990 changes only bounds and notes.
            or any(rows["R748"]["before"][k] != rows["R748"]["after"][k]
                   for k in ("bounds", "stoichiometry", "species", "groups"))
            or any(rows["R990"]["before"][k] != rows["R990"]["after"][k]
                   for k in ("name", "gpr", "stoichiometry", "species", "annotation", "groups"))
            or rows["R990"]["after"]["bounds"] != [0.0, 0.0]
            or new.get("reaction_id") != NEW_ID
            or new.get("bounds") != [0.0, 1000.0]
            or new.get("gpr") != rows["R748"]["after"]["gpr"]
            or set(new.get("stoichiometry", {})) != set(SPECIES)
            or not new.get("notes") or not new.get("annotation") or not new.get("group")
            or set(spec.get("species", {})) != set(SPECIES)
            or set(genes) != set(GENES)
            or any(not isinstance(g.get("reactions_before"), list) or not isinstance(g.get("reactions_after"), list)
                   for g in genes.values())):
        raise ValueError("Unexpected folinic acid module curation scope/status")


def _target_state(model, rid, row):
    """Return True when the target is in its curated state, False in its reviewed starting state."""
    if rid not in model.reactions:
        raise ValueError(f"{rid} missing; no edits applied")
    current = _definition(model, model.reactions.get_by_id(rid))
    if current == row["after"]:
        return True
    if current == row["before"]:
        return False
    raise ValueError(f"{rid} folinic acid module precondition differs; no edits applied")


def _new_state(model, new):
    if NEW_ID not in model.reactions:
        return False
    added = model.reactions.get_by_id(NEW_ID)
    if ((added.name, added.subsystem, list(added.bounds), added.gene_reaction_rule,
         {m.id: c for m, c in added.metabolites.items()}, dict(added.annotation), _notes(added.notes),
         _groups(model, added))
            != (new["name"], new["subsystem"], new["bounds"], new["gpr"], new["stoichiometry"], new["annotation"],
                new["notes"], [new["group"]])):
        raise ValueError(f"{NEW_ID} already exists with another definition; no edits applied")
    return True


def apply_folinic_acid_module(model, spec=None):
    """Move R748 to SHMT, close R990, add cytosolic R_THFAT; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    _check_spec(spec)
    states = {rid: _target_state(model, rid, spec["reactions"][rid]) for rid in TARGETS}
    states[NEW_ID] = _new_state(model, spec["new_reaction"])
    if len(set(states.values())) != 1:
        raise ValueError(f"Only part of the folinic acid module curation is present {states}; no edits applied")
    applied = states[NEW_ID]
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
        if "reactions_before" in row:
            expected = sorted(row["reactions_before"] + ([NEW_ID] if applied else []))
            if sorted(r.id for r in metabolite.reactions) != expected:
                raise ValueError(f"{mid} reaction context differs; no edits applied")
    for gene_id, row in spec["genes"].items():
        if gene_id not in model.genes:
            raise ValueError(f"{gene_id} missing; no edits applied")
        gene = model.genes.get_by_id(gene_id)
        expected = row["reactions_after"] if applied else row["reactions_before"]
        if dict(gene.annotation) != row["annotation"] or sorted(r.id for r in gene.reactions) != expected:
            raise ValueError(f"{gene_id} identity or reaction context differs; no edits applied")
    new = spec["new_reaction"]
    found = [g for g in model.groups if g.name == new["group"]]
    if len(found) != 1:
        raise ValueError(f"Group {new['group']!r} is missing or not unique; no edits applied")
    record = {"item": [*TARGETS, NEW_ID], "curation_id": spec["curation_id"], "evidence_status": spec["status"]}
    if applied:
        return {**record, "status": "already_correct"}

    for rid in TARGETS:
        reaction, after = model.reactions.get_by_id(rid), spec["reactions"][rid]["after"]
        context = get_context(reaction)
        if context:
            # Name, annotation and notes are not history-tracked by cobra; restore them on context exit.
            context(partial(setattr, reaction, "name", reaction.name))
            context(partial(setattr, reaction, "annotation", reaction.annotation))
            context(partial(setattr, reaction, "notes", reaction.notes))
        reaction.name = after["name"]
        reaction.bounds = tuple(after["bounds"])  # history-tracked by cobra
        reaction.gene_reaction_rule = after["gpr"]  # history-tracked by cobra
        reaction.annotation = copy.deepcopy(after["annotation"])
        reaction.notes = dict(after["notes"])
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

    if (any(_definition(model, model.reactions.get_by_id(rid)) != spec["reactions"][rid]["after"] for rid in TARGETS)
            or not _new_state(model, new)
            or any(sorted(r.id for r in model.genes.get_by_id(g).reactions) != spec["genes"][g]["reactions_after"]
                   for g in GENES)):
        raise RuntimeError("Folinic acid module did not reach the curated state")
    return {**record, "status": "applied", "added_reactions": [NEW_ID], "closed_reactions": ["R990"],
            "gpr_changed": {"R748": [spec["reactions"]["R748"]["before"]["gpr"], spec["reactions"]["R748"]["after"]["gpr"]]}}
