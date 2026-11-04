"""Opt-in: 17 reviewed reactions from Andrew McNaughton's mcnaugh_pass model, closing the lumped steps they replace."""

import copy
import json

import cobra
from cobra.core.gene import GPR

from .config import MODEL

CURATION_PATH = MODEL.curation_file("mcnaugh_reactions_20261103.json")


def _stoichiometry(reaction):
    return {m.id: c for m, c in reaction.metabolites.items()}


def _groups(model, reaction):
    return sorted(g.name for g in model.groups if reaction in g.members)


def _definition(model, reaction):
    return (reaction.name, _stoichiometry(reaction), list(reaction.bounds), reaction.gene_reaction_rule,
            reaction.annotation, reaction.notes, _groups(model, reaction))


def apply_mcnaugh_reactions(model, spec=None):
    """Add the chain-specific Kennedy-pathway and biomass reactions and close R353, R1846, R1771 and R2165.

    Every precondition is checked before any edit: on a stale, partly applied or unexpected model nothing changes.
    """
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1 or spec.get("change_id") != "MCNAUGH-REACTIONS-20261103"
            or spec.get("status") != "user_authorized" or len(spec.get("reactions", [])) != 17
            or len(spec.get("closures", [])) != 4):
        raise ValueError("Unexpected mcnaugh_pass reaction set scope/status")
    if [r.id for r in model.reactions if r.objective_coefficient] != [spec["objective_unchanged"]]:
        raise ValueError(f"Objective is not {spec['objective_unchanged']} alone; no edits applied")
    groups = {}
    for name in sorted({g for row in spec["reactions"] for g in row["groups"]}):
        found = [g for g in model.groups if g.name == name]
        if len(found) != 1:
            raise ValueError(f"Group {name!r} is missing or not unique; no edits applied")
        groups[name] = found[0]
    # A new reaction must not duplicate an existing one under another ID, written in either direction.
    existing = {}
    for r in model.reactions:
        existing.setdefault(frozenset(_stoichiometry(r).items()), set()).add(r.id)
    states = []
    for row in spec["reactions"]:
        missing = [m for m in row["stoichiometry"] if m not in model.metabolites]
        missing += [g for g in GPR.from_string(row["gpr"]).genes if g not in model.genes]
        if missing:
            raise ValueError(f"{row['id']}: species or genes missing {missing[:5]}; no edits applied")
        forward = frozenset(row["stoichiometry"].items())
        reverse = frozenset((m, -c) for m, c in row["stoichiometry"].items())
        clash = (existing.get(forward, set()) | existing.get(reverse, set())) - {row["id"]}
        if clash:
            raise ValueError(f"{row['id']}: {', '.join(sorted(clash))} has the same stoichiometry or its reverse; "
                             f"no edits applied")
        if row["id"] in model.reactions:
            if (_definition(model, model.reactions.get_by_id(row["id"]))
                    != (row["name"], row["stoichiometry"], row["bounds"], row["gpr"], row["annotation"], row["notes"],
                        row["groups"])):
                raise ValueError(f"{row['id']} exists with another definition; no edits applied")
            states.append("after")
        else:
            states.append("before")
    for row in spec["closures"]:
        if row["id"] not in model.reactions:
            raise ValueError(f"{row['id']} missing; no edits applied")
        reaction = model.reactions.get_by_id(row["id"])
        # Stoichiometry, not the equation string: the string's arrow changes with the bounds this step sets.
        if _stoichiometry(reaction) != row["stoichiometry"] or reaction.gene_reaction_rule != row["gpr"]:
            raise ValueError(f"{row['id']} stoichiometry or GPR differs from the reviewed state; no edits applied")
        bounds = list(reaction.bounds)
        states.append("before" if bounds == row["before_bounds"] else "after" if bounds == row["after_bounds"] else "stale")
    if set(states) not in ({"before"}, {"after"}):
        raise ValueError(f"mcnaugh_pass reactions partly applied or stale ({sorted(set(states))}); no edits applied")
    if set(states) == {"before"}:
        new = []
        for row in spec["reactions"]:
            # cobra writes a subsystem only through group membership, so the subsystem is the group's name.
            reaction = cobra.Reaction(row["id"], name=row["name"], subsystem=row["groups"][0] if row["groups"] else "",
                                      lower_bound=row["bounds"][0], upper_bound=row["bounds"][1])
            reaction.add_metabolites({model.metabolites.get_by_id(m): c for m, c in row["stoichiometry"].items()})
            new.append(reaction)
        model.add_reactions(new)
        for row, reaction in zip(spec["reactions"], new):
            reaction.gene_reaction_rule = row["gpr"]
            reaction.annotation = copy.deepcopy(row["annotation"])
            reaction.notes = dict(row["notes"])
            for name in row["groups"]:
                groups[name].add_members([reaction])
        for row in spec["closures"]:
            model.reactions.get_by_id(row["id"]).bounds = tuple(row["after_bounds"])
    return {"item": spec["change_id"], "status": "applied" if set(states) == {"before"} else "already_correct",
            "added": [r["id"] for r in spec["reactions"]], "closed": [r["id"] for r in spec["closures"]]}
