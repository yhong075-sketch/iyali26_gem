"""Opt-in nucleus-cytosol diffusion of ADP-ribose and diphosphate (gap-census UM7 + UM171, 2026-10-28).

The model keeps the iYli21 nuclear copies of NAD synthesis (R569, R579) and of the sirtuin step R567, but their
products ADP-ribose m818[C_nu] and diphosphate m210[C_nu] have no consumer and no transport, so the whole
nuclear NAD block is dead; the cytosolic ADP-ribose m1725[C_cy] has a consumer (R1912) and no producer. The
decision, data in ``um7_nuclear_diffusion.json``, adds two gene-less reversible diffusions in the
"Transport, nucleus" group, like the 33 nuclear diffusions the model already has:

- R_ADPRtn: ADP-ribose m1725[C_cy] <=> m818[C_nu];
- R_PPItn: diphosphate m203[C_cy] <=> m210[C_nu].

The step fails closed unless the four species are the reviewed 9d655c2b ones (identity and, for the three
dead-end ends, their reaction context, so no other reaction can already link a pair), the group exists once,
and the two reactions are either both absent or both present with this curation's definition. A rerun, also on a saved
and reloaded model, is a no-op.
"""

import copy
import json

from cobra import Reaction

from .config import MODEL

CURATION_PATH = MODEL.curation_file("um7_nuclear_diffusion.json")
CURATION_ID = "UM7-NUCLEAR-DIFFUSION-20261028"
NEW_IDS = ("R_ADPRtn", "R_PPItn")
GROUP = "Transport, nucleus"


def _stoichiometry(reaction):
    return {m.id: c for m, c in reaction.metabolites.items()}


def _groups(model, reaction):
    return sorted(g.name for g in model.groups if reaction in g.members)


def _scope_ok(spec):
    rows = spec.get("new_reactions")
    species = spec.get("species")
    if not (isinstance(rows, list) and isinstance(species, dict)):
        return False
    return (spec.get("schema_version") == 1
            and spec.get("curation_id") == CURATION_ID
            and spec.get("status") == "user_authorized"
            and spec.get("group") == GROUP
            and tuple(row.get("reaction_id") for row in rows) == NEW_IDS
            and all(row.get("bounds") == [-1000.0, 1000.0]
                    and row.get("gpr") == ""
                    and isinstance(row.get("notes"), dict) and row["notes"]
                    and isinstance(row.get("stoichiometry"), dict)
                    # One species leaves the cytosol, its nuclear copy appears: a pure transport.
                    and sorted(row["stoichiometry"].values()) == [-1.0, 1.0]
                    and set(row["stoichiometry"]) <= set(species)
                    for row in rows))


def apply_um7_nuclear_diffusion(model, spec=None):
    """Add R_ADPRtn and R_PPItn; fail closed on any unexpected model state."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if not _scope_ok(spec):
        raise ValueError("Unexpected UM7 nuclear diffusion curation scope/status")
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid} missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        if any(getattr(metabolite, key) != row[key] for key in ("name", "formula", "charge", "compartment")):
            raise ValueError(f"{mid} identity differs; no edits applied")
        # A dead end that gained a producer or consumer elsewhere would make this decision unreviewed.
        if ("reactions_before" in row
                and sorted(r.id for r in metabolite.reactions if r.id not in NEW_IDS) != row["reactions_before"]):
            raise ValueError(f"{mid} reaction context differs; no edits applied")
    found = [g for g in model.groups if g.name == GROUP]
    if len(found) != 1:
        raise ValueError(f"Group {GROUP!r} is missing or not unique; no edits applied")
    group = found[0]
    present = []
    for row in spec["new_reactions"]:
        rid = row["reaction_id"]
        ends = [model.metabolites.get_by_id(mid) for mid in row["stoichiometry"]]
        # Both ends must be the same compound, so the transport is mass and charge balanced.
        if len({(m.formula, m.charge) for m in ends}) != 1 or len({m.compartment for m in ends}) != 2:
            raise ValueError(f"{rid} ends are not one compound in two compartments; no edits applied")
        if rid in model.reactions:
            added = model.reactions.get_by_id(rid)
            if ((added.name, list(added.bounds), added.gene_reaction_rule, _stoichiometry(added), added.annotation,
                 added.notes, _groups(model, added))
                    != (row["reaction_name"], row["bounds"], row["gpr"], row["stoichiometry"], row["annotation"],
                        row["notes"], [GROUP])):
                raise ValueError(f"{rid} already exists with another definition; no edits applied")
        present.append(rid in model.reactions)
    if len(set(present)) != 1:
        raise ValueError("Only part of the UM7 nuclear diffusion curation is present; no edits applied")
    applied = not present[0]
    if applied:
        for row in spec["new_reactions"]:
            added = Reaction(row["reaction_id"], name=row["reaction_name"], subsystem=GROUP,
                             lower_bound=row["bounds"][0], upper_bound=row["bounds"][1])
            model.add_reactions([added])
            added.add_metabolites({model.metabolites.get_by_id(mid): c
                                   for mid, c in sorted(row["stoichiometry"].items())})
            added.gene_reaction_rule = row["gpr"]
            added.annotation = copy.deepcopy(row["annotation"])
            added.notes = dict(row["notes"])
            group.add_members([added])
    return {"item": "UM7 nuclear diffusion", "curation_id": spec["curation_id"],
            "status": "applied" if applied else "already_correct",
            "added_reactions": list(NEW_IDS), "evidence_status": spec["status"]}
