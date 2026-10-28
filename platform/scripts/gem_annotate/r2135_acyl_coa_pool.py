"""Opt-in R2135 acyl-CoA pool remap and R2145 1:1 rescale; stoichiometry and notes of the two reactions only."""

import json

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r2135_acyl_coa_pool.json")
REACTIONS = ("R2135", "R2145")
# The reviewed after-equations (2026-10-27); a data file that names anything else is out of scope.
AFTER = {"R2135": {"m567[C_em]": -1.0, "m1984[C_er]": 1.0},
         "R2145": {"m1993[C_cy]": -1.0, "m32[C_cy]": -1.0, "m1994[C_cy]": 1.0, "m35[C_cy]": 1.0}}


def _stoichiometry(reaction):
    return {m.id: c for m, c in reaction.metabolites.items()}


def _species_matches(model, mid, expected):
    met = model.metabolites.get_by_id(mid)
    if (met.formula, met.charge, met.compartment) != (expected["formula"], expected["charge"], expected["compartment"]):
        return False
    # The remap rests on m567 and m1984 sharing one generic acyl-CoA id; a changed id voids the review.
    if "metanetx.chemical" in expected and met.annotation.get("metanetx.chemical") != expected["metanetx.chemical"]:
        return False
    # m1984's consumers (R2134, R2168) are part of the evidence; a new consumer was not reviewed.
    return "reactions" not in expected or sorted(r.id for r in met.reactions) == expected["reactions"]


def apply_r2135_acyl_coa_pool(model, spec=None):
    """Replace acetyl-CoA by the ER-membrane acyl-CoA pool in R2135 and write R2145 1:1; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "R2135-ACYL-COA-POOL-20261027"
            or spec.get("status") != "user_authorized"
            or tuple(sorted(spec.get("reactions", {}))) != REACTIONS
            or any(spec["reactions"][rid].get("after_stoichiometry") != AFTER[rid] for rid in REACTIONS)
            or set(spec.get("species", {})) != {mid for row in spec["reactions"].values()
                                                for key in ("before_stoichiometry", "after_stoichiometry")
                                                for mid in row[key]}):
        raise ValueError("Unexpected R2135 acyl-CoA pool curation scope")
    if any(rid not in model.reactions for rid in REACTIONS) or any(mid not in model.metabolites for mid in spec["species"]):
        raise ValueError("R2135 / R2145 or one of their species missing; no edits applied")
    if not all(_species_matches(model, mid, row) for mid, row in spec["species"].items()):
        raise ValueError("A species of the R2135 / R2145 curation differs; no edits applied")
    states = {}
    for rid in REACTIONS:
        row, reaction = spec["reactions"][rid], model.reactions.get_by_id(rid)
        if (reaction.name != row["name"] or reaction.gene_reaction_rule != row["gpr"]
                or list(reaction.bounds) != row["bounds"]):
            raise ValueError(f"{rid} precondition differs; no edits applied")
        if any(k in reaction.notes and reaction.notes[k] != v for k, v in row["notes"].items()):
            raise ValueError(f"Conflicting {rid} curation notes; no edits applied")
        stoichiometry = _stoichiometry(reaction)
        if stoichiometry == row["before_stoichiometry"]:
            # This curation's notes on the old equation mean something reverted the edit.
            if any(k in reaction.notes for k in row["notes"]):
                raise ValueError(f"{rid} carries this curation's notes without its equation; no edits applied")
            states[rid] = "before"
        elif stoichiometry == row["after_stoichiometry"] and all(reaction.notes.get(k) == v for k, v in row["notes"].items()):
            states[rid] = "after"
        else:
            # Includes the after-equation without these notes: it came from somewhere else.
            raise ValueError(f"{rid} stoichiometry is neither the reviewed before- nor after-state; no edits applied")
    # The two edits were reviewed together: a half-applied state is not a state this curation knows.
    if len(set(states.values())) != 1:
        raise ValueError("Only one of R2135 / R2145 carries this curation; no edits applied")
    changed = states["R2135"] == "before"
    for rid in REACTIONS:
        row, reaction = spec["reactions"][rid], model.reactions.get_by_id(rid)
        if changed:
            metabolites = {model.metabolites.get_by_id(mid): c for mid, c in row["after_stoichiometry"].items()}
            for met in list(reaction.metabolites):
                if met.id not in row["after_stoichiometry"]:
                    reaction.subtract_metabolites({met: reaction.metabolites[met]})
            reaction.add_metabolites(metabolites, combine=False)
            reaction.notes = {**reaction.notes, **row["notes"]}
    return {"items": list(REACTIONS), "curation_id": spec["curation_id"],
            "status": "applied" if changed else "already_correct",
            "equations": {rid: model.reactions.get_by_id(rid).build_reaction_string() for rid in REACTIONS},
            "evidence_status": spec["status"]}
