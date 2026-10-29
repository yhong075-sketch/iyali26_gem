"""Opt-in R2171 PI phospholipase C written 1:1 (was 100 H2O / 100 Ins1P / 100 H+); stoichiometry and notes only."""

import json

from .config import MODEL

CURATION_PATH = MODEL.curation_file("r2171_pi_plc_1to1.json")
# The reviewed before- and after-equations (2026-10-28, KEGG R03332 with 2 H+ for charge); anything else is out
# of scope. Both pinned: the rewrite uses combine=False, so a species outside AFTER would survive it.
BEFORE = {"m1393[C_cy]": -1.0, "m32[C_cy]": -100.0, "m10[C_cy]": 100.0, "m1994[C_cy]": 1.0, "m813[C_cy]": 100.0}
AFTER = {"m1393[C_cy]": -1.0, "m32[C_cy]": -1.0, "m10[C_cy]": 2.0, "m1994[C_cy]": 1.0, "m813[C_cy]": 1.0}


def _species_matches(model, mid, expected):
    met = model.metabolites.get_by_id(mid)
    if (met.formula, met.charge, met.compartment) != (expected["formula"], expected["charge"], expected["compartment"]):
        return False
    # The 2 H+ rest on these charges and the 1:1 link on these identities; a change voids the review.
    return all(met.annotation.get(key) == expected[key] for key in ("metanetx.chemical", "kegg.compound") if key in expected)


def apply_r2171_pi_plc_1to1(model, spec=None):
    """Rewrite R2171 as PI + H2O -> DAG + Ins1P + 2 H+; fail closed otherwise."""
    if spec is None:
        spec = json.loads(CURATION_PATH.read_text())
    notes = spec.get("notes")
    if (spec.get("schema_version") != 1
            or spec.get("curation_id") != "R2171-PI-PLC-1TO1-20261028"
            or spec.get("status") != "user_authorized"
            or spec.get("reaction_id") != "R2171"
            or spec.get("before_stoichiometry") != BEFORE
            or spec.get("after_stoichiometry") != AFTER
            or set(spec.get("species", {})) != set(AFTER)
            # Empty notes would void the "after-equation without these notes" guard; libSBML drops < > & in notes.
            or not isinstance(notes, dict) or not notes
            or any(not isinstance(v, str) or any(c in v for c in "<>&") for v in notes.values())):
        raise ValueError("Unexpected R2171 PI phospholipase C curation scope")
    if "R2171" not in model.reactions or any(mid not in model.metabolites for mid in spec["species"]):
        raise ValueError("R2171 or one of its species missing; no edits applied")
    if not all(_species_matches(model, mid, row) for mid, row in spec["species"].items()):
        raise ValueError("A species of the R2171 curation differs; no edits applied")
    reaction = model.reactions.get_by_id("R2171")
    if (reaction.name != spec["name"] or reaction.gene_reaction_rule != spec["gpr"]
            or list(reaction.bounds) != spec["bounds"] or reaction.annotation.get("kegg.reaction") != spec["kegg_reaction"]):
        raise ValueError("R2171 precondition differs; no edits applied")
    if any(k in reaction.notes and reaction.notes[k] != v for k, v in spec["notes"].items()):
        raise ValueError("Conflicting R2171 curation notes; no edits applied")
    stoichiometry = {m.id: c for m, c in reaction.metabolites.items()}
    if stoichiometry == BEFORE:
        # This curation's notes on the old equation mean something reverted the edit; do not undo that.
        if any(k in reaction.notes for k in spec["notes"]):
            raise ValueError("R2171 carries this curation's notes without its equation; no edits applied")
        reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in AFTER.items()}, combine=False)
        reaction.notes = {**reaction.notes, **spec["notes"]}
        status = "applied"
    elif stoichiometry == AFTER and all(reaction.notes.get(k) == v for k, v in spec["notes"].items()):
        status = "already_correct"
    else:
        # Includes the after-equation without these notes: it came from somewhere else.
        raise ValueError("R2171 stoichiometry is neither the reviewed before- nor after-state; no edits applied")
    return {"item": "R2171", "curation_id": spec["curation_id"], "status": status,
            "equation": reaction.build_reaction_string(), "evidence_status": spec["status"]}
