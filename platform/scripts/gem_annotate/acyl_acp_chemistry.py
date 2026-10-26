"""User-authorized 2026-10-20 acyl-ACP chemistry correction (opt-in).

The mtFAS and acyl-ACP export block (R1394-R1405, R2216-R2241) was copied from BiGG iMM904 with
placeholder formulas: "ACP", read as elements; "C141"/"C161", read as 141 or 161 carbons; formulas
with "*"; or none. Its stoichiometry was also written for ionized species.

This step:
- fills 32 acyl-ACP, ACP and fatty-acid formulas;
- removes the one surplus H+ in 19 reactions;
- makes the 7 reversible cytosolic acyl-ACP hydrolases forward-only (thioester hydrolysis is
  strongly exergonic).

Transport directions and reaction membership are unchanged. The decision and its evidence are data
in ``acyl_acp_chemistry.json``.

It requires the lipoylation module (ACP1 = HSR). Every species, reaction and residual is checked
before any edit: the model must be exactly in the before or the after state, or the step fails
closed. A rerun, also on a saved and reloaded model, is a no-op.
"""

import html
import json

from .config import MODEL

SPEC_PATH = MODEL.curation_file("acyl_acp_chemistry.json")
CHANGE_ID = "ACYL-ACP-CHEMISTRY-20261020"
NOTE_KEY = "acyl_acp_chemistry"
TOL = 1e-9


def load_spec():
    return json.loads(SPEC_PATH.read_text())


def _state(reaction):
    return {"stoichiometry": {m.id: float(c) for m, c in reaction.metabolites.items()},
            "bounds": [float(b) for b in reaction.bounds], "gpr": reaction.gene_reaction_rule}


def _imbalance(reaction):
    try:
        balance = reaction.check_mass_balance()
    except ValueError:
        return None
    return {k: v for k, v in balance.items() if abs(v) > TOL}


def same_imbalance(actual, expected):
    if actual is None or expected is None:
        return actual is None and expected is None
    return set(actual) == set(expected) and all(abs(actual[k] - expected[k]) <= TOL for k in expected)


def _note_differs(obj, value):
    return NOTE_KEY not in obj.notes or html.unescape(str(obj.notes[NOTE_KEY])) != html.unescape(value)


def _verify(model, spec, applied):
    state = "after" if applied else "before"
    required = spec["requires"]
    if required["reaction"] not in model.reactions:
        raise ValueError("The acyl-ACP chemistry step requires the lipoylation module; no edits applied")
    for mid, formula in required["species_formula"].items():
        if mid not in model.metabolites or model.metabolites.get_by_id(mid).formula != formula:
            raise ValueError(f"{mid}: required formula {formula} missing (lipoylation module); no edits applied")
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid}: species missing; no edits applied")
        m = model.metabolites.get_by_id(mid)
        if (m.formula or None, m.charge, m.compartment) != (row["formula"], row["charge"], row["compartment"]):
            raise ValueError(f"{mid}: species differs from the curated record; no edits applied")
    for mid, row in spec["formula_fills"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid}: species missing; no edits applied")
        m = model.metabolites.get_by_id(mid)
        if (m.formula or None, m.charge) != (row[state]["formula"], row[state]["charge"]):
            raise ValueError(f"{mid}: formula differs from the curated {state} state; no edits applied")
        if sorted(r.id for r in m.reactions) != spec["fill_users"][mid]:
            raise ValueError(f"{mid}: reactions using it differ from the curated record; no edits applied")
    reactions = spec[f"reactions_{state}"]
    residuals = spec[f"residuals_{state}"]
    for rid in spec["block_reactions"]:
        if rid not in model.reactions:
            raise ValueError(f"{rid}: reaction missing; no edits applied")
        reaction = model.reactions.get_by_id(rid)
        if _state(reaction) != reactions[rid]:
            raise ValueError(f"{rid}: definition differs from the curated {state} state; no edits applied")
        if not same_imbalance(_imbalance(reaction), residuals[rid]):
            raise ValueError(f"{rid}: residual differs from the curated {state} state; no edits applied")
        if applied and rid in spec["changed_reactions"] and _note_differs(reaction, spec["change_id"]):
            raise ValueError(f"{rid}: acyl-ACP chemistry note missing or conflicting")
    if applied:
        for mid in spec["formula_fills"]:
            if _note_differs(model.metabolites.get_by_id(mid), spec["change_id"]):
                raise ValueError(f"{mid}: acyl-ACP chemistry note missing or conflicting")


def apply_acyl_acp_chemistry(model, spec=None):
    spec = load_spec() if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("change_id") != CHANGE_ID
            or spec.get("status") != "user_authorized_candidate" or len(spec["formula_fills"]) != 32
            or len(spec["changed_reactions"]) != 19 or len(spec["made_irreversible"]) != 7):
        raise ValueError("Unexpected acyl-ACP chemistry scope")
    fills = spec["formula_fills"]
    applied = all((model.metabolites.get_by_id(mid).formula if mid in model.metabolites else None)
                  == row["after"]["formula"] for mid, row in fills.items())
    _verify(model, spec, applied)
    record = {"change_id": spec["change_id"], "formula_fills": len(fills),
              "changed_reactions": spec["changed_reactions"], "made_irreversible": spec["made_irreversible"]}
    if applied:
        return {**record, "status": "already_correct"}

    for mid, row in fills.items():
        metabolite = model.metabolites.get_by_id(mid)
        metabolite.formula, metabolite.charge = row["after"]["formula"], row["after"]["charge"]
        metabolite.notes = {**metabolite.notes, NOTE_KEY: spec["change_id"]}
    for rid in spec["changed_reactions"]:
        reaction, target = model.reactions.get_by_id(rid), spec["reactions_after"][rid]
        current = _state(reaction)["stoichiometry"]
        delta = {mid: target["stoichiometry"].get(mid, 0.0) - current.get(mid, 0.0)
                 for mid in sorted(set(current) | set(target["stoichiometry"]))}
        reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in delta.items() if c}, combine=True)
        reaction.bounds = tuple(target["bounds"])
        reaction.notes = {**reaction.notes, NOTE_KEY: spec["change_id"]}
    try:
        _verify(model, spec, applied=True)
    except ValueError as exc:
        raise RuntimeError(f"Acyl-ACP chemistry did not reach the curated after state: {exc}") from exc
    return {**record, "status": "applied"}
