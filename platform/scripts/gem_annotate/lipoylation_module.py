"""User-authorized 2026-10-20 lipoylation module (opt-in).

Adds the mitochondrial lipoylation route the model lacks:
- a LIP2-type octanoyl transfer from octanoyl-ACP to the glycine-cleavage H-protein, which returns
  ACP1;
- lipoyl synthase (LIP5, with lumped sulfur and electron donors);
- a LIP3-type transfer of lipoyl to the E2 lipoyl domains;
- a provisional growth-dilution requirement for both lipoylated proteins.

It also fills the ACP1 formula and corrects the lipoyl-H-protein formula. The corrected formula
makes the H+ in R362 and R364 unnecessary, so it is removed there. The decision and its evidence are
data in ``lipoylation_module.json``.

Every species, reaction, gene and biomass term is checked before any edit. The model must be
exactly in the before state or the after state, or the step fails closed. A rerun, also on a saved
and reloaded model, is a no-op.
"""

import html
import json

from cobra import Metabolite, Reaction

from .config import MODEL

SPEC_PATH = MODEL.curation_file("lipoylation_module.json")
CHANGE_ID = "LIPOYLATION-MODULE-20261020"
TOL = 1e-9


def _stoichiometry(reaction):
    return {m.id: float(c) for m, c in reaction.metabolites.items()}


def _matches(reaction, state):
    return (_stoichiometry(reaction), [float(b) for b in reaction.bounds],
            reaction.gene_reaction_rule.strip()) == (state["stoichiometry"], state["bounds"], state["gpr"])


def _imbalance(reaction):
    try:
        balance = reaction.check_mass_balance()
    except ValueError:
        return None
    return {k: v for k, v in balance.items() if abs(v) > TOL}


def _same_imbalance(actual, expected):
    return (actual is not None and set(actual) == set(expected)
            and all(abs(actual[k] - expected[k]) <= TOL for k in expected))


def _species_state(metabolite):
    return {"formula": metabolite.formula or None, "charge": metabolite.charge, "compartment": metabolite.compartment}


def _reaction_notes(spec, row):
    notes = {"lipoylation_module": spec["change_id"], "lipoylation_evidence": row["evidence"],
             "lipoylation_evidence_tier": row["evidence_tier"], "lipoylation_sources": "; ".join(row["sources"])}
    return {key: value for key, value in notes.items() if value}  # empty notes do not survive SBML


def _reaction_annotation(row):
    return {key: value for key, value in (("ec-code", row["ec"]), ("rhea", row["rhea"])) if value}


def _note_differs(obj, key, value):
    # SBML reload returns XML entities (e.g. &apos;) as literal text; compare unescaped.
    return key not in obj.notes or html.unescape(str(obj.notes[key])) != html.unescape(value)


def _expected_users(spec, mid, applied):
    users = set(spec["formula_fill_users"][mid])
    if applied:
        users |= {rid for rid, row in spec["new_reactions"].items() if mid in row["stoichiometry"]}
        if mid in spec["biomass_terms"]:
            users.add(spec["biomass_reaction"])
    return users


def _verify(model, spec, applied):
    """Raise ValueError unless the model is exactly in the before (applied=False) or after state."""
    state = "after" if applied else "before"
    fills = spec["formula_fills"]
    for mid, row in spec["species"].items():
        if mid not in model.metabolites:
            raise ValueError(f"{mid}: species missing; no edits applied")
        metabolite = model.metabolites.get_by_id(mid)
        fill = fills.get(mid, {}).get(state, {})
        expected = {**row, **{k: v for k, v in fill.items() if k in ("formula", "charge")}}
        if _species_state(metabolite) != expected or ("name" in fill and metabolite.name != fill["name"]):
            raise ValueError(f"{mid}: species differs from the curated {state} state; no edits applied")
    for mid in fills:
        users = {r.id for r in model.metabolites.get_by_id(mid).reactions}
        if users != _expected_users(spec, mid, applied):
            raise ValueError(f"{mid}: reactions using it differ from the curated record; no edits applied")
    for rid, row in spec["edited_reactions"].items():
        if rid not in model.reactions or not _matches(model.reactions.get_by_id(rid), row[state]):
            raise ValueError(f"{rid}: definition differs from the curated {state} state; no edits applied")
    for mid, row in spec["new_metabolites"].items():
        if (mid in model.metabolites) != applied:
            raise ValueError(f"{mid}: presence does not match the {state} state; no edits applied")
        if applied:
            m = model.metabolites.get_by_id(mid)
            if (m.name, m.formula, m.charge, m.compartment) != (row["name"], row["formula"], row["charge"], row["compartment"]):
                raise ValueError(f"{mid}: new species differs from the curated record")
    for rid, row in spec["new_reactions"].items():
        if (rid in model.reactions) != applied:
            raise ValueError(f"{rid}: presence does not match the {state} state; no edits applied")
        if applied:
            reaction = model.reactions.get_by_id(rid)
            if reaction.name != row["name"] or not _matches(reaction, row):
                raise ValueError(f"{rid}: new reaction differs from the curated record")
            if any(_note_differs(reaction, k, v) for k, v in _reaction_notes(spec, row).items()):
                raise ValueError(f"{rid}: lipoylation notes missing or conflicting")
            if any(reaction.annotation.get(k) not in (v, [v]) for k, v in _reaction_annotation(row).items()):
                raise ValueError(f"{rid}: EC/Rhea annotation missing or conflicting")
    for gid, row in spec["new_genes"].items():
        if (gid in model.genes) != applied:
            raise ValueError(f"{gid}: presence does not match the {state} state; no edits applied")
        if applied and model.genes.get_by_id(gid).annotation.get("refseq") not in (row["refseq"], [row["refseq"]]):
            raise ValueError(f"{gid}: gene identity differs from the curated record")
    for gid in spec["context_genes"]:
        if gid not in model.genes:
            raise ValueError(f"{gid}: context gene missing; no edits applied")
    if spec["biomass_reaction"] not in model.reactions:
        raise ValueError("biomass reaction missing; no edits applied")
    biomass = model.reactions.get_by_id(spec["biomass_reaction"])
    for mid, coefficient in spec["biomass_terms"].items():
        current = (biomass.metabolites.get(model.metabolites.get_by_id(mid), 0.0)
                   if mid in model.metabolites else 0.0)
        if abs(current - (coefficient if applied else 0.0)) > TOL:
            raise ValueError(f"biomass term for {mid} differs from the curated {state} state; no edits applied")
    if not applied:
        for rid, residual in spec["acp1_user_residuals_before_fill"].items():
            if not _same_imbalance(_imbalance(model.reactions.get_by_id(rid)), residual):
                raise ValueError(f"{rid}: residual differs from the recorded before state; no edits applied")
    if applied:
        for rid in spec["balanced_after"]:
            if _imbalance(model.reactions.get_by_id(rid)) != {}:
                raise ValueError(f"{rid}: not balanced after the lipoylation module")
        for rid, residual in spec["preexisting_residuals_after_fill"].items():
            reaction = model.reactions.get_by_id(rid)
            if "acyl_acp_chemistry" in reaction.notes:
                # The later ACYL-ACP-CHEMISTRY-20261020 step corrected this reaction; use its record.
                from .acyl_acp_chemistry import CHANGE_ID as ACYL_CHANGE_ID, load_spec
                if html.unescape(str(reaction.notes["acyl_acp_chemistry"])) != ACYL_CHANGE_ID:
                    raise ValueError(f"{rid}: unexpected acyl-ACP chemistry note")
                residual = load_spec()["residuals_after"][rid]
            if not _same_imbalance(_imbalance(reaction), residual):
                raise ValueError(f"{rid}: residual differs from the recorded pre-existing imbalance")


def apply_lipoylation_module(model, spec=None):
    spec = json.loads(SPEC_PATH.read_text()) if spec is None else spec
    if (spec.get("schema_version") != 1 or spec.get("change_id") != CHANGE_ID
            or spec.get("status") != "user_authorized_candidate"
            or sorted(spec["new_reactions"]) != ["R2286", "R2287", "R2288", "R2289", "R2290"]
            or sorted(spec["formula_fills"]) != ["m179[C_mi]", "m586[C_mi]"]
            or sorted(spec["edited_reactions"]) != ["R362", "R364"]):
        raise ValueError("Unexpected lipoylation module scope")
    present = [rid for rid in spec["new_reactions"] if rid in model.reactions]
    if present and len(present) != len(spec["new_reactions"]):
        raise ValueError(f"Partial lipoylation module present ({present}); no edits applied")
    applied = bool(present)
    _verify(model, spec, applied)
    record = {"change_id": spec["change_id"], "alpha_mmol_per_gDW": spec["alpha_mmol_per_gDW"],
              "new_reactions": sorted(spec["new_reactions"]), "new_genes": sorted(spec["new_genes"])}
    if applied:
        return {**record, "status": "already_correct"}

    for mid, row in spec["formula_fills"].items():
        metabolite = model.metabolites.get_by_id(mid)
        metabolite.formula, metabolite.charge = row["after"]["formula"], row["after"]["charge"]
        if "name" in row["after"]:
            metabolite.name = row["after"]["name"]
        metabolite.notes = {**metabolite.notes, "lipoylation_module_formula": f"{spec['change_id']}: {row['reason']}"}
    for rid, row in spec["edited_reactions"].items():
        reaction = model.reactions.get_by_id(rid)
        current, target = _stoichiometry(reaction), row["after"]["stoichiometry"]
        delta = {mid: target.get(mid, 0.0) - current.get(mid, 0.0) for mid in sorted(set(current) | set(target))}
        reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in delta.items() if c}, combine=True)
        reaction.notes = {**reaction.notes, "lipoylation_module_edit": f"{spec['change_id']}: {row['reason']}"}
    model.add_metabolites([Metabolite(mid, name=row["name"], formula=row["formula"], charge=row["charge"],
                                      compartment=row["compartment"])
                           for mid, row in spec["new_metabolites"].items()])
    for rid, row in spec["new_reactions"].items():
        reaction = Reaction(rid, name=row["name"], lower_bound=row["bounds"][0], upper_bound=row["bounds"][1])
        model.add_reactions([reaction])
        reaction.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in row["stoichiometry"].items()})
        reaction.gene_reaction_rule = row["gpr"]
        reaction.annotation = _reaction_annotation(row)
        reaction.notes = _reaction_notes(spec, row)
    for gid, row in spec["new_genes"].items():
        gene = model.genes.get_by_id(gid)
        gene.annotation = {**gene.annotation, "refseq": row["refseq"]}
        gene.notes = {**gene.notes, "lipoylation_gene_function": row["function"],
                      "lipoylation_evidence_tier": row["tier"], "lipoylation_module": spec["change_id"]}
    biomass = model.reactions.get_by_id(spec["biomass_reaction"])
    biomass.add_metabolites({model.metabolites.get_by_id(mid): c for mid, c in spec["biomass_terms"].items()},
                            combine=True)
    biomass.notes = {**biomass.notes, "lipoylation_growth_dilution": (
        f"{spec['change_id']}: alpha {spec['alpha_mmol_per_gDW']} mmol/gDW each for lipoyl-H-protein and lipoyl-E2, "
        f"apo forms returned; {spec['alpha_source']}")}
    try:
        _verify(model, spec, applied=True)
    except ValueError as exc:
        raise RuntimeError(f"Lipoylation module did not reach the curated after state: {exc}") from exc
    return {**record, "status": "applied"}
