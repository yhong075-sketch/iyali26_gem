"""Dump the baseline state of the trehalose targets to before.json (read-only).

Targets: the trehalose term m294[C_cy] of biomass_C (option C) and the R143 rule (option D, separate). Context, not
to be changed: every reaction that uses trehalose (m294, m1167, m1380) or trehalose 6-phosphate (m297), every
reaction of the trehalose genes, the trehalose terms of the other (non-objective) biomass reactions, and the objective.

Usage (from the repo root): python model/reports/trehalose_biomass_20261030/scripts/before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
from html import unescape
from pathlib import Path

from cobra.io import read_sbml_model

TREHALOSE = "m294[C_cy]"
SPECIES = ("m294[C_cy]", "m297[C_cy]", "m1167[C_ex]", "m1380[C_va]")
GENES = ("YALI1E17616g", "YALI1D17620g", "YALI1E36568g", "YALI1D19108g")
TREHALOSE_MW = 342.2965  # C12H22O11, g/mol, from the model formula and standard atomic weights

model_path, out = Path(sys.argv[1]), Path(sys.argv[2])
model = read_sbml_model(model_path)  # Path object: cobra rejects long str paths on macOS


def reaction_state(r, full=True):
    state = {"id": r.id, "name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
             "compartments": sorted(r.compartments)}
    if full:
        state.update({
            "equation": r.build_reaction_string(),
            "equation_names": r.build_reaction_string(use_metabolite_names=True),
            "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
            "species": {m.id: {k: getattr(m, k) for k in ("name", "formula", "charge", "compartment")}
                        for m in sorted(r.metabolites, key=lambda m: m.id)},
            "annotation": r.annotation,
            "notes": {k: unescape(v) if isinstance(v, str) else v for k, v in r.notes.items()},
        })
    return state


objective = [r.id for r in model.reactions if r.objective_coefficient != 0]
biomass = model.reactions.get_by_id("biomass_C")
treh = model.metabolites.get_by_id(TREHALOSE)
coef = biomass.metabolites[treh]

context = set()
for s in SPECIES:
    context |= {r.id for r in model.metabolites.get_by_id(s).reactions}
for g in GENES:
    context |= {r.id for r in model.genes.get_by_id(g).reactions}
pseudo = {rid for rid in context if rid == "biomass_C" or "biomass" in model.reactions.get_by_id(rid).name.lower()
          or rid in ("xBIOMASS", "newBiom")}


def producers_consumers(met_id):
    met = model.metabolites.get_by_id(met_id)
    prod, cons = [], []
    for r in met.reactions:
        c = r.metabolites[met]
        lo, hi = r.bounds
        if (c > 0 and hi > 0) or (c < 0 and lo < 0):
            prod.append(r.id)
        if (c < 0 and hi > 0) or (c > 0 and lo < 0):
            cons.append(r.id)
    return {"can_produce": sorted(prod), "can_consume": sorted(cons)}


state = {
    "baseline_model": {"path": str(model_path), "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest()},
    "model_counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "objective": objective,
    "targets": {
        "option_C_biomass_term": {"reaction": "biomass_C", "species": TREHALOSE, "coefficient": coef,
                                  "mass_mg_per_gDW": round(-coef * TREHALOSE_MW, 4),
                                  "nmol_per_mg_DW": round(-coef * 1000, 4)},
        "option_D_rule": {"reaction": "R143", "gpr": model.reactions.get_by_id("R143").gene_reaction_rule},
    },
    "directional_producers_consumers": {s: producers_consumers(s) for s in SPECIES},
    "trehalose_in_other_biomass_reactions": {
        rid: {"name": model.reactions.get_by_id(rid).name, "bounds": list(model.reactions.get_by_id(rid).bounds),
              "coefficient": model.reactions.get_by_id(rid).metabolites[treh]}
        for rid in sorted(pseudo)},
    "biomass_C_notes_on_trehalose": {k: v for k, v in biomass.notes.items() if "trehal" in (k + str(v)).lower()},
    "reactions": {rid: reaction_state(model.reactions.get_by_id(rid), full=rid not in pseudo)
                  for rid in sorted(context)},
    "genes": {g: {"name": model.genes.get_by_id(g).name, "annotation": model.genes.get_by_id(g).annotation,
                  "notes": model.genes.get_by_id(g).notes,
                  "reactions": sorted(x.id for x in model.genes.get_by_id(g).reactions)}
              for g in GENES},
}
out.write_text(json.dumps(state, indent=1, sort_keys=True) + "\n")
print("objective", objective, "| biomass_C", TREHALOSE, coef, state["targets"]["option_C_biomass_term"])
for s, pc in state["directional_producers_consumers"].items():
    print(s, pc)
for rid, r in state["trehalose_in_other_biomass_reactions"].items():
    print("pseudo", rid, r)
for rid, r in state["reactions"].items():
    print(rid, "|", r["name"], "|", r["compartments"], "|", r["bounds"], "|", r["gpr"])
for g, d in state["genes"].items():
    print(g, d["reactions"], d["annotation"].get("refseq"), d["annotation"].get("uniprot"))
