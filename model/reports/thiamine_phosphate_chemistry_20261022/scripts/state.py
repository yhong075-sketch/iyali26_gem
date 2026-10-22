"""State of every target of the thiamine phosphate chemistry change, with per-reaction balance (read-only).

Usage: python state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

SPECIES = ["m977[C_cy]", "m973[C_ex]", "m975[C_cy]", "m971[C_ex]", "m1116[C_mi]", "m976[C_cy]", "m643[C_cy]"]

path, out = sys.argv[1], sys.argv[2]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(path)


def balance(reaction):
    if reaction.boundary:
        return "boundary"
    try:
        residual = {k: v for k, v in reaction.check_mass_balance().items() if abs(v) > 1e-9}
    except ValueError as error:
        return f"not computable: {error}"
    return residual or "balanced"


reactions = sorted({r.id for mid in SPECIES for r in model.metabolites.get_by_id(mid).reactions})
state = {
    "model": path.split("/")[-1],
    "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
    "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "species": {mid: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
                      "annotation": m.annotation, "notes": m.notes, "reactions": sorted(r.id for r in m.reactions)}
                for mid in SPECIES for m in [model.metabolites.get_by_id(mid)]},
    "reactions": {rid: {"name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
                        "stoichiometry": ({m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)}
                                          if rid != "biomass_C" else {"m975[C_cy]": r.get_coefficient("m975[C_cy]")}),
                        "annotation": r.annotation, "notes": r.notes, "balance": balance(r)}
                  for rid in reactions for r in [model.reactions.get_by_id(rid)]},
}
# Partners: every other species on these reactions, whose chemistry the balance review assumed.
state["partners"] = {m.id: [m.formula, m.charge] for rid in reactions if rid != "biomass_C"
                     for m in model.reactions.get_by_id(rid).metabolites if m.id not in SPECIES}
state["partners"] = dict(sorted(state["partners"].items()))
state["balanced_count"] = sum(v["balance"] == "balanced" for v in state["reactions"].values())
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(json.dumps({"sha256": state["sha256"], "counts": state["counts"], "balanced_count": state["balanced_count"],
                  "balance": {k: v["balance"] for k, v in state["reactions"].items()}}, indent=1))
