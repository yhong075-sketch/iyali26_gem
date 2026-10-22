"""Record the state of every target of the heme a + iron biomass change (read-only).

Usage: python before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

TARGETS = ["biomass_C"]
CONTEXT = ["R1709", "R1710", "R1189", "R2193", "R863", "R302", "R384", "R383", "R674", "R727"]
SPECIES = ["m1635[C_cy]", "m948[C_cy]", "m609[C_mi]", "m525[C_mi]", "m608[C_mi]", "m523[C_mi]", "m1273[C_ex]"]

path, out = sys.argv[1], sys.argv[2]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(path)


def reaction(rid):
    r = model.reactions.get_by_id(rid)
    return {"name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
            "annotation": r.annotation, "notes": r.notes}


state = {
    "model": path.split("/")[-1],
    "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
    "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "objective": [r.id for r in model.reactions if r.objective_coefficient],
    "targets": {rid: reaction(rid) for rid in TARGETS},
    "context": {rid: reaction(rid) if rid in model.reactions else None for rid in CONTEXT},
    "species": {mid: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
                      "reactions": sorted(r.id for r in m.reactions)}
                for mid in SPECIES for m in [model.metabolites.get_by_id(mid)]},
    "fe_containing_in_biomass_C": sorted(m.id for m in model.reactions.biomass_C.metabolites
                                         if m.formula and "Fe" in m.formula or m.id == "m948[C_cy]"),
}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(json.dumps({"sha256": state["sha256"], "counts": state["counts"], "objective": state["objective"],
                  "fe_containing_in_biomass_C": state["fe_containing_in_biomass_C"],
                  "species": {k: v["reactions"] for k, v in state["species"].items()}}, indent=1))
