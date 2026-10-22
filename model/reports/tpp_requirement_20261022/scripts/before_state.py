"""Record the state of every target of the TPP requirement change (read-only).

Usage: python before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

TARGETS = ["biomass_C", "R1342"]
CONTEXT = ["R752", "R753", "R933", "R934", "R1337", "R1338", "R754", "R408", "R413", "R645", "R970",
           "R1343", "R1344", "R1340"]
SPECIES = ["m975[C_cy]", "m976[C_cy]", "m977[C_cy]", "m974[C_cy]", "m1116[C_mi]"]
GENE = "YALI1E25443g"

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
    "objective": str(model.objective.expression),
    "targets": {rid: reaction(rid) for rid in TARGETS},
    "context": {rid: reaction(rid) for rid in CONTEXT},
    "species": {mid: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
                      "annotation": m.annotation, "reactions": sorted(r.id for r in m.reactions)}
                for mid in SPECIES for m in [model.metabolites.get_by_id(mid)]},
    "gene": {GENE: sorted(r.id for r in model.genes.get_by_id(GENE).reactions)},
}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(json.dumps({"sha256": state["sha256"], "counts": state["counts"],
                  "R1342": state["targets"]["R1342"]["bounds"],
                  "biomass_C_has_tpp": "m975[C_cy]" in state["targets"]["biomass_C"]["stoichiometry"],
                  "m975_reactions": state["species"]["m975[C_cy]"]["reactions"], "gene": state["gene"]}, indent=1))
