"""Record the state of every target of the GPI co-product change (read-only).

Usage: python before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

TARGETS = ["R373", "R375"]
CONTEXT = ["R372", "R374", "R376", "R_GPIGNT", "R_PEt_em_er", "R_DAGt_er_em", "R_ACt_er_cy", "R_DM_GPI"]
SPECIES = ["m598[C_er]", "m599[C_er]", "m600[C_er]", "m502[C_er]", "m601[C_er]", "m602[C_er]",
           "m503[C_er]", "m604[C_er]"]

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
    "targets": {rid: reaction(rid) for rid in TARGETS},
    "context": {rid: reaction(rid) if rid in model.reactions else None for rid in CONTEXT},
    "species": {mid: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
                      "reactions": sorted(r.id for r in m.reactions)}
                for mid in SPECIES for m in [model.metabolites.get_by_id(mid)]},
}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(json.dumps({"sha256": state["sha256"], "counts": state["counts"], "targets": state["targets"]}, indent=1))
