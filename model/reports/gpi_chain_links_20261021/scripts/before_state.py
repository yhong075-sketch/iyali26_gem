"""Record the baseline state of every target of the GPI chain links change (read-only).

Usage: python before_state.py BASELINE.xml OUT.json
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

NEW_REACTIONS = ["R_GPIGNT", "R_PEt_em_er", "R_DAGt_er_em", "R_ACt_er_cy", "R_DM_GPI"]
METABOLITES = ["m403[C_cy]", "m1648[C_em]", "m10[C_cy]", "m11[C_cy]", "m596[C_er]", "m1633[C_em]",
               "m599[C_er]", "m604[C_er]", "m1649[C_em]", "m597[C_er]", "m202[C_cy]", "m603[C_er]"]
CHAIN = ["R371", "R113", "R372", "R373", "R374", "R375", "R376"]

path, out = sys.argv[1], sys.argv[2]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(path)
state = {
    "model": path.split("/")[-1],
    "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
    "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "new_reaction_ids_present": {rid: rid in model.reactions for rid in NEW_REACTIONS},
    "metabolites": {mid: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
                          "reactions": sorted(r.id for r in m.reactions)}
                    for mid in METABOLITES for m in [model.metabolites.get_by_id(mid)]},
    "chain_reactions": {rid: {"name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
                              "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)}}
                        for rid in CHAIN for r in [model.reactions.get_by_id(rid)]},
}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: state[k] for k in ("sha256", "counts", "new_reaction_ids_present")}, indent=1))
