"""Record the baseline state of the R2135 / R2145 targets and their route (read-only, no solves).

Usage (from platform/): python ../model/reports/r2135_acyl_coa_pool_20261027/scripts/before_state.py MODEL.xml OUT.json
Per target: name, bounds, reversibility, GPR, stoichiometry, species formula / charge / compartment, notes,
annotations. Per route reaction: bounds and GPR. The model file's SHA-256 is recorded.
"""
import hashlib
import json
import sys
from pathlib import Path

from cobra.io import read_sbml_model

TARGETS = ("R2135", "R2145")
ROUTE = ("R2134", "R2168", "R2144", "R2169", "R2170", "R1771", "xPOOL_AC_EM", "R1726", "R1725", "R1720")
model_path, out = sys.argv[1:3]
model = read_sbml_model(model_path)


def reaction_state(r):
    return {"name": r.name, "bounds": list(r.bounds), "reversibility": r.reversibility, "gpr": r.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda kv: kv[0].id)},
            "species": {m.id: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment}
                        for m in sorted(r.metabolites, key=lambda m: m.id)},
            "notes": dict(r.notes), "annotation": dict(r.annotation)}


state = {"model": model_path, "sha256": hashlib.sha256(Path(model_path).read_bytes()).hexdigest(),
         "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
         "targets": {rid: reaction_state(model.reactions.get_by_id(rid)) for rid in TARGETS},
         "route": {rid: {"name": model.reactions.get_by_id(rid).name, "bounds": list(model.reactions.get_by_id(rid).bounds),
                         "gpr": model.reactions.get_by_id(rid).gene_reaction_rule} for rid in ROUTE},
         "species": {mid: {"name": model.metabolites.get_by_id(mid).name, "formula": model.metabolites.get_by_id(mid).formula, "charge": model.metabolites.get_by_id(mid).charge, "compartment": model.metabolites.get_by_id(mid).compartment, "annotation": dict(model.metabolites.get_by_id(mid).annotation), "reactions": sorted(r.id for r in model.metabolites.get_by_id(mid).reactions)} for mid in ("m68[C_cy]", "m567[C_em]", "m1984[C_er]", "m1993[C_cy]", "m1994[C_cy]")}}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(state["sha256"], state["counts"], {rid: state["targets"][rid]["bounds"] for rid in TARGETS})
