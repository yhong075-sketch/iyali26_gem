"""Record the baseline state of the R1726 / R1725 targets and their route (read-only, no solves).

Usage (from platform/): python ../model/reports/kennedy_closures_20261027/scripts/before_state.py MODEL.xml OUT.json
Per target: name, bounds, reversibility, GPR, stoichiometry, species formula / charge / compartment, notes,
annotations. Per route reaction: bounds and GPR. The model file's SHA-256 is recorded.
"""
import hashlib
import json
import sys
from pathlib import Path

from cobra.io import read_sbml_model

TARGETS = ("R1726", "R1725")
ROUTE = ("R212", "R213", "R1603", "R293", "R2080", "R635", "R1602", "R1042", "R1806", "R1842", "R1743", "R1776",
         "R1720", "SPHPL", "R629", "R630", "R631", "R1835", "R1845", "R1844")
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
         "gene_YALI1E31498g_in_model": "YALI1E31498g" in model.genes}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(state["sha256"], state["counts"], {rid: state["targets"][rid]["bounds"] for rid in TARGETS})
