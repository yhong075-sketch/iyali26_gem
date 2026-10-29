"""Record R2171 and its species on the baseline before the change (static read, no LP).

Usage (from platform/): python ../model/reports/r2171_pi_plc_1to1_20261028/scripts/before_state.py MODEL.xml OUT.json
The model is read through a Path object (cobra rejects string paths of 260+ characters on macOS).
"""
import hashlib
import json
import sys
import warnings
from pathlib import Path

from cobra.io import read_sbml_model

path, out = Path(sys.argv[1]), sys.argv[2]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(path)
r = model.reactions.R2171
state = {"model": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
         "R2171": {"name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule, "subsystem": r.subsystem,
                   "stoichiometry": {m.id: c for m, c in r.metabolites.items()}, "annotation": dict(r.annotation),
                   "notes": dict(r.notes), "charge_balance": sum(m.charge * c for m, c in r.metabolites.items())},
         "species": {m.id: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
                            "kegg.compound": m.annotation.get("kegg.compound"),
                            "metanetx.chemical": m.annotation.get("metanetx.chemical"),
                            "n_reactions": len(m.reactions)} for m in r.metabolites}}
json.dump(state, open(out, "w"), indent=1)
print(state["sha256"][:8], state["R2171"]["stoichiometry"], "charge", state["R2171"]["charge_balance"], file=sys.stderr)
