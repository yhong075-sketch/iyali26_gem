"""Record the baseline state of the dolichol-route targets (read-only, no solves).

Usage (from platform/): python ../model/reports/dolichol_route_20261028/scripts/before_state.py MODEL.xml OUT.json
Per target reaction (R2178 to be closed, R287 dolichol kinase, the UM22 fragment R2179 / R2180 / R2177 /
R2181 / R2182, biomass_C): name, bounds, reversibility, GPR, stoichiometry, species identity, notes,
annotations, group membership. Per species of the dolichol family and the cofactors of the new step: name,
formula, charge, compartment and the reactions it takes part in. The model file's SHA-256 is recorded.
"""
import hashlib
import json
import sys
from pathlib import Path

from cobra.io import read_sbml_model

TARGETS = ("R2178", "R287", "R2179", "R2180", "R2177", "R2181", "R2182", "biomass_C")
SPECIES = ("m2019[C_cy]", "m2020[C_cy]", "m499[C_cy]", "m501[C_cy]", "m503[C_er]", "m2018[C_cy]", "m2021[C_er]",
           "m32[C_cy]", "m122[C_cy]", "m123[C_cy]", "m40[C_cy]", "m41[C_cy]", "m10[C_cy]", "m203[C_cy]",
           "m109[C_cy]", "m82[C_cy]")
NEW_ID = "R_DOLSYN"
model_path, out = sys.argv[1:3]
model = read_sbml_model(model_path)


def reaction_state(r):
    return {"name": r.name, "bounds": list(r.bounds), "reversibility": r.reversibility, "gpr": r.gene_reaction_rule,
            "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda kv: kv[0].id)},
            "species": {m.id: {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment}
                        for m in sorted(r.metabolites, key=lambda m: m.id)},
            "notes": dict(r.notes), "annotation": dict(r.annotation),
            "groups": sorted(g.name for g in model.groups if r in g.members),
            "objective_coefficient": r.objective_coefficient}


def species_state(mid):
    m = model.metabolites.get_by_id(mid)
    return {"name": m.name, "formula": m.formula, "charge": m.charge, "compartment": m.compartment,
            "reactions": sorted(r.id for r in m.reactions)}


state = {"model": model_path, "sha256": hashlib.sha256(Path(model_path).read_bytes()).hexdigest(),
         "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes),
                    "groups": len(model.groups)},
         "objective": {r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient},
         "new_reaction_id_absent": NEW_ID not in model.reactions,
         "targets": {rid: reaction_state(model.reactions.get_by_id(rid)) for rid in TARGETS},
         "species": {mid: species_state(mid) for mid in SPECIES},
         "genes_in_model": {g: g in model.genes for g in ("YALI1E22718g", "YALI1D26639g", "YALI1C26036g", "YALI1D21543g",
                                                          "YALI1E09357g", "YALI1F26726g")}}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1)
    handle.write("\n")
print(state["sha256"], state["counts"], {rid: state["targets"][rid]["bounds"] for rid in TARGETS if rid != "biomass_C"},
      "biomass species", len(state["targets"]["biomass_C"]["stoichiometry"]), "R_DOLSYN absent", state["new_reaction_id_absent"])
