"""Dump the baseline state of the ceramide-pool targets to before.json (read-only).

Targets: biomass_C (the objective; one new pool term) and the 10 ER ceramide species that would feed the pool.
Also records that the pool metabolite and the xPOOL_CER_* reactions do not exist yet, and the producers/consumers
of each ceramide (context, not changed).

Usage (from the repo root): python model/reports/ceramide_pool_biomass_20261101/scripts/before_state.py MODEL.xml OUT.json
"""
import hashlib
import json
import sys
from html import unescape
from pathlib import Path

from cobra.io import read_sbml_model

CERAMIDES = ("m386", "m388", "m390", "m391", "m392", "m393", "m397", "m398", "m399", "m400")
POOL = "sphingolipid_ceramide_pool[C_er]"

model_path, out = Path(sys.argv[1]), Path(sys.argv[2])
model = read_sbml_model(model_path)  # Path object: cobra rejects long str paths on macOS
biomass = model.reactions.get_by_id("biomass_C")


def species(met_id):
    met = model.metabolites.get_by_id(met_id)
    reactions = []
    for r in sorted(met.reactions, key=lambda r: r.id):
        reactions.append({"id": r.id, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
                          "coefficient": r.metabolites[met], "equation": r.build_reaction_string()})
    return {"id": met.id, "name": met.name, "formula": met.formula, "charge": met.charge,
            "compartment": met.compartment, "reactions": reactions}


state = {
    "model": {"path": str(model_path), "sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
              "reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "objective": {r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient},
    "objective_direction": model.objective.direction,
    "biomass_C": {
        "bounds": list(biomass.bounds),
        "stoichiometry": {m.id: c for m, c in sorted(biomass.metabolites.items(), key=lambda x: x[0].id)},
        "note_keys": sorted(biomass.notes),
        "notes": {k: unescape(v) if isinstance(v, str) else v for k, v in sorted(biomass.notes.items())},
    },
    "ceramides": [species(f"{mid}[C_er]") for mid in CERAMIDES],
    "pool_present": POOL in model.metabolites,
    "xPOOL_CER_reactions_present": sorted(r.id for r in model.reactions if r.id.startswith("xPOOL_CER_")),
}
out.write_text(json.dumps(state, indent=1) + "\n")
print(state["model"], "pool present:", state["pool_present"], state["xPOOL_CER_reactions_present"])
for c in state["ceramides"]:
    print(c["id"], c["name"], c["formula"], len(c["reactions"]), "reactions")
