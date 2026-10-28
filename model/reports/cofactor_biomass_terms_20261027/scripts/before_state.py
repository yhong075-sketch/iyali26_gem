"""Baseline state of every cofactor-biomass target in a model (read-only).

Usage (from the worktree root): python -I before_state.py MODEL.xml OUT.json
Records the SHA-256, biomass_C (stoichiometry, bounds, notes, objective), the six cofactor species
(name, formula, charge, compartment, annotation, reaction list) in every compartment copy, and the
route reactions of TASK.md section 4 (stoichiometry, bounds, GPR).
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

TERMS = {"m479[C_cy]": 6.34e-5, "m528[C_cy]": 1e-5, "m926[C_cy]": 9.9e-4,
         "m122[C_cy]": 2.65e-3, "m41[C_cy]": 5.7e-4, "m69[C_cy]": 1.9e-4}
ROUTES = {
    "folate": ["R377", "R274", "R14", "R277", "R278", "R273", "R749", "R271", "R272", "R849", "R1127", "R811", "R991"],
    "riboflavin_fad": ["R378", "R34", "R12", "R11", "R699", "R700", "R1316", "R1315", "R1314", "R697", "R698", "R306", "R847"],
    "coa": ["R946", "R616", "R1908", "R1909", "R647", "R648", "R615", "R1092", "R1291", "R265"],
    "nad_nadp": ["R1268", "R907", "R576", "R578", "R579", "R580", "R568", "R569", "R565", "R566", "R571", "R572"],
}

model_path, out = sys.argv[1:3]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(model_path)
biomass = model.reactions.get_by_id("biomass_C")
species = {}
for mid in TERMS:
    base = mid.split("[")[0]
    for m in model.metabolites:
        if m.id.split("[")[0] == base:
            species[m.id] = {"name": m.name, "formula": m.formula or None, "charge": m.charge,
                             "compartment": m.compartment, "annotation": m.annotation,
                             "reactions": sorted(r.id for r in m.reactions)}
state = {
    "model": model_path, "sha256": hashlib.sha256(open(model_path, "rb").read()).hexdigest(),
    "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
    "objective": {r.id: r.objective_coefficient for r in model.reactions if r.objective_coefficient},
    "biomass_C": {"stoichiometry": {m.id: c for m, c in sorted(biomass.metabolites.items(), key=lambda x: x[0].id)},
                  "bounds": list(biomass.bounds), "notes": biomass.notes, "metabolite_count": len(biomass.metabolites),
                  "terms_present": {mid: biomass.metabolites.get(model.metabolites.get_by_id(mid), 0.0) for mid in TERMS}},
    "planned_terms_mmol_per_gDW": TERMS,
    "species": species,
    "routes": {name: {rid: ({"name": r.name, "stoichiometry": {m.id: c for m, c in sorted(r.metabolites.items(), key=lambda x: x[0].id)},
                             "bounds": list(r.bounds), "gpr": r.gene_reaction_rule}
                            if rid in model.reactions else None)
                      for rid in rids for r in [model.reactions.get_by_id(rid) if rid in model.reactions else None]}
               for name, rids in ROUTES.items()},
}
with open(out, "w") as handle:
    json.dump(state, handle, indent=1, default=str)
    handle.write("\n")
print(state["sha256"], "biomass_C metabolites", state["biomass_C"]["metabolite_count"],
      "terms present", state["biomass_C"]["terms_present"], "species copies", len(species))
