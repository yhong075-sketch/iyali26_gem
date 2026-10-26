"""Before-state for combining the TPP and heme/iron + R2193 lines (read-only).

Usage: python before_state.py BASE.xml TPP_PARENT.xml HEME_PARENT.xml TPP_CURATION.json OUT.json
Records each model's SHA-256 and counts, biomass_C stoichiometry and note keys, and R2193; and checks that
the TPP parent's biomass_C equals the base plus exactly the TPP curation's added term.
"""
import hashlib
import json
import sys
import warnings

from cobra.io import read_sbml_model

base_path, tpp_path, heme_path, tpp_spec_path, out = sys.argv[1:6]


def state(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = read_sbml_model(path)
    biomass, r2193 = model.reactions.biomass_C, model.reactions.R2193
    return {"file": path.split("/")[-1], "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
            "counts": {"reactions": len(model.reactions), "metabolites": len(model.metabolites), "genes": len(model.genes)},
            "objective": [r.id for r in model.reactions if r.objective_coefficient],
            "biomass_C": {"bounds": list(biomass.bounds),
                          "stoichiometry": {m.id: c for m, c in sorted(biomass.metabolites.items(), key=lambda x: x[0].id)},
                          "note_keys": sorted(biomass.notes)},
            "R2193": {"gpr": r2193.gene_reaction_rule, "bounds": list(r2193.bounds), "note_keys": sorted(r2193.notes)}}


result = {"base": state(base_path), "tpp_parent": state(tpp_path), "heme_parent": state(heme_path)}
tpp_spec = json.load(open(tpp_spec_path))
base_b, tpp_b = result["base"]["biomass_C"]["stoichiometry"], result["tpp_parent"]["biomass_C"]["stoichiometry"]
added = {k: v for k, v in tpp_b.items() if k not in base_b}
result["tpp_parent_biomass_vs_base"] = {
    "added": added, "removed": sorted(set(base_b) - set(tpp_b)),
    "changed": {k: [base_b[k], tpp_b[k]] for k in set(base_b) & set(tpp_b) if base_b[k] != tpp_b[k]},
    "added_note_keys": sorted(set(result["tpp_parent"]["biomass_C"]["note_keys"]) - set(result["base"]["biomass_C"]["note_keys"])),
    "matches_tpp_curation_after_stoichiometry": tpp_b == tpp_spec["targets"]["biomass_C"]["after_stoichiometry"]}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: {"sha256": v["sha256"][:8], "counts": v["counts"], "R2193": v["R2193"]["gpr"]} for k, v in result.items() if "sha256" in v}, indent=1))
print(json.dumps(result["tpp_parent_biomass_vs_base"], indent=1))
