"""Diagnostic LPs for the heme a + iron biomass change (read-only; models are loaded in memory).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/heme_iron_biomass_20261022/scripts/solves.py \
    BASELINE.xml NEW.xml MEDIUM.csv PROFILE.json OUT.json
BASELINE is cf0825b8, NEW adds the two biomass_C terms. Condition: Task A medium + PO1f, GLPK.
"""
import json
import math
import sys
import warnings

from cobra.flux_analysis import pfba

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

WATCH = ["R1189", "R2193", "R863", "R302", "R384", "R383", "R674", "R1709", "R727", "R1070"]
CLOSE = ["R1189", "R2193", "R302", "R383", "R1709"]
TERMS = {"m1635[C_cy]": -1e-6, "m948[C_cy]": -3.04e-5}

baseline_path, new_path, medium, profile, out = sys.argv[1:6]
lp_count = 0


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    ctx.model.solver = "glpk"
    return ctx


def solve(model):
    global lp_count
    lp_count += 1
    value = model.slim_optimize(error_value=float("nan"))
    return {"status": model.solver.status, "growth": None if math.isnan(value) else value}


def parsimonious(model):
    global lp_count
    lp_count += 2
    solution = pfba(model)
    return {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
            "fluxes": {rid: float(solution.fluxes[rid]) for rid in WATCH}}


result = {"condition": "Task A medium + PO1f, GLPK", "iron_cap": None}
for name, path in (("baseline", baseline_path), ("new", new_path)):
    ctx = context(path)
    model = ctx.model
    biomass = {m.id: c for m, c in model.reactions.biomass_C.metabolites.items()}
    result[name] = {"model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
                    "strain_profile_sha256": ctx.strain_profile_sha256,
                    "biomass_terms": {mid: biomass.get(mid, 0.0) for mid in TERMS},
                    "r1189_bounds": list(model.reactions.R1189.bounds),
                    "wt": solve(model), "pfba": parsimonious(model)}
    if name == "new":
        result["iron_cap"] = -model.reactions.R1189.lower_bound
        result[name]["closed"] = {}
        for rid in CLOSE:
            with model:
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                result[name]["closed"][rid] = solve(model)
result["wt_change"] = result["new"]["wt"]["growth"] - result["baseline"]["wt"]["growth"]
result["iron_use_fraction_of_cap"] = -result["new"]["pfba"]["fluxes"]["R1189"] / result["iron_cap"]
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for name in ("baseline", "new"):
    r = result[name]
    print(name, r["model_sha256"][:8], "WT", r["wt"], "pFBA", r["pfba"]["growth"], r["pfba"]["fluxes"])
print("closed", result["new"]["closed"])
print("WT change", result["wt_change"], "iron use / cap", result["iron_use_fraction_of_cap"], "LPs", lp_count)
