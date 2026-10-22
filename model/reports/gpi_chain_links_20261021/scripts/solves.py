"""Diagnostic LPs for the GPI chain links change (read-only; models are loaded in memory).

Usage (from platform/): python ../model/reports/gpi_chain_links_20261021/scripts/solves.py \
    BASELINE.xml NEW.xml RESEARCH_ROOT OUT.json
Condition: SD-Leu + PO1f, GLPK, objective unchanged.
"""
import json
import math
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis, pfba

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

CHAIN = ["R371", "R113", "R372", "R373", "R374", "R375", "R376"]
LINKS = ["R_GPIGNT", "R_PEt_em_er", "R_DAGt_er_em", "R_ACt_er_cy", "R_DM_GPI"]
FORCED = 1e-5

baseline_path, new_path, research, out = sys.argv[1:5]
lp_count = 0


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(
            model_path=path, media_path=f"{research}/state/media/sd_leu.csv",
            strain_profile_path=f"{research}/state/strain_profiles/po1f_sd_leu.json")
    ctx.model.solver = "glpk"
    return ctx


def solve(model):
    global lp_count
    lp_count += 1
    value = model.slim_optimize(error_value=float("nan"))
    return {"status": model.solver.status, "growth": None if math.isnan(value) else value}


def fva(model, ids):
    global lp_count
    lp_count += 2 * len(ids)
    table = flux_variability_analysis(model, ids, fraction_of_optimum=0.0, processes=1)
    return {rid: [float(table.loc[rid, "minimum"]), float(table.loc[rid, "maximum"])] for rid in ids}


result = {"condition": "SD-Leu + PO1f, GLPK", "forced_exit": FORCED}
base, new = context(baseline_path), context(new_path)
for name, ctx in (("baseline", base), ("new", new)):
    result[name] = {"model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
                    "strain_profile_sha256": ctx.strain_profile_sha256, "wt": solve(ctx.model),
                    "fva_fraction_0": fva(ctx.model, CHAIN + (LINKS if name == "new" else []))}

model = new.model
with model:
    model.reactions.R_DM_GPI.lower_bound = FORCED
    lp_count += 2
    solution = pfba(model)
    result["new"]["forced_exit_pfba"] = {
        "status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
        "fluxes": {rid: float(solution.fluxes[rid]) for rid in CHAIN + LINKS + ["R1839"]}}
    result["new"]["forced_exit_without"] = {}
    for rid in LINKS[:-1]:
        with model:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            result["new"]["forced_exit_without"][rid] = solve(model)
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: result[k] for k in ("baseline", "new")}, indent=1))
print("LPs:", lp_count)
