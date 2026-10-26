"""Diagnostic LPs for the GPI co-product change (read-only; models are loaded in memory).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/gpi_coproducts_20261021/scripts/solves.py \
    BASELINE.xml NEW.xml RESEARCH_ROOT OUT.json
BASELINE is the links candidate (89559a0f), NEW adds the co-products. Condition: SD-Leu + PO1f, GLPK.
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


def forced(model):
    global lp_count
    with model:
        model.reactions.R_DM_GPI.lower_bound = FORCED
        lp_count += 2
        solution = pfba(model)
        return {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
                "fluxes": {rid: float(solution.fluxes[rid]) for rid in CHAIN + LINKS + ["R1839"]}}


result = {"condition": "SD-Leu + PO1f, GLPK", "forced_exit": FORCED}
for name, path in (("baseline", baseline_path), ("new", new_path)):
    ctx = context(path)
    model = ctx.model
    result[name] = {"model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
                    "strain_profile_sha256": ctx.strain_profile_sha256, "wt": solve(model),
                    "fva_fraction_0": fva(model, CHAIN + LINKS), "forced_exit_pfba": forced(model)}
    if name == "new":
        result[name]["forced_exit_without"] = {}
        for rid in LINKS[:-1]:
            with model:
                model.reactions.R_DM_GPI.lower_bound = FORCED
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                result[name]["forced_exit_without"][rid] = solve(model)
for name in ("baseline", "new"):
    result[name]["forced_growth_cost"] = result[name]["wt"]["growth"] - result[name]["forced_exit_pfba"]["growth"]
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for name in ("baseline", "new"):
    r = result[name]
    print(name, r["model_sha256"][:8], "WT", r["wt"], "forced cost", r["forced_growth_cost"])
    print("  FVA max", {k: round(v[1], 6) for k, v in r["fva_fraction_0"].items()})
    print("  FVA min", {k: v[0] for k, v in r["fva_fraction_0"].items() if abs(v[0]) > 1e-9})
    print("  forced", r["forced_exit_pfba"]["status"], {k: float(f"{v:.4g}") for k, v in r["forced_exit_pfba"]["fluxes"].items()})
print("new, forced exit without each link:", result["new"]["forced_exit_without"])
print("LPs:", lp_count)
