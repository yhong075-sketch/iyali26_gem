"""Diagnostic LPs for the TPP requirement change (read-only; models are loaded in memory).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/tpp_requirement_20261022/scripts/solves.py \
    BASELINE.xml NEW.xml MEDIUM.csv PROFILE.json OUT.json
BASELINE is cf0825b8, NEW adds --tpp-requirement. Condition: Task A medium + PO1f, GLPK.
"""
import json
import math
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis, pfba

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

THI80 = "YALI1E25443g"
MAP = ["R934", "R1338", "R752", "R753", "R933", "R1337", "R1342", "R754", "R408", "R413", "R645",
       "R970", "R1343", "R1344", "R1340"]

baseline_path, new_path, media_path, profile_path, out = sys.argv[1:6]
lp_count = 0


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=media_path,
                                                strain_profile_path=profile_path)
    ctx.model.solver = "glpk"
    return ctx


def solve(model, label, knockout=False, close=(), reopen=()):
    """Growth under the given in-memory edits, with the pFBA fluxes of the thiamine map when it grows."""
    global lp_count
    with model:
        if knockout:
            model.genes.get_by_id(THI80).knock_out()
        for rid in close:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        for rid in reopen:
            model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)
        growth = model.slim_optimize(error_value=float("nan"))
        lp_count += 1
        row = {"label": label, "thi80_ko": knockout, "closed": list(close), "reopened": list(reopen),
               "status": model.solver.status, "growth": None if math.isnan(growth) else growth}
        if not math.isnan(growth) and growth > 1e-9:
            solution = pfba(model)
            lp_count += 2
            row["pfba_nonzero"] = {r: float(solution.fluxes[r]) for r in MAP if abs(solution.fluxes[r]) > 1e-12}
        return row


result = {"condition": "Task A medium + PO1f, GLPK"}
for name, path in (("baseline", baseline_path), ("new", new_path)):
    ctx = context(path)
    model = ctx.model
    table = flux_variability_analysis(model, MAP, fraction_of_optimum=0.0, processes=1)
    lp_count += 2 * len(MAP)
    rows = [solve(model, "wt"), solve(model, "thi80_ko", knockout=True),
            solve(model, "no_thiamine", close=("R1340",)),
            solve(model, "no_thiamine_thi80_ko", knockout=True, close=("R1340",))]
    if name == "new":
        # The closure is the cause: reopening R1342 in memory restores the bypass.
        rows.append(solve(model, "thi80_ko_R1342_reopened", knockout=True, reopen=("R1342",)))
    result[name] = {"model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
                    "strain_profile_sha256": ctx.strain_profile_sha256,
                    "biomass_C_tpp": model.reactions.biomass_C.get_coefficient("m975[C_cy]")
                    if model.metabolites.get_by_id("m975[C_cy]") in model.reactions.biomass_C.metabolites else 0.0,
                    "R1342_bounds": list(model.reactions.R1342.bounds),
                    "fva_fraction_0": {r: [float(table.loc[r, "minimum"]), float(table.loc[r, "maximum"])] for r in MAP},
                    "solves": {row.pop("label"): row for row in rows}}
for name in ("baseline", "new"):
    solves = result[name]["solves"]
    wt = solves["wt"]["growth"]
    for row in solves.values():
        row["ratio_to_wt"] = None if row["growth"] is None else row["growth"] / wt
result["wt_growth_change"] = result["new"]["solves"]["wt"]["growth"] - result["baseline"]["solves"]["wt"]["growth"]
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for name in ("baseline", "new"):
    print(name, result[name]["model_sha256"][:8], "TPP", result[name]["biomass_C_tpp"], "R1342", result[name]["R1342_bounds"])
    for label, row in result[name]["solves"].items():
        print(f"  {label:26s} {row['status']:8s} {row['growth']!r} ratio {row['ratio_to_wt']!r} {row.get('pfba_nonzero', '')}")
    print("  FVA max nonzero:", {r: v[1] for r, v in result[name]["fva_fraction_0"].items() if abs(v[1]) > 1e-9})
print("WT change", result["wt_growth_change"], "LPs", lp_count)
