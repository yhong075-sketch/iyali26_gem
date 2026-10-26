"""Reproduce the Claude Science TPP figure in memory on one model (read-only).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/tpp_requirement_20261022/scripts/reproduce_figure.py \
    MODEL.xml MEDIUM.csv PROFILE.json OUT.json
Condition: Task A medium + PO1f, GLPK. The TPP term and the closures are applied in memory only.
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

model_path, media_path, profile_path, out = sys.argv[1:5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=media_path,
                                            strain_profile_path=profile_path)
model = ctx.model
model.solver = "glpk"
tpp = model.metabolites.get_by_id("m975[C_cy]")
result = {"condition": "Task A medium + PO1f, GLPK", "model_sha256": ctx.canonical_model_sha256,
          "medium_sha256": ctx.medium_sha256, "strain_profile_sha256": ctx.strain_profile_sha256,
          "R1340_bounds": list(model.reactions.R1340.bounds)}
lp_count = 0


def run(label, coefficient=0.0, knockout=False, close=()):
    global lp_count
    with model:
        if coefficient:
            model.reactions.biomass_C.add_metabolites({tpp: -coefficient})
        if knockout:
            model.genes.get_by_id(THI80).knock_out()
        for rid in close:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        growth = model.slim_optimize(error_value=float("nan"))
        lp_count += 1
        row = {"tpp_coefficient": coefficient, "thi80_ko": knockout, "closed": list(close),
               "status": model.solver.status, "growth": None if math.isnan(growth) else growth}
        if not math.isnan(growth) and growth > 1e-9:
            solution = pfba(model)
            lp_count += 2
            row["pfba_nonzero"] = {r: float(solution.fluxes[r]) for r in MAP if abs(solution.fluxes[r]) > 1e-12}
        result[label] = row


run("as_supplied")
table = flux_variability_analysis(model, MAP, fraction_of_optimum=0.0, processes=1)
lp_count += 2 * len(MAP)
result["as_supplied_fva_fraction_0"] = {r: [float(table.loc[r, "minimum"]), float(table.loc[r, "maximum"])]
                                        for r in MAP}
run("thi80_ko", knockout=True)
run("tpp_1e-5", 1e-5)
run("tpp_1e-5_thi80_ko", 1e-5, knockout=True)
run("tpp_1e-5_thi80_ko_R1342_closed", 1e-5, knockout=True, close=("R1342",))
run("tpp_1e-5_R1342_closed", 1e-5, close=("R1342",))
run("tpp_1e-5_no_thiamine", 1e-5, close=("R1340",))
run("tpp_1.2e-6_R1342_closed", 1.2e-6, close=("R1342",))
run("tpp_1.2e-6_thi80_ko_R1342_closed", 1.2e-6, knockout=True, close=("R1342",))
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for key, row in result.items():
    if isinstance(row, dict) and "status" in row:
        print(f"{key:36s} {row['status']:8s} {row['growth']!r} {row.get('pfba_nonzero', '')}")
