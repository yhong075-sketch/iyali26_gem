"""Checks on the built files: 9d655c2b (control copy) vs the R2171 1:1 candidate (read-only, about 18 LPs).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/r2171_pi_plc_1to1_20261028/scripts/built_checks.py BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Pass short-path byte-identical copies: cobra on macOS rejects absolute paths of 260+ characters.
Task A medium + PO1f, Gurobi 1e-9 (checked_solve.py). Per model: R2171 equation and charge balance; WT (checked);
R2171 FVA at 99.9 % and at growth >= 0 (processes=1); no-uptake LPs with every exchange (0, 1000) and every positive
lower bound relaxed to 0 (max R_NTP1 and max biomass_C, status recorded).
"""
import hashlib
import json
import math
import sys
import warnings
from pathlib import Path

from cobra.flux_analysis import flux_variability_analysis

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

baseline_path, candidate_path, medium, profile, out = sys.argv[1:6]


def no_uptake(model, objective):
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = [[r.id, r.lower_bound] for r in model.reactions if r.lower_bound > 0]
        for rid, _ in relaxed:
            model.reactions.get_by_id(rid).lower_bound = 0.0
        model.objective = objective
        model.solver.problem.reset()
        value = model.slim_optimize(error_value=float("nan"))
        return {"objective": objective, "status": model.solver.status,
                "value": None if math.isnan(value) else value, "relaxed": relaxed}


result = {}
for label, path in (("baseline", baseline_path), ("candidate", candidate_path)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile).model
    configure(model)
    r = model.reactions.R2171
    row = {"path": path, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
           "r2171": r.build_reaction_string(), "r2171_charge_balance": sum(m.charge * c for m, c in r.metabolites.items()),
           "wt": solve(model, matrix(model), set()), "fva": {}}
    for fraction in (0.999, 0.0):
        model.solver.problem.reset()
        frame = flux_variability_analysis(model, reaction_list=["R2171"], fraction_of_optimum=fraction, processes=1)
        row["fva"][str(fraction)] = [float(frame.loc["R2171", "minimum"]), float(frame.loc["R2171", "maximum"])]
    row["no_uptake"] = [no_uptake(model, "R_NTP1"), no_uptake(model, "biomass_C")]
    result[label] = row
    print(label, row["sha256"][:8], row["r2171"], "charge", row["r2171_charge_balance"], row["wt"]["status"],
          row["wt"]["growth"], row["wt"]["passed"], row["fva"], [(x["objective"], x["status"], x["value"]) for x in row["no_uptake"]],
          file=sys.stderr, flush=True)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
