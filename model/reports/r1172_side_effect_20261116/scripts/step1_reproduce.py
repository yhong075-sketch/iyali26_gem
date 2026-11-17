"""Step 1: are the nine blocked on 7ccbb568, and does reopening R1172 in memory free them?

Usage: python step1_reproduce.py SOLVER   (gurobi or glpk)
For each condition (open, taskA) and variant (as_is, r1172_open): FVA without growth floor of the nine and the
CoA-moiety boundary reactions, plus the growth optimum (taskA only). Writes out/step1_<solver>.json.
"""
import sys

from common import BOUNDARY, LP_COUNT, NINE, R1172_OLD_BOUNDS, dump, fva, load, lp

solver = sys.argv[1]
targets = list(dict.fromkeys(NINE + BOUNDARY))
result = {"solver": solver, "targets": targets, "runs": {}}
for condition in ["open", "taskA"]:
    model = load(condition, solver)
    for variant in ["as_is", "r1172_open"]:
        with model:
            if variant == "r1172_open":
                model.reactions.R1172.bounds = R1172_OLD_BOUNDS
            run = {"R1172_bounds": list(model.reactions.R1172.bounds)}
            if condition == "taskA":
                growth, status = lp(model)
                run["growth"], run["growth_status"] = growth, status
            run["fva"] = fva(model, targets)
            run["nine_blocked"] = [r for r in NINE if run["fva"][r]["blocked"]]
            run["nine_unresolved"] = [r for r in NINE if run["fva"][r]["blocked"] is None]
            result["runs"][f"{condition}/{variant}"] = run
            print(condition, variant, "blocked of nine:", len(run["nine_blocked"]),
                  "unresolved:", run["nine_unresolved"], "growth:", run.get("growth"), flush=True)
result["lp_count"] = LP_COUNT["n"]
print("LPs", LP_COUNT["n"], "->", dump(result, f"step1_{solver}.json"))
