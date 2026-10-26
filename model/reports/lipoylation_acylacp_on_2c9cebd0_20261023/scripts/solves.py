"""WT and key fluxes on 2c9cebd0 and the new candidate (read-only; 3 LPs per model).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/lipoylation_acylacp_on_2c9cebd0_20261023/scripts/solves.py MEDIUM.csv PROFILE.json OUT.json MODEL.xml...
Copy of ../../combined_tpp_heme_r2193_20261022/scripts/solves.py with the lipoylation and acyl-ACP reactions
watched: Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9; per model a checked WT solve, then pFBA.
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

WATCH = ["R78", "R539", "R362", "R364", "R2286", "R2287", "R2288", "R2289", "R2290", "R1394", "R1395", "R1396",
         *[f"R{i}" for i in range(2216, 2242)], "R2193", "R1189"]
BIOMASS_TERMS = ("m586[C_mi]", "m2062[C_mi]", "m2065[C_mi]", "m2064[C_mi]", "m975[C_cy]", "m1635[C_cy]", "m948[C_cy]")
medium, profile, out, *paths = sys.argv[1:]
result, lp_count = {"condition": "Task A medium + PO1f, Gurobi 1e-9", "models": {}}, 0
for path in paths:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    configure(model)
    wt = solve(model, matrix(model), ())
    solution = pfba(model)
    lp_count += 3
    biomass = {m.id: c for m, c in model.reactions.biomass_C.metabolites.items()}
    result["models"][ctx.canonical_model_sha256[:8]] = {
        "file": path.split("/")[-1], "model_sha256": ctx.canonical_model_sha256,
        "medium_sha256": ctx.medium_sha256, "strain_profile_sha256": ctx.strain_profile_sha256,
        "biomass_terms": {mid: biomass.get(mid, 0.0) for mid in BIOMASS_TERMS},
        "wt": wt, "pfba": {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
                           "fluxes": {rid: float(solution.fluxes[rid]) for rid in WATCH if rid in model.reactions}}}
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for sha, row in result["models"].items():
    nonzero = {k: f"{v:.3g}" for k, v in row["pfba"]["fluxes"].items() if abs(v) > 1e-12}
    print(sha, row["file"][-60:], "WT", row["wt"]["status"], repr(row["wt"]["growth"]), "passed", row["wt"]["passed"],
          "| pFBA nonzero", nonzero)
print("LPs", lp_count)
