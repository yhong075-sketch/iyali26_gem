"""Phase 3.2 of TASK.md: checked solves on the baseline a1237821 and the built candidate (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/cofactor_biomass_terms_20261027/scripts/solves.py MEDIUM.csv PROFILE.json OUT.json MODEL.xml...
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py). Per model: checked WT,
pFBA with the route reactions and the pFBA net producers of the six cofactors, the six biomass coefficients,
WT with each vitamin uptake (R946 pantothenate, R1268 nicotinate, R1314 riboflavin, R1127 folate) closed and
with the three closed together, and the no-uptake energy check (max R_NTP1 with every exchange at (0, 1000)
and the NGAM lower bound of xMAINTENANCE relaxed to 0; see phase1_egc.py for why the census method is
infeasible here).
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

TERMS = {"m479[C_cy]": "THF", "m528[C_cy]": "FAD", "m926[C_cy]": "riboflavin", "m122[C_cy]": "NAD", "m41[C_cy]": "NADP", "m69[C_cy]": "CoA"}
VITAMINS = {"R946": "pantothenate", "R1268": "nicotinate", "R1314": "riboflavin", "R1127": "folate"}
WATCH = ["R377", "R274", "R14", "R277", "R278", "R273", "R749", "R271", "R272", "R849", "R811",
         "R378", "R34", "R12", "R11", "R699", "R700", "R1316", "R1315", "R1314", "R697", "R698", "R306", "R847",
         "R946", "R616", "R1908", "R1909", "R647", "R648", "R615", "R1092", "R1291", "R265",
         "R1268", "R907", "R576", "R578", "R579", "R580", "R568", "R569", "R565", "R566", "R571", "R572"]
medium, profile, out, *paths = sys.argv[1:]
result, lp_count = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve)", "models": {}}, 0


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


def producers(model, solution, mid, limit=10):
    met = model.metabolites.get_by_id(mid)
    rows = []
    for r in met.reactions:
        v = float(solution.fluxes[r.id])
        net = r.metabolites[met] * v
        if abs(net) > 1e-12:
            rows.append({"reaction": r.id, "name": r.name[:50], "flux": v, "net_production": net, "gpr": r.gene_reaction_rule[:60]})
    rows.sort(key=lambda row: -row["net_production"])
    return rows[:limit]


def egc_relaxed(model):
    global lp_count
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = []
        for r in model.reactions:
            if r.lower_bound > 0:
                relaxed.append([r.id, r.lower_bound])
                r.lower_bound = 0.0
        model.objective = model.reactions.get_by_id("R_NTP1")
        lp_count += 1
        value = model.slim_optimize(error_value=float("nan"))
        return {"status": model.solver.status, "max_R_NTP1_no_uptake": value, "relaxed_positive_lower_bounds": relaxed}


for path in paths:
    ctx = context(path)
    model = ctx.model
    wt = solve(model, matrix(model), ())
    solution = pfba(model)
    lp_count += 3
    biomass = {m.id: c for m, c in model.reactions.biomass_C.metabolites.items()}
    row = {"file": path.split("/")[-1], "model_sha256": ctx.canonical_model_sha256,
           "medium_sha256": ctx.medium_sha256, "strain_profile_sha256": ctx.strain_profile_sha256,
           "biomass_terms": {mid: biomass.get(mid, 0.0) for mid in TERMS}, "biomass_metabolite_count": len(biomass),
           "wt": wt,
           "pfba": {"status": solution.status, "growth": float(solution.fluxes["biomass_C"]),
                    "watched_fluxes": {rid: float(solution.fluxes[rid]) for rid in WATCH if rid in model.reactions},
                    "producers": {label: producers(model, solution, mid) for mid, label in TERMS.items()}},
           "vitamins_closed": {}, "egc_relaxed": egc_relaxed(model)}
    for rid, label in VITAMINS.items():
        m = context(path).model
        m.reactions.get_by_id(rid).lower_bound = 0.0
        row["vitamins_closed"][f"{rid} ({label})"] = solve(m, matrix(m), ())
        lp_count += 1
    m = context(path).model
    for rid in ("R946", "R1268", "R1314"):
        m.reactions.get_by_id(rid).lower_bound = 0.0
    row["vitamins_closed"]["R946+R1268+R1314"] = solve(m, matrix(m), ())
    lp_count += 1
    result["models"][ctx.canonical_model_sha256[:8]] = row
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for sha, row in result["models"].items():
    print(sha, row["file"][-70:], "WT", row["wt"]["status"], repr(row["wt"]["growth"]), "passed", row["wt"]["passed"],
          "| biomass n", row["biomass_metabolite_count"], "| EGC", row["egc_relaxed"]["status"], row["egc_relaxed"]["max_R_NTP1_no_uptake"])
    for label, r in row["vitamins_closed"].items():
        print(f"   closed {label:26s} {r['status']} {r['growth']!r}")
    nonzero = {k: f"{v:.3g}" for k, v in row["pfba"]["watched_fluxes"].items() if abs(v) > 1e-12}
    print("   pFBA watched nonzero:", nonzero)
print("LPs", lp_count)
