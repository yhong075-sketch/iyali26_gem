"""Phase 1 of TASK.md: in-memory checks of the six cofactor biomass terms before any build (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/cofactor_biomass_terms_20261027/scripts/phase1_checks.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
Task A medium + PO1f, Gurobi at FeasibilityTol/OptimalityTol 1e-9 (checked_solve.py). Solves:
  baseline WT; all six terms WT + pFBA; each term alone (6); all six with the vitamin uptakes R946 / R1268 /
  R1314 closed one at a time (3) and together (1); folate R1127 closed (1); energy-generating-cycle LP on
  the baseline and the six-term model (max ATP hydrolysis with every exchange closed, 2).
For each cofactor the pFBA net producers of the cytosolic species are listed (reaction, flux, net production).
"""
import json
import sys
import warnings

from cobra import Reaction
from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

TERMS = {"m479[C_cy]": ("THF", 6.34e-5), "m528[C_cy]": ("FAD", 1e-5), "m926[C_cy]": ("riboflavin", 9.9e-4),
         "m122[C_cy]": ("NAD", 2.65e-3), "m41[C_cy]": ("NADP", 5.7e-4), "m69[C_cy]": ("CoA", 1.9e-4)}
VITAMINS = {"R946": "pantothenate", "R1268": "nicotinate", "R1314": "riboflavin", "R1127": "folate"}
WATCH = ["R377", "R274", "R14", "R277", "R278", "R273", "R749", "R271", "R272", "R849", "R811",
         "R378", "R34", "R12", "R11", "R699", "R700", "R1316", "R1315", "R1314", "R697", "R698", "R306", "R847",
         "R946", "R616", "R1908", "R1909", "R647", "R648", "R615", "R1092", "R1291", "R265",
         "R1268", "R907", "R576", "R578", "R579", "R580", "R568", "R569", "R565", "R566", "R571", "R572"]
model_path, medium, profile, out = sys.argv[1:5]
lp_count = 0


def context():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


def add_terms(model, mids):
    model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id(mid): -TERMS[mid][1] for mid in mids})


def checked(model):
    global lp_count
    lp_count += 1
    return solve(model, matrix(model), ())


def producers(model, solution, mid, limit=12):
    met = model.metabolites.get_by_id(mid)
    rows = []
    for r in met.reactions:
        v = float(solution.fluxes[r.id])
        net = r.metabolites[met] * v
        if abs(net) > 1e-12:
            rows.append({"reaction": r.id, "name": r.name[:50], "flux": v, "net_production": net,
                         "gpr": r.gene_reaction_rule[:60]})
    rows.sort(key=lambda row: -row["net_production"])
    return rows[:limit] + ([{"note": f"{len(rows) - limit} more consumers/producers omitted"}] if len(rows) > limit else [])


def egc(model):
    """Census closure_check method: max R_NTP1 (ATP maintenance, lower bound 0) with every exchange at
    (0, 1000), i.e. no uptake; a positive value means an energy-generating cycle."""
    global lp_count
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        atp = model.reactions.get_by_id("R_NTP1")
        atp.lower_bound = 0.0
        model.objective = atp
        lp_count += 1
        value = model.slim_optimize(error_value=float("nan"))
        return {"status": model.solver.status, "objective": "R_NTP1", "max_atp_no_uptake": value}


result = {"condition": "Task A medium + PO1f, Gurobi 1e-9 (checked_solve)", "model": model_path, "terms": {k: v[1] for k, v in TERMS.items()}}
ctx = context()
result["model_sha256"], result["medium_sha256"], result["strain_profile_sha256"] = ctx.canonical_model_sha256, ctx.medium_sha256, ctx.strain_profile_sha256
base = ctx.model
result["baseline_wt"] = checked(base)
result["baseline_egc"] = egc(base)
# 1 + 3: all six terms, WT and pFBA
six = context().model
add_terms(six, TERMS)
result["six_terms_wt"] = checked(six)
sol = pfba(six)
lp_count += 2
result["six_terms_pfba"] = {"status": sol.status, "growth": float(sol.fluxes["biomass_C"]),
                            "watched_fluxes": {rid: float(sol.fluxes[rid]) for rid in WATCH if rid in six.reactions},
                            "exchanges_nonzero": {r.id: float(sol.fluxes[r.id]) for r in six.exchanges if abs(sol.fluxes[r.id]) > 1e-9},
                            "producers": {TERMS[mid][0]: producers(six, sol, mid) for mid in TERMS}}
result["six_terms_egc"] = egc(six)
# 2: each term alone
result["single_terms"] = {}
for mid, (label, coefficient) in TERMS.items():
    one = context().model
    add_terms(one, [mid])
    row = checked(one)
    row["growth_cost_vs_baseline"] = (result["baseline_wt"]["growth"] - row["growth"]) if row["growth"] is not None else None
    result["single_terms"][label] = row
# 4: vitamin uptakes closed, all six terms
result["vitamins_closed"] = {}
for rid, label in VITAMINS.items():
    m = context().model
    add_terms(m, TERMS)
    m.reactions.get_by_id(rid).lower_bound = 0.0
    result["vitamins_closed"][f"{rid} ({label})"] = checked(m)
m = context().model
add_terms(m, TERMS)
for rid in ("R946", "R1268", "R1314"):
    m.reactions.get_by_id(rid).lower_bound = 0.0
result["vitamins_closed"]["R946+R1268+R1314"] = checked(m)
if result["vitamins_closed"]["R946+R1268+R1314"]["growth"]:
    sol3 = pfba(m)
    lp_count += 2
    result["vitamins_closed_pfba_producers"] = {TERMS[mid][0]: producers(m, sol3, mid) for mid in TERMS}
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print("baseline WT", result["baseline_wt"]["status"], repr(result["baseline_wt"]["growth"]), "EGC", result["baseline_egc"])
print("six terms WT", result["six_terms_wt"]["status"], repr(result["six_terms_wt"]["growth"]), "passed", result["six_terms_wt"]["passed"], "EGC", result["six_terms_egc"])
for label, row in result["single_terms"].items():
    print(f"  {label:11s} alone: {row['status']} growth {row['growth']!r} cost {row['growth_cost_vs_baseline']!r}")
for label, row in result["vitamins_closed"].items():
    print(f"  closed {label:28s}: {row['status']} growth {row['growth']!r}")
print("pFBA producers (six terms):")
for label, rows in result["six_terms_pfba"]["producers"].items():
    print(" ", label, [(r.get("reaction"), f"{r.get('net_production', 0):.3g}") for r in rows if "reaction" in r][:8])
print("LPs", lp_count)
