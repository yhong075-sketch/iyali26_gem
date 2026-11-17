"""Step 2, LP part (Gurobi 1e-9; GLPK cross-check where noted). Writes out/step2_<solver>.json.

Usage: python step2_lp.py SOLVER [--full]
  without --full (GLPK cross-check): part A on the as-is model and part B, both conditions.
  with --full (Gurobi): parts A-D.

A. Ranges (no growth floor) of the reactions around the matrix CoA pool that the static scan names, as is and
   with R1172 reopened.
B. Net propionyl-CoA moiety import into the matrix, P = v(R2025) - v(R2113) (R2025 forward = cytosol to matrix,
   R2113 forward = matrix to cytosol). Because no matrix reaction makes or destroys the moiety (static scan),
   P equals the net moiety export through the other channels at every steady state. Max/min of P, max of R834,
   and max HMG-CoA export -v(R1172).
C. Supporting route of each of the nine: xMAINTENANCE lower bound 0, target fixed at min(1, half its maximum),
   minimum total absolute flux. Variants: R1172 reopened; CoA drain and FAD drain as in the gap census
   (7ccbb568_open_draintests.json, irreversible demands).
D. Before the closure, at the WT optimum (Task A, R1172 reopened): ranges of the nine at 99.9 % of the optimum,
   and the pFBA flux of R1172, R834 and the nine.
"""
import sys

from cobra.flux_analysis import pfba
from optlang.symbolics import Zero

from common import (LP_COUNT, NINE, R1172_OLD_BOUNDS, TOL, dump, extreme, flux_vector, fva, load, lp)

solver = sys.argv[1]
full = "--full" in sys.argv
AROUND = ["R1291", "R1292", "R265", "R1969", "R2118", "R2117", "R2116", "R2115", "R2114", "R847", "R2142",
          "R2082", "R2083", "R2084", "R2024", "R2022", "R2025", "R2113", "R2285"]
P = {"R2025": 1.0, "R2113": -1.0}
COA_DRAIN = ["m69[C_cy]", "m73[C_mi]", "m182[C_pe]", "m207[C_nu]", "m233[C_em]", "m252[C_er]", "m573[C_lp]",
             "m792[C_mm]", "m1494[C_en]"]
FAD_DRAIN = ["m528[C_cy]", "m563[C_mi]"]
# delta of matrix moiety per unit forward flux (static scan, tables/matrix_moiety_reactions.csv)
DELTA = {"R834": 1, "R1172": 1, "R1291": 1, "R2025": 1, "R2077": -1, "R2083": -1, "R2084": -1, "R2113": -1,
         "R2118": -1}


def add_drains(model, mets):
    import cobra
    for i, mid in enumerate(mets):
        d = cobra.Reaction(f"DRAIN_{i}", lower_bound=0, upper_bound=1000)
        model.add_reactions([d])
        d.add_metabolites({model.metabolites.get_by_id(mid): -1})


def min_total_flux(model, target, value):
    """Fix target at value and minimise the sum of absolute fluxes. One LP."""
    with model:
        rxn = model.reactions.get_by_id(target)
        rxn.bounds = (value, value)
        coefficients = {}
        for r in model.reactions:
            coefficients[r.forward_variable] = 1.0
            coefficients[r.reverse_variable] = 1.0
        model.objective = model.problem.Objective(Zero, direction="min", sloppy=True)  # as cobra's add_pfba
        model.objective.set_linear_coefficients(coefficients)
        total, status = lp(model)
        flux = flux_vector(model) if total is not None else None
    return total, status, flux


result = {"solver": solver, "A": {}, "B": {}, "C": {}, "D": {}}
for condition in ["open", "taskA"]:
    model = load(condition, solver)
    for variant in ["as_is", "r1172_open"]:
        with model:
            if variant == "r1172_open":
                model.reactions.R1172.bounds = R1172_OLD_BOUNDS
            key = f"{condition}/{variant}"
            if full or variant == "as_is":  # GLPK cross-checks part A on the as-is model only (LP budget)
                result["A"][key] = fva(model, AROUND)
            b = {}
            for name, obj, sense in [("P_max", P, "max"), ("P_min", P, "min"), ("R834_max", "R834", "max"),
                                     ("HMGCoA_export_max", "R1172", "min")]:
                v, s = extreme(model, obj, sense)
                b[name] = {"value": v, "status": s}
            result["B"][key] = b
            print(key, {k: (round(v["value"], 6) if v["value"] is not None else None, v["status"])
                        for k, v in b.items()}, flush=True)
    if not full:
        continue
    for variant in ["r1172_open", "coa_drain", "fad_drain"]:
        with model:
            model.reactions.xMAINTENANCE.lower_bound = 0.0
            if variant == "r1172_open":
                model.reactions.R1172.bounds = R1172_OLD_BOUNDS
            elif variant == "coa_drain":
                add_drains(model, COA_DRAIN)
            else:
                add_drains(model, FAD_DRAIN)
            routes = {}
            for target in NINE:
                hi, shi = extreme(model, target, "max")
                lo, slo = extreme(model, target, "min")
                if hi is None or lo is None:
                    routes[target] = {"max": hi, "max_status": shi, "min": lo, "min_status": slo}
                    continue
                # carry the target in its usable direction
                value = min(1.0, 0.5 * hi) if hi > TOL else max(-1.0, 0.5 * lo)
                if abs(value) <= TOL:
                    routes[target] = {"max": hi, "min": lo, "note": "blocked in this variant"}
                    continue
                total, status, flux = min_total_flux(model, target, value)
                entry = {"max": hi, "min": lo, "fixed_at": value, "total_abs_flux": total, "status": status}
                if flux is not None:
                    nz = {k: v for k, v in flux.items() if abs(v) > 1e-9}
                    entry["nonzero"] = dict(sorted(nz.items()))
                    entry["moiety_into_matrix_by_channel"] = {k: DELTA[k] * flux[k] for k in DELTA
                                                              if abs(flux[k]) > 1e-9}
                    entry["moiety_net_into_matrix"] = sum(DELTA[k] * flux[k] for k in DELTA)
                routes[target] = entry
            result["C"][f"{condition}/{variant}"] = routes
            print(condition, variant, {t: (r.get("fixed_at"), r.get("status"), r.get("moiety_into_matrix_by_channel"))
                                       for t, r in routes.items()}, flush=True)
    if condition == "taskA":
        with model:
            model.reactions.R1172.bounds = R1172_OLD_BOUNDS
            growth, status = lp(model)
            d = {"growth": growth, "growth_status": status}
            with model:
                model.reactions.biomass_C.lower_bound = 0.999 * growth
                d["fva_999"] = fva(model, NINE)
            LP_COUNT["n"] += 2  # pfba = one optimisation plus one minimisation
            sol = pfba(model)
            d["pfba_status"] = sol.status
            d["pfba"] = {k: float(sol.fluxes[k]) for k in NINE + ["R1172", "R834", "R83", "R412", "R410", "R411",
                                                                    "R82", "R2025", "R2113", "biomass_C"]}
            result["D"] = d
            print("D", d["growth"], {k: (v["min"], v["max"]) for k, v in d["fva_999"].items()}, d["pfba"])
result["lp_count"] = LP_COUNT["n"]
print("LPs", LP_COUNT["n"], "->", dump(result, f"step2_{solver}.json"))
