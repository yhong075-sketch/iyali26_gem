"""Step 1: capacity and need of the eight reactions on 7ccbb568 as is (in memory; nothing is written to SBML).

Usage: python step1_capacity.py   ->  out/step1_capacity.json

A. FVA without a growth floor of EIGHT + NEIGH, conditions open and taskA, Gurobi; GLPK repeats EIGHT.
B. taskA, Gurobi: WT growth; FVA of EIGHT with growth >= 0.999999 * WT and >= 0.999 * WT; one pFBA.
C. taskA, Gurobi: WT growth with each of EIGHT closed alone, and with R2025 + R2113 both closed.
D. taskA, Gurobi: targeted knockouts of KO_GENES as is and with R2025, R2113, or both closed.
   A knockout LP that is not optimal is re-solved once in a fresh context; if it is infeasible again, the maximum of
   xMAINTENANCE with its lower bound at 0 is recorded (benchmark.md certification items 1-3), never read as 0.
Closing a reaction that is blocked (min = max = 0 without a growth floor) cannot change any knockout, because a
knockout only removes feasible points; so D only closes the reactions that part A finds unblocked on taskA.
"""
from common import (CUTOFFS, EIGHT, KO_GENES, LP_COUNT, NEIGH, TOL, classify, dump, fva, load, lp, pfba)

result = {"A": {}, "B": {}, "C": {}, "D": {}}

# A. capacity, no growth floor
for solver in ["gurobi", "glpk"]:
    for condition in ["open", "taskA"]:
        model = load(condition, solver)
        rx = EIGHT + NEIGH if solver == "gurobi" else EIGHT
        result["A"][f"{condition}/{solver}"] = fva(model, rx)
        print("A", condition, solver, {k: (v.get("min"), v.get("max")) for k, v in
                                       result["A"][f"{condition}/{solver}"].items() if k in EIGHT}, flush=True)

model = load("taskA", "gurobi")

# B. at the WT optimum
wt, s = lp(model)
result["B"]["wt"] = {"value": wt, "status": s}
for name, frac in [("fva_999999", 0.999999), ("fva_999", 0.999)]:
    with model:
        model.reactions.biomass_C.lower_bound = frac * wt
        result["B"][name] = fva(model, EIGHT)
growth, s2, flux = pfba(model)
result["B"]["pfba"] = {"growth": growth, "status": s2,
                       "flux": {k: v for k, v in (flux or {}).items() if abs(v) > TOL}}
print("B wt", wt, {k: (v["min"], v["max"]) for k, v in result["B"]["fva_999999"].items()}, flush=True)

# C. closures
closures = {rid: [rid] for rid in EIGHT}
closures["R2025+R2113"] = ["R2025", "R2113"]
for name, rids in closures.items():
    with model:
        for rid in rids:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        v, st = lp(model)
    result["C"][name] = {"wt": v, "status": st, "change": None if v is None else v - wt}
print("C", {k: v["change"] for k, v in result["C"].items()}, flush=True)

# D. targeted knockouts
unblocked = [rid for rid in EIGHT if not result["A"]["taskA/gurobi"][rid]["blocked"]]
variants = {"as_is": []}
for rid in unblocked:
    variants[f"{rid}_closed"] = [rid]
if {"R2025", "R2113"} <= set(unblocked):
    variants["R2025+R2113_closed"] = ["R2025", "R2113"]
result["D"]["unblocked_on_taskA"] = unblocked
result["D"]["variants"] = variants


def knockout(model, gid):
    with model:
        model.genes.get_by_id(gid).knock_out()
        g, st = lp(model)
        rec = {"growth": g, "status": st}
        if g is None:
            g2, st2 = lp(model)  # fresh re-solve (basis discarded inside lp)
            rec.update({"resolve_growth": g2, "resolve_status": st2})
            if g2 is None and st2 == "infeasible":
                forced = [r.id for r in model.reactions if r.lower_bound > 0 or r.upper_bound < 0]
                rec["forced_bounds"] = forced
                m = model.reactions.xMAINTENANCE
                lb = m.lower_bound
                m.lower_bound = 0.0
                model.objective = "xMAINTENANCE"
                mx, stm = lp(model)
                rec.update({"maint_lb": lb, "maint_max": mx, "maint_status": stm})
            else:
                g = g2
                rec["growth"] = g2
    return rec


for vname, rids in variants.items():
    with model:
        for rid in rids:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        vwt, vst = lp(model)
        kos = {}
        for gid in KO_GENES:
            rec = knockout(model, gid)
            rec["ratio"] = None if (rec["growth"] is None or not vwt) else rec["growth"] / vwt
            rec["class"] = {str(c): x for c, x in classify(rec["growth"], vwt).items()}
            kos[gid] = rec
    result["D"][vname] = {"wt": vwt, "wt_status": vst, "ko": kos}
    print("D", vname, vwt, {g: (round(r["ratio"], 6) if r["ratio"] is not None else r["status"])
                            for g, r in kos.items()}, flush=True)

result["cutoffs"] = list(CUTOFFS)
result["lp_count"] = LP_COUNT["n"]
print("LPs", LP_COUNT["n"], "->", dump(result, "step1_capacity.json"))
