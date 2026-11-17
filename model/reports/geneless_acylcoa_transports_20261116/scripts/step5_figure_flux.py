"""Step 5: flux solutions for the schematic (Task A + PO1f, Gurobi 1e-9; in memory only).

Usage: python step5_figure_flux.py   ->  out/step5_figure_flux.json

- asis_route : 7ccbb568 as is, R2113 fixed at -1 (one propionyl-CoA into the matrix), xMAINTENANCE lower bound 0,
               minimum total flux. Shows the pass-through R2113 -> R2025 that the transports carry today.
- opt_route  : B0/ALL (step 4), R25 fixed at 1 and R2024 fixed at -1, xMAINTENANCE lower bound 0, minimum total
               flux. Shows the valine route inside the matrix and the cytosolic propanoate outlet after the merge.
- opt_pfba   : B0/ALL at the growth optimum (pFBA, 2 LPs).
"""
from common import LP_COUNT, TOL, dump, load, min_total_flux, pfba
from step4_options import apply_base, apply_variant


def route(model, fixed):
    with model:
        model.reactions.xMAINTENANCE.lower_bound = 0.0
        for rid, v in fixed.items():
            model.reactions.get_by_id(rid).bounds = (v, v)
        total, status, flux = min_total_flux(model)
    return {"fixed": fixed, "total": total, "status": status,
            "flux": {k: v for k, v in (flux or {}).items() if abs(v) > TOL}}


model = load("taskA", "gurobi")
out = {"asis_route": route(model, {"R2113": -1.0})}
with model:
    apply_base(model, "B0")
    edits = apply_variant(model, "ALL")
    out["opt_route"] = route(model, {"R25": 1.0, "R2024": -1.0})
    g, s, flux = pfba(model)
    out["opt_pfba"] = {"growth": g, "status": s, "flux": {k: v for k, v in (flux or {}).items() if abs(v) > TOL}}
for k, v in out.items():
    print(k, v.get("status"), v.get("total", v.get("growth")), len(v["flux"]), flush=True)
out["lp_count"] = LP_COUNT["n"]
print("LPs", LP_COUNT["n"], "->", dump(out, "step5_figure_flux.json"))
