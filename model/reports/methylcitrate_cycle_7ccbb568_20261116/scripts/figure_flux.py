"""Flux for the schematic on the built models (read only; no model file is written).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python figure_flux.py CONTROL.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json

Task A + PO1f, Gurobi 1e-9 (checked_solve.configure), solution and basis discarded before every LP.
- pFBA at the growth optimum on the control (7ccbb568) and on the candidate.
- Route on the candidate, as in the investigation's step 5b: R25 fixed at 1, xMAINTENANCE lower bound 0, biomass_C
  free, minimum total flux (sum of forward and reverse variables).
"""
import json
import sys
import warnings
from pathlib import Path

from cobra.flux_analysis import pfba
from optlang.symbolics import Zero

from checked_solve import configure, reset
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

control, candidate, medium, profile, out_json = sys.argv[1:6]
LPS = {"n": 0}


def load(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=Path(path), media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


def run_pfba(model):
    reset(model)
    LPS["n"] += 2
    sol = pfba(model)
    return {"status": sol.status, "growth": float(sol.fluxes["biomass_C"]),
            "flux": {k: float(v) for k, v in sol.fluxes.items() if abs(v) > 1e-9}}


def route(model):
    with model:
        model.reactions.R25.bounds = (1.0, 1.0)
        model.reactions.xMAINTENANCE.lower_bound = 0.0
        coefficients = {}
        for r in model.reactions:
            coefficients[r.forward_variable] = 1.0
            coefficients[r.reverse_variable] = 1.0
        model.objective = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective.set_linear_coefficients(coefficients)
        reset(model)
        LPS["n"] += 1
        model.slim_optimize(error_value=float("nan"))
        status = model.solver.status
        if status != "optimal":
            return {"status": status, "nonzero": None}
        primal = model.solver.primal_values
        flux = {r.id: primal[r.forward_variable.name] - primal[r.reverse_variable.name] for r in model.reactions}
        return {"status": status, "fixed_at": 1.0, "nonzero": {k: v for k, v in flux.items() if abs(v) > 1e-9}}


out = {}
c = load(control)
out["control"] = {"model_sha256": c.canonical_model_sha256, "pfba": run_pfba(c.model)}
k = load(candidate)
out["candidate"] = {"model_sha256": k.canonical_model_sha256, "pfba": run_pfba(k.model), "route_R25": route(k.model)}
out["lp_count"] = LPS["n"]
Path(out_json).write_text(json.dumps(out, indent=1))
for tag in ("control", "candidate"):
    p = out[tag]["pfba"]
    print(tag, out[tag]["model_sha256"][:8], p["status"], repr(p["growth"]),
          {r: p["flux"].get(r, 0.0) for r in ("R25", "R834", "R2025", "R2113", "R_OBDHm", "R1172")})
r = out["candidate"]["route_R25"]
print("route", r["status"], {x: r["nonzero"].get(x, 0.0) for x in ("R25", "R_OBDHm", "R834", "R2025", "R2113", "R2077")}
      if r["nonzero"] else None)
print("LPs", LPS["n"])
