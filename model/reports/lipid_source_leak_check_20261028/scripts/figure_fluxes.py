"""pFBA points for the leak figure on Task A + PO1f (4 LPs, inside the 300-LP cap).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/lipid_source_leak_check_20261028/scripts/figure_fluxes.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
F1 WT; F2 INO1 knockout (YALI1B05798g) with inositol uptake R1263 and R950 closed. F3, the inositol-demand maximum,
is the pFBA route already stored by source_check.py (C1|base|inositol|m814[C_cy]).
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, medium, profile, out = sys.argv[1:5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile).model
configure(model)
res = {}
model.solver.problem.reset()
res["F1"] = {"growth": model.slim_optimize(), "fluxes": {k: float(v) for k, v in pfba(model).fluxes.items() if abs(v) > 1e-12}}
with model:
    g = model.genes.get_by_id("YALI1B05798g")
    g.knock_out()
    for rid in ("R1263", "R950"):
        model.reactions.get_by_id(rid).lower_bound = 0.0
    res["F2_check"] = solve(model, matrix(model), {r.id for r in g.reactions})
    model.solver.problem.reset()
    res["F2"] = {"growth": res["F2_check"]["growth"], "fluxes": {k: float(v) for k, v in pfba(model).fluxes.items() if abs(v) > 1e-12}}
res["lp_count"] = 6  # F1 optimize + pFBA (2); F2 checked solve + pFBA (2)
json.dump(res, open(out, "w"))
print("F1", res["F1"]["growth"], "F2", res["F2_check"], file=sys.stderr)
