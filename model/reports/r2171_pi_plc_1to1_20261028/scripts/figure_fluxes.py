"""pFBA points for the before / after Escher maps on Task A + PO1f without inositol (read-only, 7 LPs).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/r2171_pi_plc_1to1_20261028/scripts/figure_fluxes.py BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Inositol uptake R1263 and R950 closed in every point. B1 baseline, INO1 intact; A1 candidate, INO1 intact (pFBA each);
A2 candidate with the INO1 knockout (checked growth only; Stage 1 found WT growth here on the baseline).
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

baseline, candidate, medium, profile, out = sys.argv[1:6]
res = {}
for label, path in (("B1", baseline), ("A1", candidate)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile).model
    configure(model)
    with model:
        for rid in ("R1263", "R950"):
            model.reactions.get_by_id(rid).lower_bound = 0.0
        check = solve(model, matrix(model), set())
        model.solver.problem.reset()
        res[label] = {"check": check, "fluxes": {k: float(v) for k, v in pfba(model).fluxes.items() if abs(v) > 1e-12}}
        if label == "A1":
            gene = model.genes.get_by_id("YALI1B05798g")
            gene.knock_out()
            res["A2_check"] = solve(model, matrix(model), {r.id for r in gene.reactions})
res["lp_count"] = 7  # B1 and A1: checked solve + pFBA (2) each; A2: checked solve
json.dump(res, open(out, "w"))
for k in ("B1", "A1"):
    f = res[k]["fluxes"]
    print(k, res[k]["check"]["growth"], {r: f.get(r, 0.0) for r in ("R2171", "R559", "R558", "R862", "R1839", "R2168")}, file=sys.stderr)
print("A2", res["A2_check"], file=sys.stderr)
