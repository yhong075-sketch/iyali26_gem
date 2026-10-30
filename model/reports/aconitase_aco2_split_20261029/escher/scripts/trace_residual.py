"""Where does the residual growth come from when YALI1D11984g is knocked out on the candidate? (read-only, 3 LPs)

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    python ../model/reports/aconitase_aco2_split_20261029/escher/scripts/trace_residual.py MEDIUM.csv PROFILE.json \
    CANDIDATE.xml OUT.json
One checked solve and one pFBA (Task A medium + PO1f, Gurobi 1e-9). Prints every medium uptake with its flux and its
ratio to the uptake bound, and every producer and consumer of 2-oxoglutarate and L-glutamate in all compartments.
pFBA is one optimum; the route is not claimed to be unique.
"""
import json
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

KO_GENE = "YALI1D11984g"
medium, profile, cand_path, out = sys.argv[1:5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=cand_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
configure(model)
result = {"candidate_sha256": ctx.canonical_model_sha256, "ko_gene": KO_GENE}
with model:
    model.genes.get_by_id(KO_GENE).knock_out()
    knocked = tuple(r.id for r in model.genes.get_by_id(KO_GENE).reactions if r.bounds == (0, 0))
    checked = solve(model, matrix(model), knocked)
    if not checked["passed"]:
        sys.exit(f"STOP: {checked}")
    sol = pfba(model)
    v = sol.fluxes
    result["checked"], result["closed"] = checked, list(knocked)
    uptakes = []
    for r in model.exchanges:
        if v[r.id] < -1e-9:
            uptakes.append({"id": r.id, "name": r.name, "flux": float(v[r.id]), "lb": r.lower_bound,
                            "at_bound": bool(abs(v[r.id] - r.lower_bound) <= 1e-7)})
    result["uptakes"] = sorted(uptakes, key=lambda u: u["flux"])
    nodes = {}
    for m in model.metabolites:
        if m.name.split("_")[0] in ("2-oxoglutarate", "L-glutamate"):
            rows = []
            for r in m.reactions:
                net = float(v[r.id] * r.metabolites[m])
                if abs(net) > 1e-9:
                    rows.append({"id": r.id, "net": net, "equation": r.build_reaction_string(use_metabolite_names=True)})
            if rows:
                nodes[m.id] = {"name": m.name, "rows": sorted(rows, key=lambda x: -x["net"])}
    result["nodes"] = nodes
json.dump(result, open(out, "w"), indent=1)
print("growth", checked["growth"])
for u in result["uptakes"]:
    print(f"uptake {u['id']:14} {u['flux']:10.4g} lb {u['lb']:10.4g} at_bound {u['at_bound']}  {u['name']}")
for mid, n in nodes.items():
    print(mid, n["name"])
    for row in n["rows"]:
        print(f"    {row['id']:12} {row['net']:10.4g}  {row['equation'][:120]}")
