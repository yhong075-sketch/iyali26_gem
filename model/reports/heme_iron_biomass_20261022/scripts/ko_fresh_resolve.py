"""Fresh re-solve of the one knockout that hit the GLPK time limit (YALI1F34051g on the new model).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/heme_iron_biomass_20261022/scripts/ko_fresh_resolve.py \
    BASELINE.xml NEW.xml MEDIUM.csv PROFILE.json OUT.json
Each case loads its own context. GLPK as in the sweep (60 s limit); Gurobi only as a labelled
cross-check of the same LP. Also records R792 (HEM4) GPR and the R409 reaction closure.
"""
import json
import math
import sys
import time
import warnings

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

GENE, REACTION = "YALI1F34051g", "R409"
baseline_path, new_path, medium, profile, out = sys.argv[1:6]
lp_count = 0


def run(path, solver, gene=None, reaction=None):
    global lp_count
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    model.solver = solver
    model.solver.configuration.timeout = 60
    if gene:
        model.genes.get_by_id(gene).knock_out()
    if reaction:
        model.reactions.get_by_id(reaction).bounds = (0.0, 0.0)
    started = time.monotonic()
    value = model.slim_optimize(error_value=float("nan"))
    lp_count += 1
    return {"model_sha256": ctx.canonical_model_sha256[:8], "solver": solver, "gene_ko": gene, "reaction_closed": reaction,
            "status": model.solver.status, "growth": None if math.isnan(value) else value,
            "seconds": round(time.monotonic() - started, 3)}


result = {"condition": "Task A medium + PO1f", "cases": [
    run(new_path, "glpk", gene=GENE),
    run(new_path, "glpk", reaction=REACTION),
    run(new_path, "gurobi", gene=GENE),
    run(baseline_path, "gurobi", gene=GENE),
    run(new_path, "gurobi"),
]}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=new_path, media_path=medium, strain_profile_path=profile)
r792 = ctx.model.reactions.R792
result["R792"] = {"name": r792.name, "gpr": r792.gene_reaction_rule, "bounds": list(r792.bounds)}
result["lp_count"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps(result, indent=1))
