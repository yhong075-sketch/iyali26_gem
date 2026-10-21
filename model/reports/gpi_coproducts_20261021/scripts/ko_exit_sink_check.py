"""Is the open GPI exit the cause of the +8.5e-7 knockout gain? (read-only, 3 LPs)

Usage (from platform/): PYTHONPATH=$PWD python ko_exit_sink_check.py BASELINE.xml NEW.xml RESEARCH_ROOT OUT.json
Knocks out YALI1E12730g (sole gene of R304) on both models, and on the new model with R_DM_GPI closed.
"""
import hashlib
import json
import sys
import warnings

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

GENE = "YALI1E12730g"
baseline_path, new_path, research, out = sys.argv[1:5]


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(
            model_path=path, media_path=f"{research}/state/media/sd_leu.csv",
            strain_profile_path=f"{research}/state/strain_profiles/po1f_sd_leu.json")
    ctx.model.solver = "glpk"
    return ctx.model


def ko(model, close_exit=False):
    with model:
        model.genes.get_by_id(GENE).knock_out()
        if close_exit:
            model.reactions.R_DM_GPI.bounds = (0.0, 0.0)
        growth = model.slim_optimize(error_value=float("nan"))
        fluxes = {rid: model.reactions.get_by_id(rid).flux for rid in ("R_DM_GPI", "R_DAGt_er_em", "R373", "R1773")}
        return {"status": model.solver.status, "growth": growth, "fluxes": fluxes,
                "gene_reactions": sorted(r.id for r in model.genes.get_by_id(GENE).reactions)}


result = {"gene": GENE,
          "baseline_sha256": hashlib.sha256(open(baseline_path, "rb").read()).hexdigest(),
          "new_sha256": hashlib.sha256(open(new_path, "rb").read()).hexdigest()}
base, new = context(baseline_path), context(new_path)
result["baseline_ko"] = ko(base)
result["new_ko"] = ko(new)
result["new_ko_exit_closed"] = ko(new, close_exit=True)
result["lp_count"] = 3
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps(result, indent=1))
