"""Why the YALI1D22124g knockout is 'optimal, growth 0' on the baseline but 'infeasible' on the new model.

Usage (from platform/): python ko_d22124g_diagnosis.py BASELINE.xml NEW.xml RESEARCH_ROOT OUT.json
Checks, in order (GLPK, SD-Leu + PO1f):
1. the KO on each model in a fresh context (no warm start from earlier knockouts);
2. the KO on the new model with the five added reactions closed;
3. whether the baseline KO flux vector, extended with zero flux on the five added reactions, satisfies
   every mass-balance row and bound of the new model (no solve needed: pure arithmetic).
"""
import json
import math
import sys
import warnings

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

GENE = "YALI1D22124g"
LINKS = ["R_GPIGNT", "R_PEt_em_er", "R_DAGt_er_em", "R_ACt_er_cy", "R_DM_GPI"]
baseline_path, new_path, research, out = sys.argv[1:5]


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(
            model_path=path, media_path=f"{research}/state/media/sd_leu.csv",
            strain_profile_path=f"{research}/state/strain_profiles/po1f_sd_leu.json")
    ctx.model.solver = "glpk"
    return ctx.model


def ko(model, close=()):
    with model:
        model.genes.get_by_id(GENE).knock_out()
        for rid in close:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        value = model.slim_optimize(error_value=float("nan"))
        status = model.solver.status
        fluxes = {r.id: r.flux for r in model.reactions} if status == "optimal" else None
        closed = sorted(r.id for r in model.genes.get_by_id(GENE).reactions)
    return {"status": status, "growth": None if math.isnan(value) else value}, fluxes, closed


base, new = context(baseline_path), context(new_path)
result = {"gene": GENE}
result["baseline_ko"], base_fluxes, closed = ko(base)
result["gene_reactions"] = closed
result["new_ko"], _, _ = ko(new)
result["new_ko_links_closed"], _, _ = ko(new, LINKS)

# 3. Feasibility of the baseline KO point in the new model.
point = dict(base_fluxes)
point.update({rid: 0.0 for rid in LINKS})
worst_row, worst_bound = (0.0, None), (0.0, None)
for metabolite in new.metabolites:
    residual = sum(c * point[r.id] for r in metabolite.reactions for c in [r.metabolites[metabolite]])
    if abs(residual) > abs(worst_row[0]):
        worst_row = (residual, metabolite.id)
with new:
    new.genes.get_by_id(GENE).knock_out()
    for reaction in new.reactions:
        v = point[reaction.id]
        excess = max(reaction.lower_bound - v, v - reaction.upper_bound, 0.0)
        if excess > abs(worst_bound[0]):
            worst_bound = (excess, reaction.id)
result["baseline_point_in_new_model"] = {"max_abs_row_residual": worst_row, "max_bound_violation": worst_bound,
                                         "objective_at_point": point["biomass_C"]}
result["lp_count"] = 3
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps(result, indent=1))
