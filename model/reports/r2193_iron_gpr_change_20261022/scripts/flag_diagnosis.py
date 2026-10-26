"""Size and location of the check failures in the stopped first sweep (7 LPs, read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/r2193_iron_gpr_change_20261022/scripts python \
    ../model/reports/r2193_iron_gpr_change_20261022/scripts/flag_diagnosis.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json GENE...
"""
import json
import sys
import warnings

import numpy as np

from checked_solve import matrix, reset_basis
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, medium, profile, out, *genes = sys.argv[1:]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
model.solver = "glpk"
stoichiometry = matrix(model)
rows = {}
for gene in genes:
    with model:
        model.genes.get_by_id(gene).knock_out()
        reset_basis(model)
        growth = model.slim_optimize(error_value=float("nan"))
        primal = model.solver.primal_values
        worst_bound, flux = (0.0, None), []
        for r in model.reactions:
            f, b = primal[r.forward_variable.name], primal[r.reverse_variable.name]
            v = f - b
            flux.append(v)
            gap = max(r.bounds[0] - v, v - r.bounds[1], -f, -b)
            if gap > worst_bound[0]:
                worst_bound = (gap, {"reaction": r.id, "flux": v, "bounds": list(r.bounds), "forward": f, "reverse": b,
                                     "gpr": r.gene_reaction_rule})
        residual = stoichiometry @ np.array(flux)
        i = int(np.argmax(np.abs(residual)))
        rows[gene] = {"status": model.solver.status, "growth": growth, "ko_over_wt": growth / 0.3124267446456700,
                      "knocked_reactions": sorted(r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)),
                      "max_bound_violation": worst_bound[0], "worst_bound": worst_bound[1],
                      "max_residual": float(abs(residual[i])), "worst_row": model.metabolites[i].id}
result = {"model_sha256": ctx.canonical_model_sha256, "genes": rows, "lp_count": len(genes)}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for g, r in rows.items():
    wb = r["worst_bound"] or {}
    print(f"{g} {r['status']} KO/WT={r['ko_over_wt']:.4f} bound_viol={r['max_bound_violation']:.2e} on {wb.get('reaction')} "
          f"flux={wb.get('flux')} bounds={wb.get('bounds')} fwd={wb.get('forward')} rev={wb.get('reverse')} | residual={r['max_residual']:.2e} on {r['worst_row']}")
