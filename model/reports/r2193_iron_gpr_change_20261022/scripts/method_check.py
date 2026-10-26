"""Method check on the known artifact (HEM4, YALI1F13039g, on 7f836c3a): 4 LPs, read-only.

Rerun after the audit with the targeted check and the knocked-out reactions supplied
(method_check_targeted.json); method_check.json is the first-version output.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/r2193_iron_gpr_change_20261022/scripts python \
    ../model/reports/r2193_iron_gpr_change_20261022/scripts/method_check.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
Case 1: WT, then the knockout without a basis reset (expected: the stale WT point, flagged by check()).
Case 2: WT, then the knockout after reset_basis() (expected: 0, passing check()).
"""
import json
import sys
import warnings

from checked_solve import matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, medium, profile, out = sys.argv[1:5]
result = {"gene": "YALI1F13039g", "cases": {}}
for label, reset in (("warm_start_no_reset", False), ("basis_reset", True)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    model.solver = "glpk"
    stoichiometry = matrix(model)
    wt = solve(model, stoichiometry, (), reset=False)
    with model:
        model.genes.YALI1F13039g.knock_out()
        knocked = {r.id for r in model.genes.YALI1F13039g.reactions if r.bounds == (0, 0)}
        ko = solve(model, stoichiometry, knocked, reset=reset)
        ko["knocked_reactions"] = sorted(knocked)
    result["cases"][label] = {"model_sha256": ctx.canonical_model_sha256, "wt": wt, "knockout": ko}
result["lp_count"] = 4
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps(result, indent=1))
