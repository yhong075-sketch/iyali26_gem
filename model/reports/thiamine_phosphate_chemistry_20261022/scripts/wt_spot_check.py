"""WT growth on the baseline and the new model, one LP each (read-only).

Usage (from platform/): PYTHONPATH=$PWD python ../model/reports/thiamine_phosphate_chemistry_20261022/scripts/wt_spot_check.py \
    BASELINE.xml NEW.xml MEDIUM.csv PROFILE.json OUT.json
Condition: Task A medium + PO1f, GLPK. Formulas, charges and annotations are not in the LP, so the two values
are expected to be identical.
"""
import json
import sys
import warnings

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

baseline_path, new_path, media_path, profile_path, out = sys.argv[1:6]
result = {"condition": "Task A medium + PO1f, GLPK"}
for name, path in (("baseline", baseline_path), ("new", new_path)):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=media_path,
                                                strain_profile_path=profile_path)
    ctx.model.solver = "glpk"
    growth = ctx.model.slim_optimize(error_value=float("nan"))
    result[name] = {"model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
                    "strain_profile_sha256": ctx.strain_profile_sha256,
                    "status": ctx.model.solver.status, "growth": growth}
result["identical"] = result["baseline"]["growth"] == result["new"]["growth"]
result["lp_count"] = 2
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(json.dumps(result, indent=1))
