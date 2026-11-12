"""Acceptance checks 4 and 5 on the built candidate (TASK.md section 8): WT, no-uptake ATP maximum, and the sweep
comparison with the Stage 1 in-memory V1 screen. No model write.

Run from platform/ with PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts
Args: CANDIDATE.xml MEDIUM.csv PROFILE.json SWEEP_BUILT.tsv SWEEP_V1.tsv OUT.json
"""
import csv
import json
import sys
import warnings

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, medium, profile, built_tsv, v1_tsv, out_path = sys.argv[1:7]
V1_WT = 0.3120526953810365
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
configure(model)
wt = solve(model, matrix(model), ())
with model:
    positive = [(r.id, r.lower_bound) for r in model.reactions if r.lower_bound > 0]
    for r in model.exchanges:
        r.bounds = (0, 1000)
    for rid, _ in positive:
        model.reactions.get_by_id(rid).lower_bound = 0
    model.reactions.get_by_id("R_NTP1").lower_bound = 0
    model.objective = "R_NTP1"
    model.solver.problem.reset()
    atp = model.slim_optimize(error_value=float("nan"))
    atp_status = model.solver.status
cutoffs = ("0.01", "0.05", "0.1", "0.15")
built = {r["gene"]: r for r in csv.DictReader(open(built_tsv), delimiter="\t")}
v1 = {r["gene"]: r for r in csv.DictReader(open(v1_tsv), delimiter="\t")}
class_diff, ko_diff = [], []
for g in sorted(set(built) | set(v1)):
    if g not in built or g not in v1:
        class_diff.append({"gene": g, "missing_in": "built" if g not in built else "v1"})
        continue
    for rule in ("A", "B"):
        for c in cutoffs:
            if built[g][f"class_{rule}_{c}"] != v1[g][f"class_{rule}_{c}"]:
                class_diff.append({"gene": g, "column": f"class_{rule}_{c}", "built": built[g][f"class_{rule}_{c}"],
                                   "v1": v1[g][f"class_{rule}_{c}"]})
    a, b = built[g]["ko"], v1[g]["ko"]
    if (a in ("", "None")) != (b in ("", "None")) or (a not in ("", "None") and abs(float(a) - float(b)) > 1e-9):
        ko_diff.append({"gene": g, "built": a, "v1": b})
result = {
    "model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
    "strain_profile_sha256": ctx.strain_profile_sha256, "r2121_bounds": list(model.reactions.R2121.bounds),
    "wt": wt, "wt_minus_v1": None if wt["growth"] is None else wt["growth"] - V1_WT,
    "criterion_4_wt_within_1e-9": wt["passed"] and abs(wt["growth"] - V1_WT) <= 1e-9,
    "energy": {"status": atp_status, "max_atp_hydrolysis_no_uptake": None if atp_status != "optimal" else atp,
               "relaxed_lower_bounds": positive},
    "criterion_4_energy_zero": atp_status == "optimal" and abs(atp) <= 1e-9,
    "sweep_genes": {"built": len(built), "v1": len(v1)},
    "class_differences": class_diff, "ko_differences_over_1e-9": ko_diff,
    "criterion_5_sweep_reproduces_v1": not class_diff and not ko_diff,
    "lp_count": 2,
}
json.dump(result, open(out_path, "w"), indent=1)
print(json.dumps({k: v for k, v in result.items() if k not in ("energy", "wt")}, indent=1))
print("wt", wt["growth"], wt["passed"], "| atp", atp_status, atp)
