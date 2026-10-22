"""Single-gene knockouts on two models under one condition, compared gene by gene (read-only).

Usage (from platform/): python ../model/reports/tpp_requirement_20261022/scripts/ko_compare.py \
    BASELINE.xml NEW.xml MEDIUM.csv PROFILE.json OUT.tsv OUT_SUMMARY.json
Condition: Task A medium + PO1f, GLPK, objective unchanged. This is a model-vs-model comparison, not a
benchmark score: no experimental labels are used. Class: essential when KO/WT is strictly below
the cutoff (model/expected/benchmark.md).
"""
import csv
import json
import math
import sys
import warnings

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
CHANGE_TOLERANCE = 1e-6

baseline_path, new_path, media_path, profile_path, out_tsv, out_json = sys.argv[1:7]


def knockouts(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(
            model_path=path, media_path=media_path, strain_profile_path=profile_path)
    model = ctx.model
    model.solver = "glpk"
    wt = model.slim_optimize(error_value=float("nan"))
    wt_status = model.solver.status
    rows = {}
    for gene in sorted(g.id for g in model.genes if g.id not in ctx.excluded_runtime_genes):
        with model:
            model.genes.get_by_id(gene).knock_out()
            growth = model.slim_optimize(error_value=float("nan"))
            rows[gene] = (model.solver.status, growth)
    return ctx, wt, wt_status, rows


def classes(wt, growth, status):
    if status != "optimal" or math.isnan(growth) or not wt:
        return ["unresolved"] * len(CUTOFFS)
    return ["essential" if growth / wt < cutoff else "non-essential" for cutoff in CUTOFFS]


base_ctx, base_wt, base_wt_status, base = knockouts(baseline_path)
new_ctx, new_wt, new_wt_status, new = knockouts(new_path)
if set(base) != set(new):
    raise ValueError("The two models have different gene sets")
changed_growth, changed_class = [], {str(c): [] for c in CUTOFFS}
with open(out_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene", "base_status", "base_ko", "base_ratio", "new_status", "new_ko", "new_ratio"]
                    + [f"base_class_{c}" for c in CUTOFFS] + [f"new_class_{c}" for c in CUTOFFS])
    for gene in sorted(base):
        (bs, bg), (ns, ng) = base[gene], new[gene]
        bc, nc = classes(base_wt, bg, bs), classes(new_wt, ng, ns)
        writer.writerow([gene, bs, bg, bg / base_wt, ns, ng, ng / new_wt] + bc + nc)
        if bs != ns or (math.isnan(bg) != math.isnan(ng)) or (not math.isnan(bg) and abs(bg - ng) > CHANGE_TOLERANCE):
            changed_growth.append({"gene": gene, "base": [bs, bg], "new": [ns, ng]})
        for cutoff, b, n in zip(CUTOFFS, bc, nc):
            if b != n:
                changed_class[str(cutoff)].append({"gene": gene, "base": b, "new": n})
summary = {
    "condition": "Task A medium + PO1f, GLPK",
    "baseline_sha256": base_ctx.canonical_model_sha256, "new_sha256": new_ctx.canonical_model_sha256,
    "medium_sha256": base_ctx.medium_sha256, "strain_profile_sha256": base_ctx.strain_profile_sha256,
    "excluded_runtime_genes": list(base_ctx.excluded_runtime_genes),
    "wt": {"baseline": [base_wt_status, base_wt], "new": [new_wt_status, new_wt]},
    "genes": len(base),
    "non_optimal": {"baseline": sorted(g for g, (s, _) in base.items() if s != "optimal"),
                    "new": sorted(g for g, (s, _) in new.items() if s != "optimal")},
    "essential_counts": {str(c): {"baseline": sum(classes(base_wt, g, s)[i] == "essential" for s, g in base.values()),
                                  "new": sum(classes(new_wt, g, s)[i] == "essential" for s, g in new.values())}
                         for i, c in enumerate(CUTOFFS)},
    "growth_changed_over_1e-6": changed_growth,
    "class_changed": changed_class,
    "lp_count": 2 * (len(base) + 1),
}
with open(out_json, "w") as handle:
    json.dump(summary, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: v for k, v in summary.items() if k != "excluded_runtime_genes"}, indent=1))
