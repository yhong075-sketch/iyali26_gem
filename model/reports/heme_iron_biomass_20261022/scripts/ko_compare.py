"""Single-gene knockouts on two models under one condition, compared gene by gene and scored (read-only).

Usage (from platform/): python ../model/reports/heme_iron_biomass_20261022/scripts/ko_compare.py \
    BASELINE.xml NEW.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.tsv OUT_SUMMARY.json
Adapted from ../gpi_chain_links_20261021/scripts/ko_compare.py: the condition is now passed in (Task A
medium + PO1f), and both models are also scored against the experimental list per
model/expected/benchmark.md. GLPK, objective unchanged. Class: essential when KO/WT is strictly below
the cutoff. Non-optimal solves are unresolved, never essential. Negatives are the user-defined class
(model genes not on the list); mitochondrially encoded genes are out of scope. The list was used in
model development, so scores are a development reference, not independent validation.
"""
import csv
import json
import math
import sys
import time
import warnings
from pathlib import Path

from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL
from tools.dfba_new_fn_essentiality import load_experimental_reference

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
CHANGE_TOLERANCE = 1e-6
# The first run (2026-10-22) hung for 34 min inside one GLPK simplex solve; a timed-out LP is
# recorded with its solver status (unresolved), never as essential.
LP_TIMEOUT_SECONDS = 60

baseline_path, new_path, medium, profile, experimental, out_tsv, out_json = sys.argv[1:8]
positives, reference_metadata = load_experimental_reference(Path(experimental))
positives = set(positives)
mtdna = set(json.loads((MODEL.root / "expected/mtdna_genes_out_of_scope.json").read_text())["genes"])


def knockouts(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    model = ctx.model
    model.solver = "glpk"
    model.solver.configuration.timeout = LP_TIMEOUT_SECONDS
    wt = model.slim_optimize(error_value=float("nan"))
    wt_status = model.solver.status
    rows = {}
    for index, gene in enumerate(sorted(g.id for g in model.genes if g.id not in ctx.excluded_runtime_genes), 1):
        started = time.monotonic()
        with model:
            model.genes.get_by_id(gene).knock_out()
            growth = model.slim_optimize(error_value=float("nan"))
            rows[gene] = (model.solver.status, growth)
        print(f"{ctx.canonical_model_sha256[:8]} {index} {gene} {rows[gene][0]} {time.monotonic() - started:.2f}s",
              file=sys.stderr, flush=True)
    return ctx, wt, wt_status, rows


def classes(wt, growth, status):
    if status != "optimal" or math.isnan(growth) or not wt:
        return ["unresolved"] * len(CUTOFFS)
    return ["essential" if growth / wt < cutoff else "non-essential" for cutoff in CUTOFFS]


def score(wt, rows):
    result = {}
    for i, cutoff in enumerate(CUTOFFS):
        counts = {"TP": 0, "FN": 0, "FP": 0, "TN": 0, "unresolved": 0}
        for gene, (status, growth) in rows.items():
            if gene in mtdna:
                continue
            call = classes(wt, growth, status)[i]
            if call == "unresolved":
                counts["unresolved"] += 1
            elif call == "essential":
                counts["TP" if gene in positives else "FP"] += 1
            else:
                counts["FN" if gene in positives else "TN"] += 1
        tp, fn, fp, tn = (counts[k] for k in ("TP", "FN", "FP", "TN"))
        denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        result[str(cutoff)] = {**counts,
                               "recall": tp / (tp + fn) if tp + fn else None,
                               "specificity": tn / (tn + fp) if tn + fp else None,
                               "precision": tp / (tp + fp) if tp + fp else None,
                               "mcc": (tp * tn - fp * fn) / denominator if denominator else None}
    return result


base_ctx, base_wt, base_wt_status, base = knockouts(baseline_path)
new_ctx, new_wt, new_wt_status, new = knockouts(new_path)
if set(base) != set(new):
    raise ValueError("The two models have different gene sets")
changed_growth, changed_class = [], {str(c): [] for c in CUTOFFS}
with open(out_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene", "experimental_positive", "mtdna_out_of_scope", "base_status", "base_ko", "base_ratio",
                     "new_status", "new_ko", "new_ratio"]
                    + [f"base_class_{c}" for c in CUTOFFS] + [f"new_class_{c}" for c in CUTOFFS])
    for gene in sorted(base):
        (bs, bg), (ns, ng) = base[gene], new[gene]
        bc, nc = classes(base_wt, bg, bs), classes(new_wt, ng, ns)
        writer.writerow([gene, gene in positives, gene in mtdna, bs, bg, bg / base_wt, ns, ng, ng / new_wt] + bc + nc)
        if bs != ns or (math.isnan(bg) != math.isnan(ng)) or (not math.isnan(bg) and abs(bg - ng) > CHANGE_TOLERANCE):
            changed_growth.append({"gene": gene, "base": [bs, bg], "new": [ns, ng]})
        for cutoff, b, n in zip(CUTOFFS, bc, nc):
            if b != n:
                changed_class[str(cutoff)].append({"gene": gene, "experimental_positive": gene in positives,
                                                   "mtdna_out_of_scope": gene in mtdna, "base": b, "new": n})
summary = {
    "condition": "Task A medium + PO1f, GLPK",
    "lp_timeout_seconds": LP_TIMEOUT_SECONDS,
    "baseline_sha256": base_ctx.canonical_model_sha256, "new_sha256": new_ctx.canonical_model_sha256,
    "medium_sha256": base_ctx.medium_sha256, "strain_profile_sha256": base_ctx.strain_profile_sha256,
    "experimental_reference": reference_metadata,
    "negative_class": "user-defined: model genes not on the consensus list",
    "validation_status": "development reference, not independent validation",
    "mtdna_out_of_scope": sorted(mtdna),
    "excluded_runtime_genes": list(base_ctx.excluded_runtime_genes),
    "wt": {"baseline": [base_wt_status, base_wt], "new": [new_wt_status, new_wt]},
    "genes": len(base),
    "in_model_positives": len(positives & set(base)),
    "out_of_scope_positives": len(positives - set(base)),
    "non_optimal": {"baseline": sorted(g for g, (s, _) in base.items() if s != "optimal"),
                    "new": sorted(g for g, (s, _) in new.items() if s != "optimal")},
    "essential_counts": {str(c): {"baseline": sum(classes(base_wt, g, s)[i] == "essential" for s, g in base.values()),
                                  "new": sum(classes(new_wt, g, s)[i] == "essential" for s, g in new.values())}
                         for i, c in enumerate(CUTOFFS)},
    "scores": {"baseline": score(base_wt, base), "new": score(new_wt, new)},
    "growth_changed_over_1e-6": changed_growth,
    "class_changed": changed_class,
    "lp_count": 2 * (len(base) + 1),
}
with open(out_json, "w") as handle:
    json.dump(summary, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: v for k, v in summary.items()
                  if k not in ("excluded_runtime_genes", "mtdna_out_of_scope", "experimental_reference")}, indent=1))
