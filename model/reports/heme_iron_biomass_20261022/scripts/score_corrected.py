"""Re-score the new model with the two fresh-solve corrections applied (no LPs).

Usage (from the worktree root): python score_corrected.py KO.tsv SUMMARY.json FRESH.json HEM4.json EXPERIMENTAL.csv OUT.json
Corrections (fresh-context solves; the sweep values are kept in ko_compare.tsv):
- YALI1F34051g (HEM3, R409): sweep time_limit -> fresh GLPK optimal 0.0 (Gurobi 0.0).
- YALI1F13039g (HEM4, R792): sweep returned WT growth with R792 above its zero bound (GLPK warm-start
  tolerance artifact) -> fresh GLPK pFBA optimal 0.0 (Gurobi 0.0 at default and 1e-9 tolerance).
"""
import csv
import json
import math
import sys

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
ko_tsv, summary_path, fresh_path, hem4_path, experimental, out = sys.argv[1:7]
summary = json.load(open(summary_path))
fresh = {c["gene_ko"]: c for c in json.load(open(fresh_path))["cases"]
         if c["solver"] == "glpk" and c["gene_ko"] and c["model_sha256"] == summary["new_sha256"][:8]}
hem4 = json.load(open(hem4_path))
corrections = {"YALI1F34051g": ("optimal", fresh["YALI1F34051g"]["growth"], "fresh GLPK slim_optimize"),
               "YALI1F13039g": (hem4["status"], hem4["growth"], "fresh GLPK pFBA")}
assert hem4["model_sha256"] == summary["new_sha256"]
positives = {row["gene_id"] for row in csv.DictReader(open(experimental, encoding="utf-8-sig"))}
mtdna = set(summary["mtdna_out_of_scope"])
wt = summary["wt"]["new"][1]
rows = {}
for row in csv.DictReader(open(ko_tsv), delimiter="\t"):
    rows[row["gene"]] = (row["new_status"], float(row["new_ko"]))
for gene, (status, growth, _) in corrections.items():
    rows[gene] = (status, growth)


def call(status, growth, cutoff):
    if status != "optimal" or math.isnan(growth):
        return "unresolved"
    return "essential" if growth / wt < cutoff else "non-essential"


scores = {}
for cutoff in CUTOFFS:
    counts = {"TP": 0, "FN": 0, "FP": 0, "TN": 0, "unresolved": 0}
    for gene, (status, growth) in rows.items():
        if gene in mtdna:
            continue
        c = call(status, growth, cutoff)
        key = "unresolved" if c == "unresolved" else (
            ("TP" if gene in positives else "FP") if c == "essential" else ("FN" if gene in positives else "TN"))
        counts[key] += 1
    tp, fn, fp, tn = (counts[k] for k in ("TP", "FN", "FP", "TN"))
    d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    scores[str(cutoff)] = {**counts, "recall": tp / (tp + fn), "specificity": tn / (tn + fp),
                           "precision": tp / (tp + fp), "mcc": (tp * tn - fp * fn) / d}
result = {"new_sha256": summary["new_sha256"], "corrections": {g: {"status": s, "growth": v, "source": src}
                                                               for g, (s, v, src) in corrections.items()},
          "scores_new_corrected": scores, "scores_baseline": summary["scores"]["baseline"],
          "negative_class": summary["negative_class"], "validation_status": summary["validation_status"]}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for c in CUTOFFS:
    b, n = summary["scores"]["baseline"][str(c)], scores[str(c)]
    print(c, "baseline", {k: b[k] for k in ("TP", "FN", "FP", "TN", "unresolved")}, "mcc %.4f" % b["mcc"],
          "| new corrected", {k: n[k] for k in ("TP", "FN", "FP", "TN", "unresolved")}, "mcc %.4f" % n["mcc"])
