"""Step 6, after the sweep: compare the option sweep (B0/ALL) with the as-is baseline, gene by gene. No LP.

Usage: python step6_compare.py [CONTROL.tsv]
  baseline  sweep684_candidate.tsv (353cefdc; step6_identity.py shows it is the same problem as 7ccbb568)
  option    $ART/out/sweep_B0ALL.tsv
  control   optional: the as-is (B0/K) re-run of the changed genes in the step-6 harness; its rows must equal the
            baseline's
Writes tables/sweep_class_changes.csv, tables/sweep_scores.csv and tables/step6_compare.json.
Ratios use the unrounded KO/WT of each sweep's own WT; Rule B scores a certified infeasible knockout as 0.
"""
import csv
import json
import sys
from pathlib import Path

from common import ART, KO_GENES, LIST, RR, sha256

HERE = Path(__file__).resolve().parents[1]
FOLINIC = HERE.parent / "folinic_acid_module_20261112"
BASE_TSV = RR / "artifacts/tasks/folinic_acid_module_20261112/sweep684_candidate.tsv"
BASE_SHA = "49b4f27be05beffe32b206b75e7fe571b69d07af6d13a48b568f8e39f07e91db"
BASE_SUM = FOLINIC / "sweep684_candidate_summary.json"
OPT_TSV, OPT_SUM = ART / "out/sweep_B0ALL.tsv", HERE / "tables/sweep_B0ALL_summary.json"
CUTS = ("0.01", "0.05", "0.1", "0.15")

if sha256(BASE_TSV) != BASE_SHA:
    sys.exit("baseline TSV hash mismatch")


def rows(path):
    with open(path) as f:
        return {r["gene"]: r for r in csv.DictReader(f, delimiter="\t")}


def growth_b(r):
    """Final growth under Rule B: a passed point's growth, 0 if certified, else None."""
    if r["ko"] not in ("", "None"):
        return float(r["ko"])
    return 0.0 if r["cert_certified"] == "True" else None


base, opt = rows(BASE_TSV), rows(OPT_TSV)
sb, so = json.load(open(BASE_SUM)), json.load(open(OPT_SUM))
wt_b, wt_o = sb["wt"]["growth"], so["wt"]["growth"]
if set(base) != set(opt):
    sys.exit(f"gene sets differ: {sorted(set(base) ^ set(opt))}")

changes, ratio_moves = [], []
for g in sorted(base):
    b, o = base[g], opt[g]
    gb, go = growth_b(b), growth_b(o)
    rb = None if gb is None else gb / wt_b
    ro = None if go is None else go / wt_o
    cls = {f"{rule}_{c}": (b[f"class_{rule}_{c}"], o[f"class_{rule}_{c}"]) for rule in "AB" for c in CUTS}
    if any(x != y for x, y in cls.values()):
        row = {"gene": g, "experimental_positive": b["experimental_positive"], "status_as_is": b["status"],
               "status_option": o["status"], "ratio_as_is": rb, "ratio_option": ro}
        for k, (x, y) in cls.items():
            row[f"class_{k}_as_is"], row[f"class_{k}_option"] = x, y
        changes.append(row)
    elif rb is not None and ro is not None and abs(rb - ro) > 1e-6:
        ratio_moves.append({"gene": g, "ratio_as_is": rb, "ratio_option": ro})

cols = (["gene", "experimental_positive", "status_as_is", "status_option", "ratio_as_is", "ratio_option"]
        + [f"class_{rule}_{c}_{s}" for rule in "BA" for c in CUTS for s in ("as_is", "option")])
with open(HERE / "tables/sweep_class_changes.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(changes)

with open(HERE / "tables/sweep_scores.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["model", "rule_scope", "cutoff", "TP", "FN", "FP", "TN", "unresolved", "recall", "specificity",
                "precision", "mcc"])
    for name, s in (("as_is_353cefdc=7ccbb568", sb), ("option_B0/ALL", so)):
        for key, by_cut in s["scores"].items():
            for c, v in by_cut.items():
                w.writerow([name, key, c] + [v[k] for k in ("TP", "FN", "FP", "TN", "unresolved", "recall",
                                                             "specificity", "precision", "mcc")])

listed = {r["gene_id"] for r in csv.DictReader(open(LIST))}
out = {
    "wt_as_is": wt_b, "wt_option": wt_o, "wt_change": wt_o - wt_b,
    "genes": len(base), "class_changes": len(changes), "class_change_genes": [r["gene"] for r in changes],
    "ratio_moves_over_1e-6_without_class_change": ratio_moves,
    "listed_genes": len(listed), "listed_not_in_model": len(listed - set(base)),
    "option_certified": so["certified"], "option_not_certified_infeasible": so["not_certified_infeasible"],
    "option_fresh_resolves": so["fresh_resolves"], "option_failing_check": so["sweep_points_failing_check"],
    "option_magnitudes": so["magnitudes_in_sweep"], "option_lp_count": so["lp_count"],
    "baseline_certified_same": sb["certified"] == so["certified"],
}
# Cross-check: the step-4 targeted knockouts of B0/ALL must equal the sweep's rows for the same genes.
s4 = json.load(open(ART / "out/step4_options.json"))["variants"]["B0/ALL"]["ko"]
out["step4_crosscheck_max_abs_ratio_diff"] = max(
    abs((s4[g]["ratio"] or 0) - ((growth_b(opt[g]) or 0) / wt_o)) for g in KO_GENES if g in opt)
out["step4_crosscheck_missing"] = [g for g in KO_GENES if g not in opt]
if len(sys.argv) > 1:  # the as-is control in the step-6 harness must reproduce the baseline rows
    ctl = rows(sys.argv[1])
    sc = json.load(open(Path(sys.argv[1]).with_suffix(".summary.json")))
    out["control"] = {
        "wt": sc["wt"]["growth"], "wt_diff_vs_baseline": sc["wt"]["growth"] - wt_b,
        "genes": sorted(ctl),
        "max_abs_ratio_diff_vs_baseline": max(
            abs(((growth_b(ctl[g]) or 0) / sc["wt"]["growth"]) - ((growth_b(base[g]) or 0) / wt_b)) for g in ctl),
        "class_mismatch_vs_baseline": sorted(
            g for g in ctl if any(ctl[g][f"class_{r}_{c}"] != base[g][f"class_{r}_{c}"] for r in "AB" for c in CUTS)),
    }
(HERE / "tables/step6_compare.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "option_magnitudes"}, indent=1, default=str))
