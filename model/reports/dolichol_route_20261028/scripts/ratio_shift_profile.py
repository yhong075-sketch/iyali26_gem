"""Profile the knockout growth shifts behind the stop (read-only).

Usage: python -I ratio_shift_profile.py ko_compare_summary.json stop_condition.json OUT.json
For each gene outside the route whose KO/WT ratio moved by more than 1e-6: the relative change of its absolute KO
growth and of its ratio, compared with the WT relative change; counts and extremes, plus the largest ratio shift.
"""
import json
import sys

summary, stop, out = sys.argv[1:4]
s, c = json.load(open(summary)), json.load(open(stop))
wt_rel = c["wt_relative_change"]
rows = []
for e in c["ratio_changes_outside_route"]:
    gb, gc = e["baseline_ratio"] * c["wt"]["baseline"], e["candidate_ratio"] * c["wt"]["candidate"]
    rows.append({"gene": e["gene"], "experimental_positive": e["experimental_positive"],
                 "baseline_ratio": e["baseline_ratio"], "ratio_change": e["candidate_ratio"] - e["baseline_ratio"],
                 "ko_growth_relative_change": (gc - gb) / gb})
rows.sort(key=lambda r: -abs(r["ratio_change"]))
result = {"wt_relative_change": wt_rel, "genes": len(rows),
          "max_abs_ratio_change": max(abs(r["ratio_change"]) for r in rows),
          "ko_growth_relative_change_range": [min(r["ko_growth_relative_change"] for r in rows), max(r["ko_growth_relative_change"] for r in rows)],
          "ratio_range_of_affected_genes": [min(r["baseline_ratio"] for r in rows), max(r["baseline_ratio"] for r in rows)],
          "nearest_to_a_cutoff": sorted(({"gene": r["gene"], "ratio": r["baseline_ratio"], "cutoff": k, "distance": abs(r["baseline_ratio"] - k)}
                                         for r in rows for k in (0.01, 0.05, 0.10, 0.15)), key=lambda x: x["distance"])[:5],
          "experimental_positives_among_them": sum(r["experimental_positive"] == "True" for r in rows),
          "rows": rows}
json.dump(result, open(out, "w"), indent=1)
open(out, "a").write("\n")
print({k: v for k, v in result.items() if k != "rows"})
print("top 5:", [(r["gene"], round(r["baseline_ratio"], 6), f"{r['ratio_change']:+.2e}", f"{r['ko_growth_relative_change']:+.2e}") for r in rows[:5]])
