"""Apply the TASK.md section 7 stop condition to the two-model sweep (read-only).

Usage: python -I stop_condition_check.py ko_compare_summary.json ko_compare.tsv OUT.json
Absolute knockout growth moves with the WT shift, so growth change is judged on KO/WT ratios (each model's own
WT). Stop if any class changes other than YALI1E22718g non-essential -> essential, or any |ratio change| > 1e-6
for a gene outside the route (route gene: YALI1E22718g, the only gene whose reactions changed state).
"""
import csv
import json
import sys

summary_path, tsv_path, out = sys.argv[1:4]
s = json.load(open(summary_path))
wt = {k: v["growth"] for k, v in s["wt"].items()}
rows = list(csv.DictReader(open(tsv_path), delimiter="\t"))
ROUTE_GENES = {"YALI1E22718g"}
ratio_changes, unresolved_mismatch = [], []
for r in rows:
    b, c = r["baseline_ko"], r["candidate_ko"]
    if (b == "") != (c == ""):
        unresolved_mismatch.append(r["gene"])
        continue
    if b == "":
        continue
    rb, rc = float(b) / wt["baseline"], float(c) / wt["candidate"]
    if abs(rc - rb) > 1e-6:
        ratio_changes.append({"gene": r["gene"], "experimental_positive": r["experimental_positive"], "baseline_ratio": rb,
                              "candidate_ratio": rc, "route_gene": r["gene"] in ROUTE_GENES})
classes = s["candidate_vs_baseline"]["class_changed"]
allowed = lambda e: e["gene"] in ROUTE_GENES and e["baseline"] == "non-essential" and e["candidate"] == "essential"
bad_class = {c: [e for e in v if not allowed(e)] for c, v in classes.items()}
result = {"wt": wt, "wt_relative_change": (wt["candidate"] - wt["baseline"]) / wt["baseline"],
          "class_changed": classes, "class_changes_violating": bad_class,
          "ratio_changes_over_1e-6": ratio_changes,
          "ratio_changes_outside_route": [e for e in ratio_changes if not e["route_gene"]],
          "unresolved_mismatch": unresolved_mismatch}
result["stop"] = bool(any(bad_class.values()) or result["ratio_changes_outside_route"] or unresolved_mismatch)
json.dump(result, open(out, "w"), indent=1)
open(out, "a").write("\n")
print("WT rel change", result["wt_relative_change"], "| class changes", {c: [(e["gene"], e["baseline"], e["candidate"]) for e in v] for c, v in classes.items()})
print("ratio changes >1e-6:", [(e["gene"], round(e["baseline_ratio"], 6), round(e["candidate_ratio"], 6)) for e in ratio_changes])
print("unresolved mismatch:", unresolved_mismatch, "| STOP:", result["stop"])
