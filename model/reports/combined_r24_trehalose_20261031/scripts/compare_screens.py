"""Compare the certified knockout screens of 164542fb and the combined candidate (read-only).

Usage: python compare_screens.py BASE.tsv BASE_SUMMARY.json NEW.tsv NEW_SUMMARY.json OUT.json
Reports WT growth of both, every gene whose class differs under rule A or rule B at any cutoff, the largest KO
growth change, and both summaries' score blocks side by side.
"""
import csv
import json
import math
import sys

base_tsv, base_json, new_tsv, new_json, out = sys.argv[1:6]
read = lambda p: {r["gene"]: r for r in csv.DictReader(open(p), delimiter="\t")}
base, new = read(base_tsv), read(new_tsv)
assert base.keys() == new.keys(), "gene sets differ"
class_cols = [c for c in next(iter(base.values())) if c.startswith("class_")]

changes = []
for g in sorted(base):
    diff = {c: [base[g][c], new[g][c]] for c in class_cols if base[g][c] != new[g][c]}
    if diff or base[g]["status"] != new[g]["status"]:
        changes.append({"gene": g, "status": [base[g]["status"], new[g]["status"]],
                        "ko": [base[g]["ko"], new[g]["ko"]], "classes": diff})


def ratio_delta(g):
    a, b = base[g]["ko"], new[g]["ko"]
    try:
        return abs(float(a) - float(b))
    except ValueError:
        return math.nan


deltas = {g: ratio_delta(g) for g in base}
finite = {g: d for g, d in deltas.items() if not math.isnan(d)}
top = sorted(finite, key=finite.get, reverse=True)[:10]
sb, sn = json.load(open(base_json)), json.load(open(new_json))
result = {
    "genes": len(base),
    "wt": {"base": sb["wt"], "new": sn["wt"]},
    "class_changes": changes,
    "largest_ko_growth_changes": [{"gene": g, "base": base[g]["ko"], "new": new[g]["ko"], "abs_delta": finite[g]}
                                  for g in top],
    "lp_count": {"base": sb["lp_count"], "new": sn["lp_count"]},
    "sweep_points_failing_check": {"base": sb["sweep_points_failing_check"], "new": sn["sweep_points_failing_check"]},
    "certified": {"base": sb["certified"], "new": sn["certified"]},
    "not_certified_infeasible": {"base": sb["not_certified_infeasible"], "new": sn["not_certified_infeasible"]},
    "scores": {"base": sb["scores"], "new": sn["scores"]},
}
json.dump(result, open(out, "w"), indent=1)
print(json.dumps({k: result[k] for k in ("genes", "wt", "class_changes", "lp_count", "sweep_points_failing_check")},
                 indent=1))
