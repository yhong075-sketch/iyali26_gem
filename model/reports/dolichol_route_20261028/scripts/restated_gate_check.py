"""Check the restated knockout gate (TASK.md section 11) on the two-model sweep TSV (read-only).

Usage: python -I restated_gate_check.py ko_compare.tsv OUT.json
Gate: no knockout gains absolute growth above 1e-6 (the class half is in stop_condition.json).
"""
import csv
import json
import sys

tsv, out = sys.argv[1:3]
rows = list(csv.DictReader(open(tsv), delimiter="\t"))
gains = [(r["gene"], float(r["candidate_ko"]) - float(r["baseline_ko"])) for r in rows if r["baseline_ko"] and r["candidate_ko"]]
largest = max(gains, key=lambda t: t[1])
result = {"restated_gate": "no class change other than YALI1E22718g; no knockout gains absolute growth above 1e-6",
          "knockouts_compared": len(gains), "largest_absolute_gain": {"gene": largest[0], "gain": largest[1]},
          "gains_over_1e-6": [g for g in gains if g[1] > 1e-6]}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print(result)
