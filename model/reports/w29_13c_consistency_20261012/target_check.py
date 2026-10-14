"""Task C: check the pre-registered elastic-fit totals against the user's amended targets (Amendment 1).

Reads elastic_fit.json (f = 1 runs) and writes target_check.json. Run from the repository root:
    .venv/bin/python -B model/reports/w29_13c_consistency_20261012/target_check.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGETS = {("V2", "overlay"): 6.50, ("V1", "overlay"): 5.47, ("V2", "overlay_no_malate"): 2.54,
           ("V1", "overlay_no_malate"): 2.06, ("V2", "R825_rev"): 1.44, ("V1", "R825_rev"): 1.21,
           ("V2", "R825_rev_no_malate"): 0.90, ("V1", "R825_rev_no_malate"): 0.68}
fit = json.loads((HERE / "elastic_fit.json").read_text())
checks = []
for (medium, setup), target in TARGETS.items():
    run = next(r for r in fit["runs"] if r["medium"] == medium and r["setup"] == setup
               and r["fraction_of_measured_growth"] == 1.0)
    rel = run["total_violation"] / target - 1
    checks.append({"medium": medium, "setup": setup, "target": target, "observed": run["total_violation"],
                   "relative_difference": rel, "within_2pct": abs(rel) <= 0.02})
out = {"source": "elastic_fit.json", "targets": "user Amendment 1, 2026-10-13", "checks": checks,
       "stop_fired": not all(c["within_2pct"] for c in checks)}
(HERE / "target_check.json").open("x").write(json.dumps(out, indent=1) + "\n")
for c in checks:
    print(c["medium"], c["setup"], c["target"], round(c["observed"], 4), f'{c["relative_difference"]:+.4f}')
print("stop fired:", out["stop_fired"])
