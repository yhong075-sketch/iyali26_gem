"""Task D step 4 (diagnostic): which capacity of the matrix-NADH outlet R533 (reverse) makes complex I essential?

Added after the tier-1 and step-5 results. Cap: R533 >= -cap x glucose uptake / 100 (cap in per 100 glucose), the
only matrix-NADH sink besides complex I (mechanism.json). Scanned with no 13C floor (R825 as supplied) and with
PDH + CS + G6PDH floors (R825 ASSUMPTION; the accepted step-5 set with the lowest KO/WT), configurations A and V2.
WT and KO get the same cap. ASSUMPTION (diagnostic): no measured capacity exists.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/shuttle_cap.py
Writes shuttle_cap.json and shuttle_cap.csv here. Never overwrites.
Optional: `shuttle_cap.py fine` runs a finer grid around the thresholds found by the audit (written to
shuttle_cap_fine.json / .csv); the default grid is unchanged. The fine grid gained the boundary points the
2026-10-14 audit identified (0.2357, 0.449; 86.66-87.2).
"""
import csv
import io
import json
import sys
import time
from datetime import datetime, timezone

from common import CONFIGS, CUTOFF, HERE, Counter_, check_pins, growth, ko, load, mfa_rows, r825_assumption, write_new
from ratios import FLOORS, add_floors

CAPS = (0.0, 0.5, 1.0, 2.0, 2.5, 3.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, None)
FINE = {"no_floor_R825_as_supplied": (0.1, 0.2, 0.2357, 0.236, 0.25, 0.3, 0.4, 0.449, 0.45, 0.5),
        "PDH+CS+G6PDH_R825_assumption": (85.0, 86.0, 86.5, 86.66, 86.67, 86.7, 86.8, 86.9, 87.0, 87.1, 87.2,
                                         88.0, 90.0, 95.0)}
SETTINGS = {"no_floor_R825_as_supplied": ([], False), "PDH+CS+G6PDH_R825_assumption": (["PDH", "CS", "G6PDH"], True)}


def main(tag=""):
    inputs = check_pins()
    started = time.monotonic()
    mfa = mfa_rows()
    fractions = {n: float(mfa[row]["glucose_LB"]) / 100.0 for n, (_, row) in FLOORS.items()}
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "fractions": fractions,
              "label": "ASSUMPTION (diagnostic): cap on R533 reverse", "rows": []}
    for config in CONFIGS:
        for setting, (floors, r825) in SETTINGS.items():
            wt, _ = load(config)
            kom, _ = load(config)
            ko(kom)
            for m in (wt, kom):
                if r825:
                    r825_assumption(m)
                add_floors(m, floors, fractions)
            w_open = growth(wt)[0]
            for cap in (FINE[setting] if tag else CAPS):
                res = {}
                for label, m in (("WT", wt), ("KO", kom)):
                    with m:
                        if cap is not None:
                            uptake = -m.reactions.R1070.flux_expression
                            m.add_cons_vars([m.problem.Constraint(m.reactions.R533.flux_expression + cap / 100.0 * uptake,
                                                                  lb=0, name="cap_R533_reverse")])
                        res[label] = growth(m)
                w, k = res["WT"][0], res["KO"][0]
                row = {"config": config, "setting": setting, "cap_per100": cap, "WT": w,
                       "WT_cost": None if w is None else 1 - w / w_open, "KO": k, "KO_status": res["KO"][1],
                       "KO_over_WT": None if not w or k is None else k / w}
                # benchmark.md: a KO with no solution is unresolved, not lethal (fixed after the 2026-10-14 audit).
                row["essential_15"] = (None if w is None else "unresolved" if k is None
                                       else row["KO_over_WT"] < CUTOFF)
                record["rows"].append(row)
                print(config, setting, cap, None if w is None else round(w, 5),
                      None if row["KO_over_WT"] is None else round(row["KO_over_WT"], 4), res["KO"][1], flush=True)
    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started)
    write_new(HERE / f"shuttle_cap{'_' + tag if tag else ''}.json", json.dumps(record, indent=1, default=str) + "\n")
    buffer = io.StringIO()
    w_ = csv.writer(buffer, lineterminator="\n")
    w_.writerow(["config", "setting", "cap_R533_reverse_per100", "WT", "WT_cost", "KO", "KO_status", "KO_over_WT",
                 "essential_15"])
    for r in record["rows"]:
        w_.writerow([r["config"], r["setting"], r["cap_per100"], r["WT"], r["WT_cost"], r["KO"], r["KO_status"],
                     r["KO_over_WT"], r["essential_15"]])
    write_new(HERE / f"shuttle_cap{'_' + tag if tag else ''}.csv", buffer.getvalue())


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "")
