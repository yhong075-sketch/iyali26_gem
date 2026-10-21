"""Diagnostic added after the screen (not pre-registered): complex I KO/WT vs the lysine uptake cap R1220.

Under (a)/(b) the complex I KO grows on medium lysine alone (R1220 at its cap), so the KO/WT ratio scales with
that cap. Bench context only (as loaded, no overlay), Task B context() + Task D close(), same solver settings.

Run from the repository root:
    .venv/bin/python -B model/reports/complexI_shuttle_assumption_screen_20261014/lysine_sweep.py
Writes lysine_sweep.json and lysine_sweep.csv here. Never overwrites.
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_screen as rs  # noqa: E402

FACTORS = (0.0, 0.5, 1.0, 1.5, 2.0, 2.1, 2.2, 2.5, 3.0, 5.0, 10.0)


def main():
    rs.td.check_pins()
    _, model = rs.ORIGINAL_CONTEXT("candidate", "glucose")
    cap = -model.reactions.R1220.lower_bound
    assert abs(cap - 0.00595) < 1e-9, cap
    rows = []
    for variant, closures in rs.VARIANTS.items():
        for factor in FACTORS:
            with model:
                rs.td.close(model, closures)
                model.reactions.R1220.lower_bound = -cap * factor
                wt, wt_status = rs.tb.growth(model)
                model.reactions.R2062.bounds = (0.0, 0.0)
                ko, ko_status = rs.tb.growth(model)
            rows.append({"variant": variant, "lysine_factor": factor, "lysine_cap": cap * factor, "WT": wt,
                         "WT_status": wt_status, "KO": ko, "KO_status": ko_status,
                         "KO_over_WT": None if ko is None or not wt else ko / wt})
            print(rows[-1], flush=True)
    here = Path(__file__).resolve().parent
    with (here / "lysine_sweep.csv").open("x") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (here / "lysine_sweep.json").open("x").write(json.dumps(
        {"base_cap": cap, "lps": rs.tb.Counter_.lps, "code": rs.sha(Path(__file__)), "rows": rows}, indent=1) + "\n")


if __name__ == "__main__":
    main()
