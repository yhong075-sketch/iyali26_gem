"""Task C: reconstruct the setup behind the recalled V1 elastic-fit totals (6.94 / 1.14 / 2.72 / 0.68).

Diagnostic only, run after the stop rule fired (TASK.md). Grid: model {parent 0ac2addf with R132,
Task B candidate a56f3649 without R132} x glycine-cleavage row {module (pre-registered), -R361, +R361}
x the four setups, V1 medium, holistic overlay, f = 1. "+R361" maps R361 directly onto the MFA R51-net
interval [-0.208, 0.876], i.e. caps R361 at 0.876 per 100 glucose (Task B's "R361 looser").

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/reproduction_reconstruction.py
Writes reproduction_reconstruction.json and .csv here. Never overwrites.
"""
import csv
import json
import time
from datetime import datetime, timezone

import common
from common import HERE, MU, PINS, Counter_, apply_overlay, check_pins, load, sha
from elastic_fit import MAPPING, RECALLED_V1, SETUPS, ranges
from reproduction_diagnostics import PARENT, total

GCV_FORMS = {"module": {"R357": -1, "R360": -1, "R362": -1}, "minus_R361": {"R361": -1}, "plus_R361": {"R361": 1}}


def main():
    inputs = check_pins()
    assert sha(PARENT[0]) == PARENT[1], "parent model changed"
    started = time.monotonic()
    rng = ranges()
    results = []
    candidate = PINS["model"]
    for label, pin in (("parent_R132_present", PARENT), ("taskB_candidate_R132_removed", candidate)):
        PINS["model"] = pin
        try:
            sim, model = load("V1")
        finally:
            PINS["model"] = candidate
        for form, coeffs in GCV_FORMS.items():
            mapping = [(k, r, n, coeffs if k == "GCV" else c) for k, r, n, c in MAPPING]
            for setup, rev, with_malate in SETUPS:
                keys = [k for k, *_ in mapping if with_malate or k != "malate_import"]
                with model:
                    apply_overlay(model)
                    if rev:
                        common.r825_reversible(model)
                    with model:
                        t, v = total(model, rng, mapping, keys, MU)
                exp = RECALLED_V1[setup]
                results.append({"model": label, "gcv_form": form, "setup": setup, "total": t, "recalled_V1": exp,
                                "relative_difference": None if t is None else t / exp - 1, "violations": v})
                print(label, form, setup, None if t is None else round(t, 3), v, flush=True)
    by = {}
    for r in results:
        by.setdefault((r["model"], r["gcv_form"]), []).append(r["relative_difference"] is not None
                                                             and abs(r["relative_difference"]) <= 0.02)
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
              "parent_model": {"path": str(PARENT[0]), "sha256": PARENT[1]}, "medium": "V1", "overlay": True,
              "results": results,
              "all_four_within_2pct": {f"{m}|{g}": all(v) for (m, g), v in by.items()},
              "lps": Counter_.lps, "elapsed_seconds": time.monotonic() - started}
    (HERE / "reproduction_reconstruction.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "reproduction_reconstruction.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["model", "gcv_form", "setup", "recalled_V1_total", "total", "relative_difference", "violations_per100"])
        for r in results:
            w.writerow([r["model"], r["gcv_form"], r["setup"], r["recalled_V1"],
                        None if r["total"] is None else round(r["total"], 4),
                        None if r["relative_difference"] is None else round(r["relative_difference"], 4),
                        json.dumps(r["violations"])])
    print(record["all_four_within_2pct"])


if __name__ == "__main__":
    main()
