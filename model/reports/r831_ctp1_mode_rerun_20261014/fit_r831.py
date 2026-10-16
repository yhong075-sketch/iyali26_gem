"""C1 (DIAGNOSTIC): Task C's W29 13C elastic fit with R831 restricted in memory (asis, off, fwd0, rev0).

Imports Task C's elastic_fit.py (fit(), MAPPING, SETUPS, FRACTIONS, ranges()) and its common.py unchanged, by path
and SHA. Only Task C's INFO tuple is extended in memory, so the face spans also cover every mitochondrial malate /
citrate carrier. The variant closure is applied after the overlay (and R825), before the WT solve, exactly where
Task C's setup ends. `asis` must reproduce Task C's 24 totals (TASK.md).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/r831_ctp1_mode_rerun_20261014/fit_r831.py
Writes fit_r831_totals.csv, fit_r831_summary.csv and fit_r831_checks.json here; the full record goes to
$IYALI26_RESEARCH_ROOT/artifacts/tasks/r831_ctp1_mode_rerun_20261014/fit_r831.json. Never overwrites.
"""
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
TASK_C = REPO / "model/reports/w29_13c_consistency_20261012"
CODE_PINS = {TASK_C / "common.py": "487c95d9d32102d78d81815d4abc9e700cab0b09b02e3d3833ef0156649cc839",
             TASK_C / "elastic_fit.py": "c14d42ddf99f0bdec54df5c294c5d1c44fc6f6c6902f927afbb0a904e1bd9dea",
             TASK_C / "elastic_fit.json": "8d7d3bf690d09642159fdfcf9161e3d44caf629a168cfb7b7f54c4ec1229697e"}
ART = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve() / "artifacts/tasks" / HERE.name
sys.path.insert(0, str(TASK_C))
import elastic_fit as ec  # noqa: E402
import common as cc  # noqa: E402  (Task C's common, the one elastic_fit imports)

VARIANTS = {"asis": None, "off": (0.0, 0.0), "fwd0": (-1000.0, 0.0), "rev0": (0.0, 1000.0)}
EXTRA_INFO = ("R807", "R831", "R832", "R833", "R817", "R902", "R911", "R928", "R929", "R1386")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def restrict(model, bounds):
    """Intersect R831's bounds with the variant's (never opens anything)."""
    r = model.reactions.R831
    lb, ub = max(r.lower_bound, bounds[0]), min(r.upper_bound, bounds[1])
    assert lb <= ub
    r.bounds = (lb, ub)


def main():
    for path, digest in CODE_PINS.items():
        assert sha(path) == digest, f"Task C code changed: {path.name}"
    inputs = cc.check_pins()
    assert ec.INFO[-1] == "R2062" and not set(EXTRA_INFO) & set(ec.INFO)
    ec.INFO = ec.INFO + EXTRA_INFO
    reference = {(r["medium"], r["setup"], r["fraction_of_measured_growth"]): r.get("total_violation")
                 for r in json.loads((TASK_C / "elastic_fit.json").read_text())["runs"]}
    assert len(reference) == 24
    rng = ec.ranges()
    started = time.monotonic()
    record = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "command": "python -B model/reports/r831_ctp1_mode_rerun_20261014/fit_r831.py",
              "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
              "inputs": inputs, "code": {p.name: d for p, d in CODE_PINS.items()}, "fit_r831_sha": sha(__file__),
              "variants": VARIANTS, "info_added": EXTRA_INFO, "ranges_per_100_glucose": rng, "runs": []}
    for variant, bounds in VARIANTS.items():
        for medium in ("V2", "V1"):
            _, model = cc.load(medium)
            for name, rev, with_malate in ec.SETUPS:
                keys = [k for k, *_ in ec.MAPPING if with_malate or k != "malate_import"]
                with model:
                    cc.apply_overlay(model)
                    if rev:
                        cc.r825_reversible(model)
                    if bounds:
                        restrict(model, bounds)
                    r831 = list(model.reactions.R831.bounds)
                    wt, wt_status = cc.growth(model)
                    if wt is None:
                        raise SystemExit(f"STOP: WT infeasible {variant} {medium} {name}")
                    for f in ec.FRACTIONS:
                        with model:
                            out = ec.fit(model, rng, keys, f * cc.MU)
                        record["runs"].append({"variant": variant, "medium": medium, "setup": name,
                                               "fraction_of_measured_growth": f, "R831_bounds": r831, "wt_max": wt,
                                               "wt_status": wt_status, "ranges_used": keys, **out})
                        print(variant, medium, name, f, round(wt, 5), out.get("status"),
                              out.get("total_violation"), "LPs", cc.Counter_.lps, flush=True)
            if cc.Counter_.lps > 21000:
                raise SystemExit("STOP: C1 share of the LP budget exceeded")
        if variant == "asis":
            checks = []
            for run in (r for r in record["runs"] if r["variant"] == "asis"):
                key = (run["medium"], run["setup"], run["fraction_of_measured_growth"])
                exp, got = reference[key], run.get("total_violation")
                ok = (exp is None and got is None) or (exp is not None and got is not None and (
                    abs(got - exp) <= 1e-9 if exp == 0 else abs(got / exp - 1) <= 1e-6))
                checks.append({"medium": key[0], "setup": key[1], "fraction": key[2], "task_c": exp, "got": got,
                               "ok": ok})
            record["reproduction_checks"] = checks
            (HERE / "fit_r831_checks.json").open("x").write(json.dumps(checks, indent=1) + "\n")
            if not all(c["ok"] for c in checks):
                ART.mkdir(parents=True, exist_ok=True)
                (ART / "fit_r831_partial.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
                raise SystemExit("STOP: asis does not reproduce Task C's totals")
    record.update(lps=cc.Counter_.lps, elapsed_seconds=time.monotonic() - started)
    ART.mkdir(parents=True, exist_ok=True)
    (ART / "fit_r831.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "fit_r831_totals.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["variant", "medium", "setup", "fraction", "R831_bounds", "wt_max", "status", "total_violation",
                    "violated_ranges_point"])
        for run in record["runs"]:
            violated = ";".join(f"{k}{'+' if run['slack_point'][f'sp_{k}'] > 1e-6 else '-'}"
                                f"{max(run['slack_point'][f'sp_{k}'], run['slack_point'][f'sm_{k}']):.2f}"
                                for k in run.get("ranges_used", []) if run.get("status") == "optimal"
                                and max(run["slack_point"][f"sp_{k}"], run["slack_point"][f"sm_{k}"]) > 1e-6)
            w.writerow([run["variant"], run["medium"], run["setup"], run["fraction_of_measured_growth"],
                        run["R831_bounds"], round(run["wt_max"], 6), run.get("status"),
                        None if run.get("total_violation") is None else round(run["total_violation"], 6), violated])
    with (HERE / "fit_r831_summary.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["variant", "medium", "setup", "fraction", "item", "point", "face_min", "face_max"])
        for run in record["runs"]:
            if run.get("status") != "optimal":
                continue
            for item, span in run["face"].items():
                if item.startswith(("flux_", "info_")):
                    point = run["flux_point"].get(item[5:]) if item.startswith("flux_") else None
                    w.writerow([run["variant"], run["medium"], run["setup"], run["fraction_of_measured_growth"], item,
                                None if point is None else round(point, 4),
                                *(None if x is None else round(x, 4) for x in span)])
    print("LPs", cc.Counter_.lps, "elapsed", round(record["elapsed_seconds"]), "s")


if __name__ == "__main__":
    main()
