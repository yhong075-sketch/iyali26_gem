"""C4 (DIAGNOSTIC): gene benchmark with R831 off or reverse-closed, benchmark context (as loaded, no overlay).

Reuses today's screen runner (complexI_shuttle_assumption_screen_20261014/run_screen.py, SHA-pinned): Task B's
screen()/metrics() unchanged, Task D's close(); only the variant table is replaced in memory. fwd0 is that screen's
option (b) and is not re-run. The base run must reproduce 82/33/229/697/16 at 15 %.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/r831_ctp1_mode_rerun_20261014/screen_r831.py
Writes screen_r831.json and screen_r831_changes.csv here; per-gene tables go to
$IYALI26_RESEARCH_ROOT/artifacts/tasks/r831_ctp1_mode_rerun_20261014/screen/. Never overwrites.
"""
import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RUNNER = REPO / "model/reports/complexI_shuttle_assumption_screen_20261014/run_screen.py"
RUNNER_SHA = "c53c8c21cebef614de903c4d4787746c3c26c9ec58ea6c8e110c1c7159f6a173"
OUT = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve() / "artifacts/tasks" / HERE.name / "screen"
assert hashlib.sha256(RUNNER.read_bytes()).hexdigest() == RUNNER_SHA
sys.path.insert(0, str(RUNNER.parent))
import run_screen as rs  # noqa: E402

VARIANTS = {"base": {}, "off": {"R831": "whole"}, "rev0": {"R831": "rev"}}
EXPECTED_BASE = dict(TP=82, FP=33, FN=229, TN=697, unresolved=16)


def main():
    rs.td.check_pins()
    rs.VARIANTS = VARIANTS
    OUT.mkdir(parents=True, exist_ok=False)
    positives, _ = rs.tb.load_experimental_reference(rs.tb.INPUTS["experimental"][0])
    positives = set(positives)
    mtdna = json.loads(rs.tb.INPUTS["mtdna_scope"][0].read_text())["genes"]
    started = time.monotonic()
    record = {"runner_sha": RUNNER_SHA, "screen_r831_sha": rs.sha(Path(__file__)), "variants": VARIANTS, "runs": {}}
    screens = {}
    for variant in VARIANTS:
        applied = {}
        rs.tb.context = rs.wrapped_context("bench", variant, applied)
        result, rows = rs.tb.screen("candidate", positives, mtdna)
        screens[variant] = rows
        run = {k: v for k, v in result.items() if k != "route"}
        run["applied"] = applied
        run["route"] = {k: {"growth": r["growth"], "status": r["status"], "fva": r["fva"]}
                        for k, r in result["route"].items()}
        record["runs"][variant] = run
        with (OUT / f"gene_screen_{variant}.tsv").open("x") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        primary = next(s for s in run["scores"] if s["cutoff"] == 0.15 and s["view"] == "benchmark rule")
        print(variant, run["WT"], {k: v["ratio"] for k, v in run["reaction_ko"].items()},
              {k: primary[k] for k in ("TP", "FP", "FN", "TN", "unresolved", "mcc")}, flush=True)
        if variant == "base":
            got = {k: primary[k] for k in EXPECTED_BASE}
            record["base_check"] = {"expected": EXPECTED_BASE, "got": got, "ok": got == EXPECTED_BASE}
            if not record["base_check"]["ok"]:
                (HERE / "screen_r831.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
                raise SystemExit("STOP: base does not reproduce the benchmark")
    base = {r["gene_id"]: r for r in screens["base"]}
    changes = []
    for variant in ("off", "rev0"):
        for r in screens[variant]:
            b = base[r["gene_id"]]
            if (b["status"] == r["status"] and (b["ratio"] is None) == (r["ratio"] is None)
                    and (r["ratio"] is None or abs(r["ratio"] - b["ratio"]) <= 1e-6)):
                continue
            row = {"variant": variant, "gene_id": r["gene_id"], "positive": r["positive"], "mtdna": r["mtdna"],
                   "reactions": r["reactions"], "base_status": b["status"], "variant_status": r["status"],
                   "base_ratio": b["ratio"], "variant_ratio": r["ratio"]}
            for c in rs.tb.CUTOFFS:
                def call(x, c=c):
                    return "unresolved" if x is None else ("E" if x < c else "NE")
                row[f"call_change_{c:g}"] = (f"{call(b['ratio'])}->{call(r['ratio'])}"
                                             if call(b["ratio"]) != call(r["ratio"]) else "")
            changes.append(row)
    record["gene_differences"] = changes
    record.update(lps=rs.tb.Counter_.lps, elapsed_seconds=time.monotonic() - started)
    (HERE / "screen_r831.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "screen_r831_changes.csv").open("x") as f:
        fields = list(changes[0]) if changes else ["variant", "gene_id"]
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(changes)
    print("gene differences", len(changes), "LPs", rs.tb.Counter_.lps)


if __name__ == "__main__":
    main()
