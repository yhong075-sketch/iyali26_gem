"""Task D: whole-reaction closure cost of the "capacity" steps (TASK.md part rules; run after the audit found it missing).

A capacity step is a reaction the KO uses in a direction the WT also requires (no escape direction). TASK.md says
it is not closed in the main test and its whole-reaction closure cost is reported. Here each one is closed whole
(bounds 0, 0), alone, from the open state, in WT and the complex I KO, configurations A and V2.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/capacity_costs.py
Writes capacity_costs.json here. Never overwrites.
"""
import json
import time
from datetime import datetime, timezone

from common import CONFIGS, HERE, Counter_, check_pins, close, growth, ko, load, write_new

# Capacity steps named in parts_pool.json (round 0) and in the discovery rounds (KO above WT, WT-required direction).
STEPS = ("R134", "R330", "R691", "R361", "R534", "R641", "R660", "R659", "R740", "R929", "R490", "R315", "R1385",
         "R570", "R325")


def main():
    inputs = check_pins()
    started = time.monotonic()
    rows = []
    for config in CONFIGS:
        wt, _ = load(config)
        kom, _ = load(config)
        ko(kom)
        w_open = growth(wt)[0]
        for rid in STEPS:
            res = {}
            for label, m in (("WT", wt), ("KO", kom)):
                with m:
                    close(m, {rid: "whole"})
                    res[label] = growth(m)
            w, k = res["WT"][0], res["KO"][0]
            rows.append({"config": config, "reaction": rid, "name": wt.reactions.get_by_id(rid).name, "WT": w,
                         "WT_status": res["WT"][1], "WT_cost": None if w is None else 1 - w / w_open, "KO": k,
                         "KO_status": res["KO"][1], "KO_over_WT": None if not w or k is None else k / w})
            print(config, rid, None if w is None else round(1 - w / w_open, 4),
                  None if rows[-1]["KO_over_WT"] is None else round(rows[-1]["KO_over_WT"], 4), res["KO"][1], flush=True)
    write_new(HERE / "capacity_costs.json", json.dumps(
        {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "steps": STEPS, "rows": rows,
         "lps": Counter_.lps, "elapsed_seconds": time.monotonic() - started}, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
