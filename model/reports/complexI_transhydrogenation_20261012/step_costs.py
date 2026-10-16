"""Task D: single-closure WT cost and KO/WT of candidate WT-free steps for one discovery round (configuration A).

The round's state is discovery round r (every part with round_added < r closed, WT and KO). Each candidate
closure "reaction:direction" is then added alone. WT cost is given against the round's WT and against the
all-open WT. Used to pick, per route, the WT-free step with the smallest WT cost (TASK.md).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/step_costs.py <r> R511:fwd,R328:fwd,...
Writes step_costs_round<r>.json here. Never overwrites.
"""
import json
import sys
import time
from datetime import datetime, timezone

from common import HERE, Counter_, check_pins, close, growth, ko, load, sha, write_new
from discovery import merged_closures


def main(round_, steps):
    inputs = check_pins()
    started = time.monotonic()
    definition_path = HERE / "parts_definition.json"
    definition = json.loads(definition_path.read_text())
    closures = merged_closures([p for p in definition["parts"] if p["round_added"] < round_])
    wt, _ = load("A")
    w_open = growth(wt)[0]
    kom, _ = load("A")
    ko(kom)
    for m in (wt, kom):
        close(m, closures)
    w_round, k_round = growth(wt)[0], growth(kom)[0]
    rows = []
    for step in steps:
        rid, how = step.split(":")
        res = {}
        for label, m in (("WT", wt), ("KO", kom)):
            with m:
                close(m, {rid: how})
                res[label] = growth(m)
        w, k = res["WT"][0], res["KO"][0]
        rows.append({"step": step, "name": wt.reactions.get_by_id(rid).name, "WT": w,
                     "WT_cost_vs_round": None if w is None else 1 - w / w_round,
                     "WT_cost_vs_open": None if w is None else 1 - w / w_open,
                     "KO": k, "KO_status": res["KO"][1], "KO_over_WT": None if not w or k is None else k / w})
        print(f"{step:14s} WT cost {rows[-1]['WT_cost_vs_round']} KO/WT {rows[-1]['KO_over_WT']}", flush=True)
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "round": round_,
              "parts_definition_sha256": sha(definition_path), "closures": closures, "WT_open": w_open,
              "WT_round": w_round, "KO_round": k_round, "steps": rows, "lps": Counter_.lps,
              "elapsed_seconds": time.monotonic() - started}
    write_new(HERE / f"step_costs_round{round_}.json", json.dumps(record, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2].split(","))
