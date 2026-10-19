"""Task D step 5: W29 13C flux ratios as an ASSUMPTION overlay on complex I KO/WT (diagnostic only).

Floors = the 95 % lower bounds of the pinned MFA file (glucose_LB / 100) x actual glucose uptake:
PDH R693 (MFA R19), CS R236 (R20), SDH R740 (R26 net), G6PDH R325 (R11). All 16 subsets, configurations A and
V2, R825 as supplied and with the R825 ASSUMPTION file, WT and KO with the same floors. WT cost is against the WT
with no floor in the same R825 state (and against the as-supplied WT). A floor set is accepted only if its WT
cost is <= 2 %. For the KO, FVA at the optimum of the matrix-NADH outlet (R533) and the malate exporter (R831).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/complexI_transhydrogenation_20261012/ratios.py
Writes ratios.json and ratios.csv here. Never overwrites.
"""
import csv
import io
import itertools
import json
import time
from datetime import datetime, timezone

from common import (CONFIGS, CUTOFF, HERE, Counter_, check_pins, fva, growth, ko, load, mfa_rows, per100,
                    r825_assumption, write_new)

FLOORS = {"PDH": ("R693", "R19"), "CS": ("R236", "R20"), "SDH": ("R740", "R26 net"), "G6PDH": ("R325", "R11")}
MAX_WT_COST = 0.02


def add_floors(model, names, fractions):
    uptake = -model.reactions.R1070.flux_expression
    model.add_cons_vars([model.problem.Constraint(model.reactions.get_by_id(FLOORS[n][0]).flux_expression
                                                  - fractions[n] * uptake, lb=0, name=f"floor_{n}") for n in names])


def main():
    inputs = check_pins()
    started = time.monotonic()
    mfa = mfa_rows()
    fractions = {n: float(mfa[row]["glucose_LB"]) / 100.0 for n, (_, row) in FLOORS.items()}
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "fractions": fractions,
              "label": "ASSUMPTION (diagnostic): W29 13C ratio floors; rows with r825='assumption' also use "
                       "R825_REVERSIBLE_ASSUMPTION", "rows": []}
    names = list(FLOORS)
    for config in CONFIGS:
        reference = None
        for r825 in ("as_supplied", "assumption"):
            wt, _ = load(config)
            kom, _ = load(config)
            ko(kom)
            if r825 == "assumption":
                r825_assumption(wt)
                r825_assumption(kom)
            base_wt = None
            for n in range(len(names) + 1):
                for combo in itertools.combinations(names, n):
                    res = {}
                    for label, m in (("WT", wt), ("KO", kom)):
                        n_constraints = len(m.constraints)
                        with m:
                            add_floors(m, combo, fractions)
                            res[label] = growth(m)
                            if label == "KO" and res[label][0]:
                                res["KO_fva"] = {rid: [round(per100(x), 3) for x in v]
                                                 for rid, v in fva(m, ["R533", "R831", "R693", "R570"], 1.0).items()}
                        assert len(m.constraints) == n_constraints, "floor constraints were not removed"
                    w, k = res["WT"][0], res["KO"][0]
                    if not combo:
                        base_wt = w
                        if r825 == "as_supplied":
                            reference = w
                    row = {"config": config, "r825": r825, "floors": list(combo), "WT": w, "WT_status": res["WT"][1],
                           "KO": k, "KO_status": res["KO"][1],
                           "KO_over_WT": None if not w or k is None else k / w,
                           "WT_cost": None if w is None else 1 - w / base_wt,
                           "WT_cost_vs_as_supplied": None if w is None else 1 - w / reference,
                           "KO_fva_per100": res.get("KO_fva")}
                    row["accepted"] = row["WT_cost"] is not None and row["WT_cost"] <= MAX_WT_COST
                    row["complex_I_essential_15"] = w is not None and (k is None or row["KO_over_WT"] < CUTOFF)
                    record["rows"].append(row)
                    print(config, r825, "+".join(combo) or "none", None if w is None else round(w, 5),
                          None if row["KO_over_WT"] is None else round(row["KO_over_WT"], 4), res["KO"][1],
                          "accepted" if row["accepted"] else "rejected", flush=True)
    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started)
    write_new(HERE / "ratios.json", json.dumps(record, indent=1, default=str) + "\n")
    buffer = io.StringIO()
    w = csv.writer(buffer, lineterminator="\n")
    w.writerow(["config", "r825", "floors", "WT", "WT_cost", "WT_cost_vs_as_supplied", "accepted", "KO", "KO_status",
                "KO_over_WT", "R533_KO_fva_per100", "R831_KO_fva_per100"])
    for r in record["rows"]:
        f = r["KO_fva_per100"] or {}
        w.writerow([r["config"], r["r825"], "+".join(r["floors"]) or "none", r["WT"], r["WT_cost"],
                    r["WT_cost_vs_as_supplied"], r["accepted"], r["KO"], r["KO_status"], r["KO_over_WT"],
                    f.get("R533"), f.get("R831")])
    write_new(HERE / "ratios.csv", buffer.getvalue())


if __name__ == "__main__":
    main()
