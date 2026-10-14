"""Task C2: holistic test of the matrix acetyl-unit export module (prompt rule 1), diagnostic only.

Phenotype: complex I (R2062) KO/WT under PDH >= 0.866629 x glucose uptake. Parts (k = 5, all 2^5
closure combinations, one LP each for WT and KO):
  P1 R825_export        acetylcarnitine export: open = R825 reversible; closed = as supplied (import only)
  P2 R1172              HMG-CoA transport matrix <-> cytosol (no GPR)
  P3 R1996              acetoacetate transport matrix <-> cytosol (no GPR)
  P4 R1894              ATP-citrate lyase (citrate route to cytosolic acetyl-CoA)
  P5 isoprenoid_exch    R945, R1106, R1108, R1124, R1235, R1369, R1388, R1390 (mevalonate and 7 sterol
                        exchanges; grouped: same mechanism, secretion of an isoprenoid-pathway intermediate)
Conditions: V2, holistic overlay, caps {none, module}, PDH floor {on, off}.
"All open" (reference) = R825 reversible with every other part as supplied.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/c2_acetyl_module.py
Writes c2_acetyl_module.json and .csv here. Never overwrites.
"""
import csv
import itertools
import json
import time
from datetime import datetime, timezone

from common import HERE, Counter_, add_caps, apply_overlay, check_pins, growth, load, load_caps, pdh_floor, r825_reversible

PARTS = ("R825_export", "R1172", "R1996", "R1894", "isoprenoid_exch")
ISOPRENOID_EXCHANGES = ("R945", "R1106", "R1108", "R1124", "R1235", "R1369", "R1388", "R1390")
CUTOFF = 0.15


def close(model, part):
    if part == "R825_export":
        return  # closed = as supplied (lb 0)
    if part == "isoprenoid_exch":
        for rid in ISOPRENOID_EXCHANGES:
            r = model.reactions.get_by_id(rid)
            assert tuple(r.bounds) == (0.0, 1000.0), rid
            r.bounds = (0.0, 0.0)
        return
    model.reactions.get_by_id(part).bounds = (0.0, 0.0)


def main():
    inputs = check_pins()
    caps = load_caps()
    started = time.monotonic()
    sim, model = load("V2")
    rows = []
    for capped in (False, True):
        for floor in (True, False):
            for closed in itertools.chain.from_iterable(itertools.combinations(PARTS, n) for n in range(len(PARTS) + 1)):
                out = {}
                for label, ko in (("WT", None), ("KO", "R2062")):
                    with model:
                        apply_overlay(model)
                        if capped:
                            add_caps(model, caps, "module")
                        if floor:
                            pdh_floor(model)
                        if "R825_export" not in closed:
                            r825_reversible(model)
                        for p in closed:
                            close(model, p)
                        if ko:
                            model.reactions.get_by_id(ko).bounds = (0.0, 0.0)
                        out[label] = growth(model)
                wt, ko = out["WT"][0], out["KO"][0]
                rows.append({"caps": "module" if capped else None, "pdh_floor": floor, "closed": list(closed),
                             "WT": wt, "WT_status": out["WT"][1], "KO": ko, "KO_status": out["KO"][1],
                             "KO_over_WT": None if not wt or ko is None else ko / wt})
    # Minimal cut sets and leave-one-out, per condition.
    analysis = {}
    for capped in (None, "module"):
        for floor in (True, False):
            sub = [r for r in rows if r["caps"] == capped and r["pdh_floor"] == floor]
            ref = next(r for r in sub if not r["closed"])
            hit = lambda r: r["WT"] is not None and (r["KO"] is None or r["KO_over_WT"] < CUTOFF)
            hits = [set(r["closed"]) for r in sub if hit(r)]
            minimal = [s for s in hits if not any(o < s for o in hits)]
            full = next(r for r in sub if len(r["closed"]) == len(PARTS))
            loo = {p: next(r for r in sub if set(r["closed"]) == set(PARTS) - {p}) for p in PARTS}
            key = f"caps={capped}|floor={floor}"
            analysis[key] = {
                "reference_all_open": ref,
                "minimal_cut_sets": [{"closed": sorted(s),
                                      "WT": next(r for r in sub if set(r["closed"]) == s)["WT"],
                                      "WT_cost": 1 - next(r for r in sub if set(r["closed"]) == s)["WT"] / ref["WT"],
                                      "KO_over_WT": next(r for r in sub if set(r["closed"]) == s)["KO_over_WT"],
                                      "KO_status": next(r for r in sub if set(r["closed"]) == s)["KO_status"]}
                                     for s in minimal],
                "full_set": full,
                "leave_one_out": {p: {"KO_over_WT": r["KO_over_WT"], "KO_status": r["KO_status"], "WT": r["WT"]}
                                  for p, r in loo.items()}}
            print(key, "ref", round(ref["WT"], 4), ref["KO_over_WT"], "| minimal", analysis[key]["minimal_cut_sets"],
                  flush=True)
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs, "parts": PARTS,
              "isoprenoid_exchanges": ISOPRENOID_EXCHANGES, "cutoff": CUTOFF, "rows": rows, "analysis": analysis,
              "note": "KO with no solution counts as a hit here only to locate the constraint conflict; under the "
                      "benchmark rule it is unresolved, not lethal.",
              "lps": Counter_.lps, "elapsed_seconds": time.monotonic() - started}
    (HERE / "c2_acetyl_module.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "c2_acetyl_module.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["caps", "pdh_floor", "closed", "WT", "WT_status", "KO", "KO_status", "KO_over_WT"])
        for r in rows:
            w.writerow([r["caps"], r["pdh_floor"], "+".join(r["closed"]) or "none", r["WT"], r["WT_status"], r["KO"],
                        r["KO_status"], r["KO_over_WT"]])


if __name__ == "__main__":
    main()
