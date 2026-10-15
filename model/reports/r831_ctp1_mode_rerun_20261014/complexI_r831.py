"""C2 + C3 (DIAGNOSTIC): complex I KO/WT with R831 restricted in memory (asis, off, fwd0, rev0).

Imports Task D's common.py, mechanism.py (load_variant) and ratios.py (FLOORS, add_floors) unchanged, by path and
SHA. A KO with no solution is unresolved, never lethal (benchmark.md).
C2: robustness grid (A/V2 x overlay on/off x caps on/off) per variant; FVA at the WT and KO optima (overlay + caps
    on); lysine check for KOs below 15 % in A; floor row (R825 ASSUMPTION + PDH + CS + G6PDH).
C3: WHAT-IF with R831 off (overlay + caps on): W1 R902 reverse, W2 R1386 reverse, W3 both (lb -1000).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/r831_ctp1_mode_rerun_20261014/complexI_r831.py
Writes complexI_r831.json and complexI_r831.csv here. Never overwrites.
"""
import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
TASK_D = REPO / "model/reports/complexI_transhydrogenation_20261012"
CODE_PINS = {TASK_D / "common.py": "bfe0b2ced3d85cab77fae376cc89dcb220651f747ce25aaf5666707ef70b8ab3",
             TASK_D / "mechanism.py": "0a939031e73eee724dc20123ef1ba6d3df1bc0645f7bc0e7039f2b1aea92d527",
             TASK_D / "mechanism.json": "a26b8f7430abfe118da25873f9d9df81141dcecc7f619d4bcb1aaffe205d6952",
             TASK_D / "ratios.py": "d00481bf77f019f3083fc5fbdfa153e53f4f1a29d2519a35e868419647307471",
             TASK_D / "ratios.json": "fd7e409dc229432b31b988429fb46968ea5d93af353cfa82ecabc91afa381ce4"}
sys.path.insert(0, str(TASK_D))
import common as dc  # noqa: E402  (Task D's common)
from mechanism import load_variant  # noqa: E402
from ratios import add_floors  # noqa: E402

VARIANTS = {"asis": None, "off": (0.0, 0.0), "fwd0": (-1000.0, 0.0), "rev0": (0.0, 1000.0)}
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
FVA_RXNS = ("R831", "R832", "R833", "R817", "R807", "R902", "R1386", "R928", "R929", "R911", "R533", "R398", "R1220",
            "R570", "R2062")
WHATIF = {"W1_R902_reverse": ("R902",), "W2_R1386_reverse": ("R1386",), "W3_both": ("R902", "R1386")}
FLOOR_SET = ["PDH", "CS", "G6PDH"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def restrict(model, bounds):
    if bounds:
        r = model.reactions.R831
        lb, ub = max(r.lower_bound, bounds[0]), min(r.upper_bound, bounds[1])
        assert lb <= ub
        r.bounds = (lb, ub)


def calls(ratio, ko_status):
    if ratio is None:
        return {f"{c:g}": "unresolved" for c in CUTOFFS}
    return {f"{c:g}": ratio < c for c in CUTOFFS}


def pair(wt, kom, setup):
    """WT and KO growth under `setup(model)`, inside model contexts."""
    with wt:
        setup(wt)
        w, ws = dc.growth(wt)
    with kom:
        setup(kom)
        k, ks = dc.growth(kom)
    ratio = None if not w or k is None else k / w
    return {"WT": w, "WT_status": ws, "KO": k, "KO_status": ks, "KO_over_WT": ratio, "calls": calls(ratio, ks)}


def close_enough(a, b):
    if a is None or b is None:
        return a is None and b is None
    return abs(a - b) <= 1e-9 if b == 0 else abs(a / b - 1) <= 1e-6


def main():
    for path, digest in CODE_PINS.items():
        assert sha(path) == digest, f"Task D code changed: {path.name}"
    inputs = dc.check_pins()
    mech = json.loads((TASK_D / "mechanism.json").read_text())
    ratios_ref = json.loads((TASK_D / "ratios.json").read_text())
    fractions = ratios_ref["fractions"]
    started = time.monotonic()
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
              "code": {p.name: d for p, d in CODE_PINS.items()}, "complexI_r831_sha": sha(__file__),
              "variants": VARIANTS, "grid": [], "fva": {}, "lysine": [], "floors": [], "whatif": [], "checks": []}
    rows = []

    # C2 robustness grid.
    for config in dc.CONFIGS:
        for overlay in (True, False):
            for caps in (True, False):
                wt = load_variant(config, overlay, caps)
                kom = load_variant(config, overlay, caps)
                dc.ko(kom)
                w_asis = dc.growth(wt)[0]
                if w_asis is None:
                    raise SystemExit(f"STOP: WT infeasible {config} overlay={overlay} caps={caps}")
                ref = next(r for r in mech["robustness"] if r["config"] == config and r["overlay"] == overlay
                           and r["caps"] == caps)
                for variant, bounds in VARIANTS.items():
                    res = pair(wt, kom, lambda m, b=bounds: restrict(m, b))
                    res.update(config=config, overlay=overlay, caps=caps, variant=variant,
                               WT_cost=None if res["WT"] is None else 1 - res["WT"] / w_asis)
                    record["grid"].append(res)
                    rows.append(res)
                    label = {"asis": "none", "fwd0": "G3"}.get(variant)
                    if label:
                        ok = close_enough(res["KO_over_WT"], ref[label]["KO_over_WT"])
                        record["checks"].append({"what": f"mechanism.json {config} overlay={overlay} caps={caps} "
                                                         f"{label}", "expected": ref[label]["KO_over_WT"],
                                                 "got": res["KO_over_WT"], "ok": ok})
                    print(config, overlay, caps, variant, res["WT"], res["KO_over_WT"], res["KO_status"], flush=True)
                    # Lysine check for KOs below 15 % in A.
                    if config == "A" and res["KO_over_WT"] is not None and res["KO_over_WT"] < 0.15:
                        with kom:
                            restrict(kom, bounds)
                            cap = -kom.reactions.R1220.lower_bound
                            span = dc.fva(kom, ["R1220"], 1.0)["R1220"]
                            kom.reactions.R1220.lower_bound = -1000.0
                            k_open = dc.growth(kom)[0]
                        with wt:
                            restrict(wt, bounds)
                            wt.reactions.R1220.lower_bound = -1000.0
                            w_open = dc.growth(wt)[0]
                        record["lysine"].append({"config": config, "overlay": overlay, "caps": caps,
                                                 "variant": variant, "R1220_cap": cap, "R1220_fva_at_KO_opt": span,
                                                 "at_cap": all(abs(x + cap) < 1e-9 for x in span),
                                                 "KO_over_WT_lysine_open": None if not w_open or k_open is None
                                                 else k_open / w_open})
                if overlay and caps:
                    for variant, bounds in VARIANTS.items():
                        for label, m in (("WT", wt), ("KO", kom)):
                            with m:
                                restrict(m, bounds)
                                if dc.growth(m)[0]:
                                    record["fva"][f"{config}|{variant}|{label}"] = {
                                        rid: [round(dc.per100(x), 4) for x in v]
                                        for rid, v in dc.fva(m, list(FVA_RXNS), 1.0).items()}

    # C2 floor row (R825 ASSUMPTION + PDH + CS + G6PDH), overlay + caps on (Task D step 5).
    for config in dc.CONFIGS:
        wt, _ = dc.load(config)
        kom, _ = dc.load(config)
        dc.ko(kom)
        dc.r825_assumption(wt)
        dc.r825_assumption(kom)
        ref = next(r for r in ratios_ref["rows"] if r["config"] == config and r["r825"] == "assumption"
                   and r["floors"] == FLOOR_SET)
        for variant, bounds in VARIANTS.items():
            res = pair(wt, kom, lambda m, b=bounds: (add_floors(m, FLOOR_SET, fractions), restrict(m, b)))
            res.update(config=config, variant=variant, floors=FLOOR_SET, r825="assumption")
            record["floors"].append(res)
            if variant == "asis":
                ok = close_enough(res["KO_over_WT"], ref["KO_over_WT"]) and close_enough(res["WT"], ref["WT"])
                record["checks"].append({"what": f"ratios.json {config} R825 assumption PDH+CS+G6PDH",
                                         "expected": ref["KO_over_WT"], "got": res["KO_over_WT"], "ok": ok})
            print("floors", config, variant, res["WT"], res["KO_over_WT"], res["KO_status"], flush=True)

    # C3 WHAT-IF (R831 off, overlay + caps on).
    for config in dc.CONFIGS:
        wt, _ = dc.load(config)
        kom, _ = dc.load(config)
        dc.ko(kom)

        for label, opened in WHATIF.items():
            def setup(m, opened=opened):
                restrict(m, VARIANTS["off"])
                for rid in opened:
                    r = m.reactions.get_by_id(rid)
                    assert r.lower_bound == 0.0, (rid, r.bounds)
                    r.lower_bound = -1000.0
            res = pair(wt, kom, setup)
            res.update(config=config, whatif=label, label="WHAT-IF (diagnostic; not a proposal)")
            record["whatif"].append(res)
            print("whatif", config, label, res["WT"], res["KO_over_WT"], res["KO_status"], flush=True)

    record.update(lps=dc.Counter_.lps, elapsed_seconds=time.monotonic() - started,
                  all_checks_ok=all(c["ok"] for c in record["checks"]))
    dc.write_new(HERE / "complexI_r831.json", json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "complexI_r831.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["block", "config", "overlay", "caps", "variant_or_whatif", "WT", "WT_cost_vs_asis", "KO",
                    "KO_status", "KO_over_WT", "essential_1", "essential_5", "essential_10", "essential_15"])
        for r in record["grid"]:
            w.writerow(["grid", r["config"], r["overlay"], r["caps"], r["variant"], r["WT"], r["WT_cost"], r["KO"],
                        r["KO_status"], r["KO_over_WT"], *r["calls"].values()])
        for r in record["floors"]:
            w.writerow(["floors_R825_PDH_CS_G6PDH", r["config"], True, True, r["variant"], r["WT"], None, r["KO"],
                        r["KO_status"], r["KO_over_WT"], *r["calls"].values()])
        for r in record["whatif"]:
            w.writerow(["whatif_R831_off", r["config"], True, True, r["whatif"], r["WT"], None, r["KO"],
                        r["KO_status"], r["KO_over_WT"], *r["calls"].values()])
    print("checks ok", record["all_checks_ok"], "LPs", dc.Counter_.lps, "elapsed", round(record["elapsed_seconds"]))
    if not record["all_checks_ok"]:
        print("STOP: a reproduction check failed", [c for c in record["checks"] if not c["ok"]])


if __name__ == "__main__":
    main()
