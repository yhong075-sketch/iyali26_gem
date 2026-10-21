"""C2 follow-up (DIAGNOSTIC, added after the run; not pre-registered).

The pre-registered lysine check reported KO/WT with lysine uptake R1220 opened to 1000 in both WT and KO. At 1000,
lysine also becomes a carbon source, so WT growth itself rises and that ratio is not interpretable. This records
absolute growth for WT and KO with R1220 at its cap, at 2x and 5x the cap, and at 1000, for asis / off / fwd0,
config A (overlay + caps on and both off), using the same Task D helpers as complexI_r831.py.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/r831_ctp1_mode_rerun_20261014/lysine_detail.py
Writes lysine_detail.json and lysine_detail.csv here. Never overwrites.
"""
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import complexI_r831 as c2  # noqa: E402  (pins Task D code on import of its helpers)

FACTORS = (1.0, 2.0, 5.0, None)  # None = lb -1000


def main():
    for path, digest in c2.CODE_PINS.items():
        assert c2.sha(path) == digest
    c2.dc.check_pins()
    rows = []
    for overlay_caps in (True, False):
        wt = c2.load_variant("A", overlay_caps, overlay_caps)
        kom = c2.load_variant("A", overlay_caps, overlay_caps)
        c2.dc.ko(kom)
        cap = -wt.reactions.R1220.lower_bound
        for variant in ("asis", "off", "fwd0"):
            for factor in FACTORS:
                def setup(m, factor=factor, variant=variant):
                    c2.restrict(m, c2.VARIANTS[variant])
                    m.reactions.R1220.lower_bound = -1000.0 if factor is None else -cap * factor
                res = c2.pair(wt, kom, setup)
                rows.append({"config": "A", "overlay_and_caps": overlay_caps, "variant": variant,
                             "lysine_lb": "-1000" if factor is None else -cap * factor,
                             "WT": res["WT"], "KO": res["KO"], "KO_status": res["KO_status"],
                             "KO_over_WT": res["KO_over_WT"]})
                print(rows[-1], flush=True)
    (HERE / "lysine_detail.json").open("x").write(json.dumps({"lps": c2.dc.Counter_.lps, "rows": rows},
                                                             indent=1) + "\n")
    with (HERE / "lysine_detail.csv").open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
