"""iYali21 vs the newest candidate b3b9fe22 (2026-10-21), both under the diagnostic overlay
(assumed caps) in the labelled lethal view, PO1f + CSM-Leu with glucose 2.414. Same layout as
make_figure_overlay_lethal.py, with b3b9fe22 in place of a56f3649.

No optimization. Counts are read from the two overlay summary.json files (hashes checked).

Usage (from the repository root):
    python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure_b3b9fe22_overlay_lethal.py
"""
import json

from make_figure import HERE, REPORTS, at_primary, draw, sha
from make_figure_overlay_lethal import SOURCES as A56_SOURCES, panel

SOURCES = {
    "iyali21_overlay": A56_SOURCES["iyali21_overlay"],
    "b3b9fe22": (REPORTS / "screen_b3b9fe22_overlay_w29_20261021/summary.json",
                 "7700fcaa209c6b11e7aeeb78356e76696d6bff25586f9dc347db14334464ef3d"),
}


def main():
    for key, (path, digest) in SOURCES.items():
        assert sha(path) == digest, f"Source summary changed: {key}"
    load = lambda key: json.loads(SOURCES[key][0].read_text())
    newest = load("b3b9fe22")["conditions"]["overlay"]
    panels = [panel("iYali21 (iyli21.xml, 6974b758)", load("iyali21_overlay"), [74, 246, 30, 692]),
              panel("iYali26 newest candidate (b3b9fe22)", newest, [113, 211, 43, 694])]
    stem = HERE / "confusion_matrices_15pct_b3b9fe22_overlay_lethal"
    draw(panels, stem)
    record = dict(sources={k: dict(path=str(p.relative_to(REPORTS.parent)), sha256=d) for k, (p, d) in SOURCES.items()},
                  panels=[dict(title=p["title"], condition=p["condition"],
                               counts={k: p["metrics"][k] for k in ("TP", "FN", "FP", "TN")}) for p in panels],
                  figure={stem.with_suffix(s).name: sha(stem.with_suffix(s)) for s in (".png", ".pdf")})
    (HERE / "figure_b3b9fe22_overlay_lethal.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["panels"], indent=1))


if __name__ == "__main__":
    main()
