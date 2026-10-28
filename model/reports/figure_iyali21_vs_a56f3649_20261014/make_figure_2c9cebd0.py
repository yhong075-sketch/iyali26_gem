"""iYali21 vs the newest candidate 2c9cebd0 (2026-10-23), both on the Task A medium + PO1f without the
overlay, benchmark rule (no-solution knockouts left out). Same layout as make_figure_like_for_like.py.

No optimization. Counts are read from the two screen summary.json files (hashes checked).

Usage (from the repository root):
    python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure_2c9cebd0.py
"""
import json

from make_figure import HERE, REPORTS, SOURCES as MOCKUP_SOURCES, at_primary, draw, sha

SOURCES = {
    "iyali21": MOCKUP_SOURCES["iyali21"],
    "2c9cebd0": (REPORTS / "screen_2c9cebd0_20261022/summary.json",
                 "d549108ec37aff1d44f195b1850ce3980f914a7e3fe2264eb546864625e04337"),
}


def main():
    for key, (path, digest) in SOURCES.items():
        assert sha(path) == digest, f"Source summary changed: {key}"
    m_left = at_primary(json.loads(SOURCES["iyali21"][0].read_text())["models"]["iyali21"]["metrics_mtdna_out_of_scope"])
    m_right = at_primary(json.loads(SOURCES["2c9cebd0"][0].read_text())["metrics_mtdna_out_of_scope"]["benchmark_rule"])
    # The counts reported to the user; stop if the sources say otherwise.
    assert [m_left[k] for k in ("TP", "FN", "FP", "TN")] == [58, 255, 25, 694]
    assert [m_right[k] for k in ("TP", "FN", "FP", "TN")] == [93, 220, 36, 697]
    panels = [
        dict(title="iYali21 (iyli21.xml, 6974b758)", metrics=m_left,
             condition=f"No overlay · {m_left['unresolved']} no-solution KOs left out"),
        dict(title="iYali26 newest candidate (2c9cebd0)", metrics=m_right,
             condition=f"No overlay · {m_right['unresolved']} no-solution KOs left out"),
    ]
    stem = HERE / "confusion_matrices_15pct_2c9cebd0"
    draw(panels, stem)
    record = dict(sources={k: dict(path=str(p.relative_to(REPORTS.parent)), sha256=d) for k, (p, d) in SOURCES.items()},
                  panels=[dict(title=p["title"], condition=p["condition"],
                               counts={k: p["metrics"][k] for k in ("TP", "FN", "FP", "TN")}) for p in panels],
                  figure={stem.with_suffix(s).name: sha(stem.with_suffix(s)) for s in (".png", ".pdf")})
    (HERE / "figure_2c9cebd0.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["panels"], indent=1))


if __name__ == "__main__":
    main()
