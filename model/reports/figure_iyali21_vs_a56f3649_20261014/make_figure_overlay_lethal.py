"""Companion to make_figure.py: both models under the a56f3649 panel's conditions.

Left: iYali21; right: a56f3649. Both under the in-memory diagnostic overlay (assumed caps R1384 <= 8,
R860 <= 1; adapted to iYali21 in screen_iyali21_overlay_w29_20261014), PO1f + CSM-Leu with glucose
2.414, labelled lethal view (no-solution knockouts counted as lethal). Same layout as the mock-up.

No optimization. Counts are read from the two overlay summary.json files (hashes checked).

Usage (from the repository root):
    python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure_overlay_lethal.py
"""
import json

from make_figure import HERE, REPORTS, SOURCES as MOCKUP_SOURCES, at_primary, draw, sha

SOURCES = {
    "iyali21_overlay": (REPORTS / "screen_iyali21_overlay_w29_20261014/summary.json",
                        "3394996746d0e47dea9a4297e705e58d9e45c0d0226970de222ccbb53794c92d"),
    "a56f3649_overlay": MOCKUP_SOURCES["a56f3649_overlay"],
}


def panel(title, summary, expected):
    views = summary["metrics_mtdna_out_of_scope"]
    m = at_primary(views["lethal_view"])
    unresolved = at_primary(views["benchmark_rule"])["unresolved"]
    # The counts reported to the user; stop if the source says otherwise.
    assert [m[k] for k in ("TP", "FN", "FP", "TN")] == expected
    return dict(title=title, metrics=m,
                condition=f"Diagnostic overlay (assumed caps) · {unresolved} no-solution KOs counted as lethal")


def main():
    for key, (path, digest) in SOURCES.items():
        assert sha(path) == digest, f"Source summary changed: {key}"
    load = lambda key: json.loads(SOURCES[key][0].read_text())
    panels = [panel("iYali21 (iyli21.xml, 6974b758)", load("iyali21_overlay"), [74, 246, 30, 692]),
              panel("iYali26 newest candidate (a56f3649)", load("a56f3649_overlay"), [102, 220, 40, 695])]
    stem = HERE / "confusion_matrices_15pct_overlay_lethal"
    draw(panels, stem)
    record = dict(sources={k: dict(path=str(p.relative_to(REPORTS.parent)), sha256=d) for k, (p, d) in SOURCES.items()},
                  panels=[dict(title=p["title"], condition=p["condition"],
                               counts={k: p["metrics"][k] for k in ("TP", "FN", "FP", "TN")}) for p in panels],
                  figure={stem.with_suffix(s).name: sha(stem.with_suffix(s)) for s in (".png", ".pdf")})
    (HERE / "figure_overlay_lethal.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["panels"], indent=1))


if __name__ == "__main__":
    main()
