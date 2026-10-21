"""Like-for-like companion to make_figure.py: both models under iYali21's conditions.

Left: iYali21; right: a56f3649. Both without the overlay, on PO1f + CSM-Leu with glucose 2.414,
scored with the benchmark rule (no-solution knockouts left out). Same layout as the mock-up figure.

No optimization. Counts are read from the committed confusion-matrix summary.json (hash checked).

Usage (from the repository root):
    python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure_like_for_like.py
"""
import json

from make_figure import HERE, REPORTS, SOURCES, at_primary, draw, sha

SOURCE = SOURCES["iyali21"]  # confusion_matrices_csmleu_w29_20261014/summary.json holds both models


def load_panels():
    path, digest = SOURCE
    assert sha(path) == digest, "Source summary changed"
    models = json.loads(path.read_text())["models"]
    m_left = at_primary(models["iyali21"]["metrics_mtdna_out_of_scope"])
    m_right = at_primary(models["a56f3649"]["metrics_mtdna_out_of_scope"])
    # The panels the user asked for; stop if the source says otherwise.
    assert [m_left[k] for k in ("TP", "FN", "FP", "TN")] == [58, 255, 25, 694]
    assert [m_right[k] for k in ("TP", "FN", "FP", "TN")] == [82, 229, 33, 697]
    return [
        dict(title="iYali21 (iyli21.xml, 6974b758)", metrics=m_left,
             condition=f"No overlay · {m_left['unresolved']} no-solution KOs left out"),
        dict(title="iYali26 newest candidate (a56f3649)", metrics=m_right,
             condition=f"No overlay · {m_right['unresolved']} no-solution KOs left out"),
    ]


def main():
    panels = load_panels()
    stem = HERE / "confusion_matrices_15pct_like_for_like"
    draw(panels, stem)
    record = dict(source=dict(path=str(SOURCE[0].relative_to(REPORTS.parent)), sha256=SOURCE[1]),
                  panels=[dict(title=p["title"], condition=p["condition"],
                               counts={k: p["metrics"][k] for k in ("TP", "FN", "FP", "TN")}) for p in panels],
                  figure={stem.with_suffix(s).name: sha(stem.with_suffix(s)) for s in (".png", ".pdf")})
    (HERE / "figure_like_for_like.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["panels"], indent=1))


if __name__ == "__main__":
    main()
