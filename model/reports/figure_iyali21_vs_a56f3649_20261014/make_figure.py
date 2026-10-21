"""Paired 15 % confusion matrices as in the user's mock-up (2026-10-14).

Left: iYali21, no overlay, benchmark rule (no-solution knockouts left out).
Right: a56f3649 under Task A's diagnostic overlay, labelled lethal view (no-solution knockouts
counted as lethal). The two panels use different conditions and scoring; each panel says which.

No optimization. Counts are read from the two committed summary.json files (hashes checked).

Usage (from the repository root):
    python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure.py
"""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
REPORTS = HERE.parent
PRIMARY = 0.15
SOURCES = {
    "iyali21": (REPORTS / "confusion_matrices_csmleu_w29_20261014/summary.json",
                "0472de4103f4e6a30d18065c673a9ed854842210fe80088842eec5b191e5f393"),
    "a56f3649_overlay": (REPORTS / "screen_a56f3649_overlay_w29_20261014/summary.json",
                         "68bf6d1daaaab6d1d725877811b241089ff134eba0b60d9b730237043f794f6b"),
}
# Reference sequential blue ramp, steps 100-700 (dataviz skill, references/palette.md).
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SURFACE, TEXT_PRIMARY, TEXT_SECONDARY = "#fcfcfb", "#0b0b0b", "#52514e"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def text_color(fill):
    # Pick the ink with the higher WCAG contrast against the cell fill.
    def luminance(hex_color):
        channel = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in to_rgb(hex_color)]
        return 0.2126 * channel[0] + 0.7152 * channel[1] + 0.0722 * channel[2]
    fill_l = luminance(fill)
    dark = (fill_l + 0.05) / (luminance(TEXT_PRIMARY) + 0.05)
    light = 1.05 / (fill_l + 0.05)
    return TEXT_PRIMARY if dark >= light else "#ffffff"


def at_primary(metrics):
    return next(m for m in metrics if m["cutoff"] == PRIMARY)


def load_panels():
    for key, (path, digest) in SOURCES.items():
        assert sha(path) == digest, f"Source summary changed: {key}"
    left = json.loads(SOURCES["iyali21"][0].read_text())["models"]["iyali21"]
    right = json.loads(SOURCES["a56f3649_overlay"][0].read_text())
    m_left = at_primary(left["metrics_mtdna_out_of_scope"])
    m_right = at_primary(right["metrics_mtdna_out_of_scope"]["lethal_view"])
    unresolved_right = at_primary(right["metrics_mtdna_out_of_scope"]["benchmark_rule"])["unresolved"]
    # The panels the user approved; stop if the sources say otherwise.
    assert [m_left[k] for k in ("TP", "FN", "FP", "TN")] == [58, 255, 25, 694]
    assert [m_right[k] for k in ("TP", "FN", "FP", "TN")] == [102, 220, 40, 695]
    return [
        dict(title="iYali21 (iyli21.xml, 6974b758)", metrics=m_left,
             condition=f"No overlay · {m_left['unresolved']} no-solution KOs left out"),
        dict(title="iYali26 newest candidate (a56f3649)", metrics=m_right,
             condition=f"Diagnostic overlay (assumed caps) · {unresolved_right} no-solution KOs counted as lethal"),
    ]


def draw(panels, path_stem):
    cmap = LinearSegmentedColormap.from_list("seq_blue", RAMP)
    fig = plt.figure(figsize=(12, 6.3), facecolor=SURFACE)
    for i, panel in enumerate(panels):
        m = panel["metrics"]
        ax = fig.add_axes([0.11 + i * 0.48, 0.10, 0.34, 0.66])
        ax.set_facecolor(SURFACE)
        cells = [[("TP", m["TP"]), ("FN", m["FN"])], [("FP", m["FP"]), ("TN", m["TN"])]]
        for row, pair in enumerate(cells):
            total = sum(v for _, v in pair)
            for col, (name, value) in enumerate(pair):
                share = value / total
                fill = matplotlib.colors.to_hex(cmap(share))
                ax.add_patch(Rectangle((col, 1 - row), 1, 1, facecolor=fill,
                                       edgecolor=SURFACE, linewidth=2))
                ink = text_color(fill)
                ax.text(col + 0.5, 1 - row + 0.58, f"{value:,}", ha="center", va="center",
                        fontsize=22, fontweight="bold", color=ink)
                ax.text(col + 0.5, 1 - row + 0.32, f"{name} · {share:.1%} of row",
                        ha="center", va="center", fontsize=10, color=ink)
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 2)
        ax.set_aspect("equal")
        positives, negatives = m["TP"] + m["FN"], m["FP"] + m["TN"]
        ax.set_xticks([0.5, 1.5], ["Essential", "Non-essential"])
        ax.set_yticks([1.5, 0.5], [f"Essential\n(n = {positives})", f"Non-essential\n(n = {negatives})"])
        ax.set_xlabel("Model prediction", color=TEXT_SECONDARY, fontsize=10)
        ax.set_ylabel("Experiment (consensus list)", color=TEXT_SECONDARY, fontsize=10)
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.tick_params(length=0, colors=TEXT_PRIMARY, labelsize=10)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(panel["title"], fontsize=13, fontweight="bold", color=TEXT_PRIMARY, pad=54)
        ax.annotate(panel["condition"], xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 38),
                    textcoords="offset points", ha="center", va="bottom", fontsize=9.5, color=TEXT_SECONDARY)
    footer = ("Both panels: Task A medium + PO1f (glucose uptake 2.414); static FBA; essential = KO/WT < 15 %; "
              "16 mtDNA genes out of scope; negatives = model genes not on the consensus list (user-defined class).")
    fig.text(0.04, 0.025, footer, ha="left", va="bottom", fontsize=8, color=TEXT_SECONDARY)
    for suffix in (".png", ".pdf"):
        fig.savefig(path_stem.with_suffix(suffix), dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main():
    panels = load_panels()
    stem = HERE / "confusion_matrices_15pct_paired"
    draw(panels, stem)
    record = dict(sources={k: dict(path=str(p.relative_to(REPORTS.parent)), sha256=d) for k, (p, d) in SOURCES.items()},
                  panels=[dict(title=p["title"], condition=p["condition"],
                               counts={k: p["metrics"][k] for k in ("TP", "FN", "FP", "TN")}) for p in panels],
                  figure={stem.with_suffix(s).name: sha(stem.with_suffix(s)) for s in (".png", ".pdf")})
    (HERE / "figure.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["panels"], indent=1))


if __name__ == "__main__":
    main()
