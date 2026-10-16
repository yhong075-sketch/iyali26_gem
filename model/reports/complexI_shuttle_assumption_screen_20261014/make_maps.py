"""Maps for the (a)/(b) screen (drawing only, no LP).

map_matrix_nadh_outlet.png: benchmark context (as loaded, no overlay); panels WT base, complex I KO base, KO + (a)
R533 reductive closed, KO + (b) R831 forward closed, plus the lysine-cap sweep. Edge values are per 100 mol glucose
taken up in that panel, from one pFBA vector each (one optimum; judgments rest on the FVA in run.json). Dashed
grey = closed in that panel. Also writes map_edges.csv.

Run from the repository root with a Python that has matplotlib:
    IYALI26_RESEARCH_ROOT=... python -B model/reports/complexI_shuttle_assumption_screen_20261014/make_maps.py
"""
import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
ART = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve() / "artifacts/tasks" / HERE.name

NODES = {
    "MED": (0.8, 5.2, "medium\nlysine"), "LYSm": (0.8, 1.0, "lysine\nsynthesis (mito)"),
    "NADHm": (3.4, 1.0, "NADH (mito)"), "MALm": (3.4, 3.0, "malate (mito)"), "MALc": (6.4, 3.0, "malate (cyt)"),
    "NADHc": (6.4, 5.2, "NADH (cyt)"), "CoQ": (8.8, 1.0, "CoQ"), "BIO": (3.4, 5.2, "biomass\nlysine"),
}
# (label, from, to, {reaction: sign}, closed_in)
EDGES = [
    ("R1220 lysine uptake", "MED", "BIO", {"R1220": -1}, ()),
    ("R398 homoisocitrate DH", "LYSm", "NADHm", {"R398": 1}, ()),
    ("R533 matrix MDH, reductive", "NADHm", "MALm", {"R533": -1}, ("a",)),
    ("R831 malate out / citrate in", "MALm", "MALc", {"R831": 1}, ("b",)),
    ("R807 malate out / oxoadipate in", "MALm", "MALc", {"R807": 1}, ()),
    ("R534 cyt MDH, oxidative", "MALc", "NADHc", {"R534": 1}, ()),
    ("R2062 complex I", "NADHm", "CoQ", {"R2062": 1}, ("KO",)),
    ("R570 NDH2 (external)", "NADHc", "CoQ", {"R570": 1}, ()),
]
PANELS = [("bench_base", "WT", "WT, as supplied", ()), ("bench_base", "complex_I_KO", "complex I KO", ("KO",)),
          ("bench_a", "complex_I_KO", "complex I KO + (a)", ("KO", "a")),
          ("bench_b", "complex_I_KO", "complex I KO + (b)", ("KO", "b"))]


def main():
    fig, axes = plt.subplots(1, 5, figsize=(30, 6.5), gridspec_kw={"width_ratios": [1, 1, 1, 1, 0.8]})
    rows = []
    for ax, (run, state, title, closed) in zip(axes, PANELS):
        v = json.loads((ART / f"pfba_route_{run}.json").read_text())[state]
        uptake = -v["R1070"]
        ax.set_xlim(0, 9.8)
        ax.set_ylim(0, 6.2)
        ax.axis("off")
        ax.axhspan(0, 3.6, xmin=0, xmax=0.5, color="#eef3fb", zorder=0)
        ax.text(0.1, 3.75, "matrix", fontsize=9, color="#4a6fa5")
        ax.set_title(f"{title}\ngrowth {v['biomass_C']:.4f} h$^{{-1}}$, glucose {uptake:.3f}", fontsize=11)
        for key, (x, y, name) in NODES.items():
            ax.text(x, y, name, ha="center", va="center", fontsize=9,
                    bbox=dict(boxstyle="round", fc="white", ec="#555"), zorder=3)
        seen = {}
        for label, a, b, rxns, closed_in in EDGES:
            value = 100 * sum(s * v.get(r, 0.0) for r, s in rxns.items()) / uptake
            if label.startswith("R1220"):
                value = 100 * -v["R1220"] / uptake
            is_closed = any(c in closed for c in closed_in)
            offset = seen.get((a, b), 0)
            seen[(a, b)] = offset + 1
            (x1, y1), (x2, y2) = NODES[a][:2], NODES[b][:2]
            if value < -1e-9:
                (x1, y1), (x2, y2) = (x2, y2), (x1, y1)
            rad = 0.25 * offset
            ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), connectionstyle=f"arc3,rad={rad}",
                                         arrowstyle="-|>", mutation_scale=14, shrinkA=22, shrinkB=22,
                                         lw=1 if is_closed or abs(value) < 1e-6 else 1 + min(abs(value), 600) / 80,
                                         color="#aaa" if is_closed or abs(value) < 1e-6 else "#c0392b",
                                         linestyle="--" if is_closed else "-", zorder=2))
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2 + 0.35 * offset + 0.15
            ax.text(mx, my, f"{label.split(' ')[0]} {'closed' if is_closed else f'{abs(value):.1f}'}",
                    ha="center", fontsize=8, color="#333", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
            rows.append({"panel": title, "edge": label, "per_100_glucose": value, "closed": is_closed})
    sweep = json.loads((HERE / "lysine_sweep.json").read_text())["rows"]
    ax = axes[4]
    for variant, colour in (("base", "#555"), ("a", "#c0392b"), ("b", "#2471a3")):
        pts = [(r["lysine_factor"], r["KO_over_WT"]) for r in sweep if r["variant"] == variant]
        ax.plot(*zip(*pts), marker="o", color=colour, label={"base": "as supplied", "a": "(a) R533 rev closed",
                                                              "b": "(b) R831 fwd closed"}[variant],
                linestyle="--" if variant == "b" else "-")
    ax.axhline(0.15, color="#999", lw=0.8)
    ax.text(9.8, 0.16, "15 % cutoff", ha="right", fontsize=8, color="#666")
    ax.axvline(1.0, color="#999", lw=0.8, linestyle=":")
    ax.text(1.05, 0.7, "Task A cap\n(0.00595)", fontsize=8, color="#666")
    ax.set_xlabel("lysine uptake cap R1220 (× Task A medium)")
    ax.set_ylabel("complex I KO / WT")
    ax.set_title("complex I KO/WT vs lysine cap\n(benchmark context)", fontsize=11)
    ax.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(HERE / "map_matrix_nadh_outlet.png", dpi=110)
    with (HERE / "map_edges.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
