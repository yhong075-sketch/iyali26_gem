"""Map for the R831 rerun (drawing only, no LP).

map_r831_matrix_malate.png: config A (Task A medium, overlay + caps), complex I KO for R831 asis / off / rev0 and WT
asis, from the FVA ranges at the optimum in complexI_r831.json (per 100 mol glucose at the 2.414 bound, Task D's
per100). Each edge shows [min, max]; values above 30,000 are shown as "loop", because the R831/R1386/R817 cycle and the separate
R833 cycle reach the 1,000 flux bound. Dashed grey = closed in that panel. fwd0 is identical to off (complexI_r831.csv).
Also writes map_edges.csv.

Run: python -B model/reports/r831_ctp1_mode_rerun_20261014/make_maps.py   (any Python with matplotlib)
"""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
NODES = {"LYS": (0.9, 1.0, "lysine synthesis\n(R398 HIDH)"), "NADHm": (3.2, 1.0, "NADH (mito)"),
         "MALm": (3.2, 3.2, "malate / OAA\n(mito)"), "CITm": (1.2, 4.6, "citrate (mito)"),
         "MALc": (6.4, 3.2, "malate (cyt)"), "CITc": (6.4, 4.6, "citrate (cyt)"), "OGc": (6.4, 1.6, "NADH (cyt)"),
         "CoQ": (8.6, 1.0, "CoQ / NDH2")}
# (label, from, to, reaction, sign: +1 = forward is from->to)
EDGES = [("R398", "LYS", "NADHm", "R398", 1),
         ("R533 MDH (− = reductive)", "MALm", "NADHm", "R533", 1),
         ("R831 malate out / citrate in", "MALm", "MALc", "R831", 1),
         ("R807 malate in / 2-OA out", "MALc", "MALm", "R807", 1),
         ("R1386 malate in / 2-OG out", "MALc", "MALm", "R1386", 1),
         ("R902 malate in / Pi out", "MALc", "MALm", "R902", 1),
         ("R817 citrate out / 2-OG in", "CITm", "CITc", "R817", 1),
         ("R2062 complex I", "NADHm", "CoQ", "R2062", 1),
         ("R570 NDH2", "OGc", "CoQ", "R570", 1)]
PANELS = [("A|asis|WT", "WT, R831 as supplied", ()), ("A|asis|KO", "complex I KO, R831 as supplied", ("R2062",)),
          ("A|off|KO", "complex I KO, R831 off (= fwd0)", ("R2062", "R831")),
          ("A|rev0|KO", "complex I KO, R831 rev0", ("R2062",))]


def fmt(span):
    lo, hi = span
    f = lambda x: "loop" if abs(x) > 30000 else f"{x:.1f}"  # noqa: E731
    return f"[{f(lo)}, {f(hi)}]" if abs(hi - lo) > 0.05 else f"{lo:.1f}"


def main():
    d = json.loads((HERE / "complexI_r831.json").read_text())
    growth = {f"{r['config']}|{r['variant']}": r for r in d["grid"] if r["overlay"] and r["caps"]}
    fig, axes = plt.subplots(1, 4, figsize=(26, 6.5))
    rows = []
    for ax, (key, title, closed) in zip(axes, PANELS):
        fva = d["fva"][key]
        g = growth["|".join(key.split("|")[:2])]
        ax.set_xlim(0, 9.6)
        ax.set_ylim(0, 5.4)
        ax.axis("off")
        ax.axvspan(0, 4.6, color="#eef3fb", zorder=0)
        ax.text(0.15, 5.15, "matrix", color="#4a6fa5", fontsize=9)
        ax.text(4.8, 5.15, "cytosol", color="#777", fontsize=9)
        grow = g["WT"] if key.endswith("WT") else g["KO"]
        extra = "" if key.endswith("WT") else f", KO/WT {g['KO_over_WT']:.3f}"
        ax.set_title(f"{title}\ngrowth {grow:.4f} h$^{{-1}}${extra}", fontsize=11)
        for name, (x, y, text) in NODES.items():
            ax.text(x, y, text, ha="center", va="center", fontsize=8.5,
                    bbox=dict(boxstyle="round", fc="white", ec="#555"), zorder=3)
        seen = {}
        for label, a, b, rid, _ in EDGES:
            span = fva[rid]
            is_closed = rid in closed
            carries = max(abs(x) for x in span) > 1e-6
            off = seen.get(frozenset((a, b)), 0)
            seen[frozenset((a, b))] = off + 1
            (x1, y1), (x2, y2) = NODES[a][:2], NODES[b][:2]
            if span[1] <= 1e-9 and span[0] < -1e-9:
                (x1, y1), (x2, y2) = (x2, y2), (x1, y1)
            ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), connectionstyle=f"arc3,rad={0.18 * off}",
                                         arrowstyle="-|>", mutation_scale=13, shrinkA=24, shrinkB=24,
                                         lw=2.2 if carries and not is_closed else 1,
                                         color="#aaa" if is_closed or not carries else "#c0392b",
                                         linestyle="--" if is_closed else "-", zorder=2))
            mx, my = (x1 + x2) / 2, (y1 + y2) / 2 + 0.28 * off
            ax.text(mx, my, f"{rid} {'closed' if is_closed else fmt(span)}", ha="center", fontsize=7.5, zorder=4,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
            rows.append({"panel": title, "edge": label, "reaction": rid, "fva_min_per100": span[0],
                         "fva_max_per100": span[1], "closed": is_closed})
    fig.text(0.5, 0.01, "Config A (Task A medium, PO1f, holistic overlay + 13C caps; diagnostic). FVA at the optimum, "
             "per 100 glucose at the 2.414 bound. 'loop' = flux at the 1,000 bound from the R831/R1386/R817 or R833 cycle. "
             "", ha="center", fontsize=8.5, color="#555")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(HERE / "map_r831_matrix_malate.png", dpi=110)
    with (HERE / "map_edges.csv").open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
