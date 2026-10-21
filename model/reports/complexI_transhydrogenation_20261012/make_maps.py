"""Task D maps (diagnostic; values per 100 mol glucose actually taken up in each panel, from pFBA vectors, one optimum
per panel; the uptake is printed in each info box).

1. map_complexI_cutsets.png: rows A / V2; columns WT, complex I KO, KO + G5 (R533 reverse closed), KO + G3 (R831
   forward closed).
2. map_13c_floors.png: WT and KO under PDH + CS + G6PDH floors with the R825 ASSUMPTION (A, V2).
3. map_discovery_rounds.png: the complex I KO (A) at discovery rounds 0-4.
Also writes map_edges.csv (every drawn edge and its value per panel). Judgments rest on the FVA in the JSON outputs.
Dashed grey edges are closed in that panel (every reaction of the edge, or the stated direction). The PNGs and
map_edges.csv are drawing outputs regenerated from map_vectors_pfba.json; a re-run overwrites them.

Drawing only (no LP). Run from the repository root with a Python that has matplotlib:
    IYALI26_RESEARCH_ROOT=... python -B model/reports/complexI_transhydrogenation_20261012/make_maps.py
"""
import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
ART = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve() / "artifacts/tasks" / HERE.name
UPTAKE = 2.414

NODES = {
    "GLC": (0.7, 9.7, "glucose"), "G6P": (0.7, 8.3, "G6P"), "OXPPP": (2.5, 9.7, "oxPPP"),
    "NADPHc": (4.7, 9.7, "NADPH (cyt)"), "MTHFD": (3.9, 8.0, "MTHFD pair"), "GDH": (5.5, 8.0, "GDH pair"),
    "OTHER": (7.0, 8.0, "other cyt\nNADH routes"), "NADHc": (8.4, 9.7, "NADH (cyt)"), "NDH2": (9.5, 7.8, "NDH2\n(external)"),
    "PYRc": (0.7, 6.6, "pyruvate (cyt)"), "ACAL": (2.5, 7.3, "acetaldehyde /\nacetate"), "GLX": (3.7, 5.9, "glyoxylate\nshunt"),
    "MALc": (6.3, 6.2, "malate (cyt)"), "OUT": (1.6, 5.5, "secreted"), "CoQ": (8.9, 5.0, "CoQ"),
    "CIII": (9.8, 3.3, "complex\nIII / IV"), "PYRm": (0.7, 3.8, "pyruvate (mito)"), "ACCOAm": (2.4, 2.6, "acetyl-CoA (mito)"),
    "TCA": (4.4, 1.6, "TCA cycle"), "MALm": (6.3, 3.8, "malate (mito)"), "NADHm": (6.7, 1.6, "NADH (mito)"),
    "CI": (8.3, 2.3, "complex I"), "LYS": (4.4, 0.2, "lysine synthesis"),
}
EDGES = [
    ("R1070 glucose in (=100)", "GLC", "G6P", {"R1070": -1}),
    ("R325 G6PDH", "G6P", "OXPPP", {"R325": 1}),
    ("R325+R639 NADPH", "OXPPP", "NADPHc", {"R325": 1, "R639": 1}),
    ("R694 PYK", "G6P", "PYRc", {"R694": 1}),
    ("R344 GAPDH", "G6P", "NADHc", {"R344": 1}),
    ("R691 PDC", "PYRc", "ACAL", {"R691": 1}),
    ("R134 ALD (NADP)", "ACAL", "NADPHc", {"R134": 1}),
    ("R1378 ALD (NAD)", "ACAL", "NADHc", {"R1378": 1}),
    ("R550 rev (NADP)", "NADPHc", "MTHFD", {"R550": -1}),
    ("R549 (NAD)", "MTHFD", "NADHc", {"R549": 1}),
    ("R330 (NADP)", "NADPHc", "GDH", {"R330": 1}),
    ("R329 (NAD)", "GDH", "NADHc", {"R329": 1}),
    ("R641 PGDH, R1931, R2119, ...", "OTHER", "NADHc",
     {"R641": 1, "R1931": 1, "R2119": 1, "R2064": 1, "R133": 1, "R2284": 1, "R45": 1}),
    ("R534 cyt MDH", "MALc", "NADHc", {"R534": 1}),
    ("R570 NDH2", "NADHc", "NDH2", {"R570": 1}),
    ("R570 to CoQ", "NDH2", "CoQ", {"R570": 1}),
    ("R90 ACS", "ACAL", "GLX", {"R90": 1}),
    ("R536+R537 MS", "GLX", "MALc", {"R536": 1, "R537": 1}),
    ("R1148 glyoxylate", "GLX", "OUT", {"R1148": 1}),
    ("R1310 pyruvate", "PYRc", "OUT", {"R1310": 1}),
    ("R1311-R842", "PYRc", "PYRm", {"R1311": 1, "R842": -1}),
    ("R693 PDH", "PYRm", "ACCOAm", {"R693": 1}),
    ("R693 PDH NADH", "ACCOAm", "NADHm", {"R693": 1}),
    ("R236 CS", "ACCOAm", "TCA", {"R236": 1}),
    ("R487+R361 NADH", "TCA", "NADHm", {"R487": 1, "R361": 1}),
    ("R740 SDH", "TCA", "CoQ", {"R740": 1}),
    ("R533 matrix MDH", "MALm", "NADHm", {"R533": 1}),
    ("R831 malate out", "MALm", "MALc", {"R831": 1}),
    ("R902+R1386+R807 malate in", "MALc", "MALm", {"R902": 1, "R1386": 1, "R807": 1}),
    ("R398 HICDH", "LYS", "NADHm", {"R398": 1}),
    ("R2062 complex I", "NADHm", "CI", {"R2062": 1}),
    ("R2062 to CoQ", "CI", "CoQ", {"R2062": 1}),
    ("R305 CIII", "CoQ", "CIII", {"R305": 1}),
]


def value(vec, coeffs):
    """Per 100 mol glucose actually taken up in that vector (falls back to the 2.414 bound at zero uptake)."""
    uptake = -vec.get("R1070", 0.0)
    return 100.0 * sum(c * vec.get(r, 0.0) for r, c in coeffs.items()) / (uptake if uptake > 1e-9 else UPTAKE)


def draw(ax, vec, title, info, removed=()):
    ax.add_patch(Rectangle((-0.2, 5.0), 10.4, 5.4, color="#eef4fb", zorder=0))
    ax.add_patch(Rectangle((-0.2, -0.4), 10.4, 5.4, color="#fbf1e6", zorder=0))
    ax.text(-0.1, 10.3, "cytosol", fontsize=7, color="#456", va="top")
    ax.text(-0.1, 4.9, "matrix", fontsize=7, color="#754", va="top")
    for key, (x, y, label) in NODES.items():
        ax.text(x, y, label, ha="center", va="center", fontsize=6.3, zorder=3,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#555", lw=0.6))
    rows = []
    for label, a, b, coeffs in EDGES:
        v = None if vec is None else value(vec, coeffs)
        is_removed = all(r in removed for r in coeffs)
        src, dst = (a, b) if v is None or v >= -1e-9 else (b, a)
        width = 0.6 if v is None else 0.6 + min(abs(v), 300) / 50
        color = "#bbb" if v is None or abs(v) < 0.05 else ("#1f6fb4" if v > 0 else "#c0392b")
        if is_removed:
            color = "#888"
        (x1, y1, _), (x2, y2, _) = NODES[src], NODES[dst]
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=8, lw=width, color=color,
                                     connectionstyle="arc3,rad=0.12", shrinkA=13, shrinkB=13, zorder=2,
                                     linestyle="--" if is_removed else "-"))
        mx, my = (x1 + x2) / 2 + 0.06 * (y2 - y1), (y1 + y2) / 2 - 0.06 * (x2 - x1)
        txt = f"{label}{' [closed]' if is_removed else ''}\n{'n/a' if v is None else f'{abs(v):.1f}'}"
        ax.text(mx, my, txt, fontsize=4.8, ha="center", va="center", zorder=4,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
        rows.append((label, v))
    ax.set_title(title, fontsize=8)
    ax.text(10.1, -0.35, info, fontsize=5.6, ha="right", va="bottom", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", fc="#fffbe0", ec="#aa8"))
    ax.set_xlim(-0.3, 10.3)
    ax.set_ylim(-0.5, 10.5)
    ax.axis("off")
    return rows


def info_text(meta, vec, wt_growth, extra):
    g = meta["growth"]
    ratio = "" if wt_growth is None or g is None else f", KO/WT {g / wt_growth:.3f}"
    text = f"growth {g:.4f}{ratio}" if g is not None else f"NO SOLUTION ({meta['status']})"
    if vec:
        text += f"; glucose uptake {-vec.get('R1070', 0.0):.3f}"
        text += (f"\ncomplex I {value(vec, {'R2062': 1}):.1f}, NDH2 {value(vec, {'R570': 1}):.1f}, "
                 f"R533 {value(vec, {'R533': 1}):.1f}")
    return text + ("\n" + extra if extra else "")


def main():
    data = json.loads((ART / "maps/map_vectors_pfba.json").read_text())
    vectors, meta = data["vectors"], data["meta"]
    definition = json.loads((HERE / "parts_definition.json").read_text())
    edge_rows = []

    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    columns = (("WT", "WT", (), ""), ("KO", "complex I KO", ("R2062",), "R2062 closed"),
               ("KO+G5", "KO + G5: R533 reverse closed", ("R2062", "R533"), "R533 lb 0 (forward kept)"),
               ("KO+G3", "KO + G3: R831 forward closed", ("R2062", "R831"), "R831 ub 0 (reverse kept)"))
    for i, config in enumerate(("A", "V2")):
        wt = meta[f"{config}|WT"]["growth"]
        for j, (suffix, title, removed, extra) in enumerate(columns):
            key = f"{config}|{suffix}"
            medium = "Task A medium" if config == "A" else "W29 V2 defined"
            info = info_text(meta[key], vectors[key], None if suffix == "WT" else wt, f"{medium}; overlay + 13C caps; {extra}")
            rows = draw(axes[i, j], vectors[key], f"{config}: {title}", info, removed)
            edge_rows += [("cutsets", key, lbl, v) for lbl, v in rows]
    fig.suptitle("Task D: complex I escape and the two single-group cut sets (pFBA, per 100 glucose; blue = drawn "
                 "direction, red = reversed; dashed = closed)", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "map_complexI_cutsets.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(24, 6.5))
    for j, key in enumerate(("A|WT|floors", "A|KO|floors", "V2|WT|floors", "V2|KO|floors")):
        config = key.split("|")[0]
        wt = meta[f"{config}|WT|floors"]["growth"]
        state = key.split("|")[1]
        info = info_text(meta[key], vectors[key], None if state == "WT" else wt,
                         "ASSUMPTION: PDH, CS, G6PDH >= 13C lower bounds; R825 reversible")
        rows = draw(axes[j], vectors[key], f"{config}: {state} under 13C floors", info,
                    ("R2062",) if state == "KO" else ())
        edge_rows += [("floors", key, lbl, v) for lbl, v in rows]
    fig.suptitle("Task D step 5: complex I KO under the W29 13C ratio floors (ASSUMPTION overlay; pFBA, per 100 glucose)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "map_13c_floors.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 5, figsize=(30, 6.5))
    for j, r in enumerate((0, 1, 2, 3, 4)):
        key = "A|KO" if r == 0 else f"A|KO|round{r}"
        wt = meta["A|WT"]["growth"] if r == 0 else json.loads((HERE / f"discovery_round{r}.json").read_text())["WT"]
        closed = sorted({rid for p in definition["parts"] if p["round_added"] < r for rid in p["closures"]})
        info = info_text(meta[key], vectors[key], wt, f"parts closed: {len(closed)} reactions")
        rows = draw(axes[j], vectors[key], f"A: complex I KO, discovery round {r}", info, ["R2062"] + closed)
        edge_rows += [("discovery", key, lbl, v) for lbl, v in rows]
    fig.suptitle("Task D discovery: where the complex I KO escape moves as parts are closed (configuration A; pFBA, "
                 "per 100 glucose; KO/WT against the WT with the same closures)", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "map_discovery_rounds.png", dpi=130)
    plt.close(fig)

    with (HERE / "map_edges.csv").open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["map", "panel", "edge", "value_per_100_glucose"])
        for r in edge_rows:
            w.writerow([r[0], r[1], r[2], None if r[3] is None else round(r[3], 3)])
    print("maps written")


if __name__ == "__main__":
    main()
