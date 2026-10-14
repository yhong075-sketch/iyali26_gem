"""Task C maps (diagnostic; values per 100 glucose from pFBA vectors, i.e. one optimum per panel).

1. map_c2_acetyl_routes.png: matrix acetyl-unit routes under PDH >= 0.866629 x uptake, V2 + overlay;
   rows = caps {none, module}; columns = WT as supplied, WT R825 reversible, complex I KO as supplied,
   complex I KO R825 reversible. Vectors from map_vectors.py
   ($IYALI26_RESEARCH_ROOT/artifacts/tasks/<task>/maps/map_vectors_pfba.json).
2. map_c1_citrate_malate_fit.png: elastic-fit representative points (V2, overlay, 22 ranges; pFBA on the
   optimal face, from c2/pfba_vectors.json), R825 as supplied vs reversible, with MFA 95 % intervals.
Also writes map_edges.csv (every drawn edge and its value per panel). FVA ranges for judgments are in
c2_runs.json / elastic_fit.json, not in these maps.

Drawing only (no LP). Run from the repository root with a Python that has matplotlib (imports nothing from the platform):
    IYALI26_RESEARCH_ROOT=... python -B model/reports/w29_13c_consistency_20261012/make_maps.py
"""
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

import os  # noqa: E402
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
ART = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve() / "artifacts/tasks" / HERE.name
UPTAKE = 2.414
ISOPRENOID_EXCHANGES = ("R945", "R1106", "R1108", "R1124", "R1235", "R1369", "R1388", "R1390")  # c2_acetyl_module.py

# Node positions: cytosol y > 5, matrix y < 5.
C2_NODES = {
    "PYRc": (0.8, 8.6, "pyruvate"), "ACALc": (2.6, 9.6, "acetaldehyde"), "ACc": (5.6, 9.6, "acetate"),
    "ACCOAc": (4.6, 7.6, "acetyl-CoA"), "ACARc": (6.8, 7.6, "acetylcarnitine"), "CITc": (2.6, 6.4, "citrate"),
    "MS": (4.6, 5.9, "malate synthase\n(glyoxylate shunt)"), "HMGc": (8.6, 6.3, "HMG-CoA"),
    "MEVc": (8.6, 8.1, "mevalonate ->\nisoprenoids"), "EXiso": (8.6, 9.7, "isoprenoid/sterol\nexchanges (8)"),
    "PYRm": (0.8, 3.4, "pyruvate"), "ACCOAm": (4.6, 2.6, "acetyl-CoA"), "CITm": (2.6, 1.0, "citrate"),
    "ACARm": (6.8, 2.6, "acetylcarnitine"), "AACOAm": (6.8, 0.9, "acetoacetyl-CoA"), "HMGm": (8.6, 0.9, "HMG-CoA"),
    "ACALm": (2.6, 4.3, "acetaldehyde"), "ACm": (5.6, 4.3, "acetate"),
}
# (label, from, to, {reaction: coefficient}); value = sum(coef * flux); negative values are drawn reversed.
C2_EDGES = [
    ("R691 PDC", "PYRc", "ACALc", {"R691": 1}),
    ("R1378+R134 ALD", "ACALc", "ACc", {"R1378": 1, "R134": 1}),
    ("R90 ACS", "ACc", "ACCOAc", {"R90": 1}),
    ("R997", "ACALm", "ACALc", {"R997": 1}),
    ("R135+R136+R2070 ALD", "ACALm", "ACm", {"R135": 1, "R136": 1, "R2070": 1}),
    ("R1708 SCACT", "ACm", "ACCOAm", {"R1708": 1}),
    ("R1311-R842", "PYRc", "PYRm", {"R1311": 1, "R842": -1}),
    ("R693 PDH", "PYRm", "ACCOAm", {"R693": 1}),
    ("R236 CS", "ACCOAm", "CITm", {"R236": 1}),
    ("cit export net", "CITm", "CITc", {"R817": 1, "R831": -1, "R832": -1, "R833": -1}),
    ("R1894 ACL", "CITc", "ACCOAc", {"R1894": 1}),
    ("R191 CAT", "ACCOAc", "ACARc", {"R191": 1}),
    ("R825+R1275", "ACARc", "ACARm", {"R825": 1, "R1275": 1}),
    ("R193 CAT", "ACARm", "ACCOAm", {"R193": 1}),
    ("R83 thiolase (x2 AcCoA)", "ACCOAm", "AACOAm", {"R83": 1}),
    ("R412 HMGS", "AACOAm", "HMGm", {"R412": 1}),
    ("R1172", "HMGc", "HMGm", {"R1172": 1}),
    ("R410 HMGR", "HMGc", "MEVc", {"R410": 1}),
    ("isoprenoid exch (sum)", "MEVc", "EXiso", {r: 1 for r in ISOPRENOID_EXCHANGES}),
    ("R536+R537 MS", "ACCOAc", "MS", {"R536": 1, "R537": 1}),
]
C1_NODES = {
    "PYRc": (0.8, 8.8, "pyruvate"), "OAAc": (2.8, 7.2, "oxaloacetate"), "MALc": (4.8, 8.8, "malate"),
    "CITc": (6.8, 7.0, "citrate"), "ACCOAc": (8.8, 8.8, "acetyl-CoA"), "ACARc": (8.8, 6.6, "acetylcarnitine"),
    "PYRm": (0.8, 2.2, "pyruvate"), "ACCOAm": (2.8, 1.0, "acetyl-CoA"), "OAAm": (4.8, 1.0, "oxaloacetate"),
    "MALm": (4.8, 3.2, "malate"), "CITm": (6.8, 2.2, "citrate"), "ACARm": (8.8, 2.4, "acetylcarnitine"),
}
C1_EDGES = [
    ("R1311-R842 [MFA 109.9-123.7]", "PYRc", "PYRm", {"R1311": 1, "R842": -1}),
    ("R693 PDH [86.7-100.7]", "PYRm", "ACCOAm", {"R693": 1}),
    ("R236 CS [28.2-51.2]", "ACCOAm", "CITm", {"R236": 1}),
    ("R533 MDH [28.2-51.2]", "MALm", "OAAm", {"R533": 1}),
    ("R690 PC", "PYRc", "OAAc", {"R690": 1}),
    ("R534 MDH cyt", "MALc", "OAAc", {"R534": 1}),
    ("R817 cit/2-OG", "CITm", "CITc", {"R817": 1}),
    ("R831 cit out / mal in", "CITm", "CITc", {"R831": -1}),
    ("R833 cit/icit", "CITm", "CITc", {"R833": -1}),
    ("R832 cit/PEP", "CITm", "CITc", {"R832": -1}),
    ("R831 (same flux): mal in", "MALc", "MALm", {"R831": -1}),
    ("R902+R1386+R807 mal in", "MALc", "MALm", {"R902": 1, "R1386": 1, "R807": 1}),
    ("R1894 ACL", "CITc", "ACCOAc", {"R1894": 1}),
    ("R191 CAT", "ACCOAc", "ACARc", {"R191": 1}),
    ("R825+R1275", "ACARc", "ACARm", {"R825": 1, "R1275": 1}),
    ("R193 CAT", "ACARm", "ACCOAm", {"R193": 1}),
]


def value(vec, coeffs):
    return 100.0 * sum(c * vec.get(r, 0.0) for r, c in coeffs.items()) / UPTAKE


def draw(ax, nodes, edges, vec, title, info, removed=()):
    ax.add_patch(Rectangle((-0.2, 5.0), 10.4, 5.4, color="#eef4fb", zorder=0))
    ax.add_patch(Rectangle((-0.2, -0.4), 10.4, 5.2, color="#fbf1e6", zorder=0))
    ax.text(-0.1, 10.2, "cytosol", fontsize=7, color="#456", va="top")
    ax.text(-0.1, 4.7, "matrix", fontsize=7, color="#754", va="top")
    for key, (x, y, label) in nodes.items():
        ax.text(x, y, label, ha="center", va="center", fontsize=6.5, zorder=3,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#555", lw=0.6))
    seen = {}
    rows = []
    for label, a, b, coeffs in edges:
        v = None if vec is None else value(vec, coeffs)
        pair = tuple(sorted((a, b)))
        k = seen.get(pair, 0)
        seen[pair] = k + 1
        rad = [0.0, 0.25, -0.25, 0.45, -0.45][k % 5]
        src, dst = (a, b) if v is None or v >= -1e-9 else (b, a)
        width = 0.6 if v is None else 0.6 + min(abs(v), 200) / 40
        color = "#bbb" if v is None or abs(v) < 0.05 else ("#1f6fb4" if v > 0 else "#c0392b")
        is_removed = any(r in label for r in removed)
        if is_removed:
            color = "#999"
        (x1, y1, _), (x2, y2, _) = nodes[src], nodes[dst]
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=8, lw=width, color=color,
                                     connectionstyle=f"arc3,rad={rad}", shrinkA=14, shrinkB=14, zorder=2,
                                     linestyle="--" if is_removed else "-"))
        mx, my = (x1 + x2) / 2 + rad * (y2 - y1) * 0.5, (y1 + y2) / 2 - rad * (x2 - x1) * 0.5
        my += [0.0, 0.7, -0.7, 1.4, -1.4][k % 5]  # keep labels of parallel edges apart
        txt = f"{label}\n{'n/a' if v is None else f'{abs(v):.1f}'}"
        ax.text(mx, my, txt, fontsize=5.2, ha="center", va="center", zorder=4,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
        rows.append((label, v))
    ax.set_title(title, fontsize=8)
    ax.text(10.1, -0.3, info, fontsize=5.8, ha="right", va="bottom", zorder=5,
            bbox=dict(boxstyle="round,pad=0.3", fc="#fffbe0", ec="#aa8"))
    ax.set_xlim(-0.3, 10.3)
    ax.set_ylim(-0.5, 10.5)
    ax.axis("off")
    return rows


def main():
    data = json.loads((ART / "maps/map_vectors_pfba.json").read_text())
    vectors, meta = data["vectors"], data["meta"]
    edge_rows = []
    fig, axes = plt.subplots(2, 4, figsize=(22, 11))
    titles = {"WT|asis": "WT, R825 as supplied (import only)", "WT|rev": "WT, R825 reversible (diagnostic)",
              "KO|asis": "complex I KO, R825 as supplied", "KO|rev": "complex I KO, R825 reversible"}
    for i, capped in enumerate(("none", "module")):
        for j, panel in enumerate(("WT|asis", "WT|rev", "KO|asis", "KO|rev")):
            key = f"caps={capped}|{panel}"
            vec, m = vectors[key], meta[key]
            info = (f"V2 + overlay; PDH >= 0.867 x uptake; 13C caps: {capped}\n"
                    f"growth {m['growth']:.4f}" if m["growth"] else
                    f"V2 + overlay; PDH >= 0.867 x uptake; 13C caps: {capped}\nNO SOLUTION ({m['status']}): "
                    f"unresolved, not lethal")
            if vec:
                info += (f"\ncomplex I R2062 {value(vec, {'R2062': 1}):.1f}, NDH2 R570 {value(vec, {'R570': 1}):.1f}"
                         f"\nR825 lb {'-1000' if 'rev' in panel else '0 (export not allowed)'}")
            title = titles[panel] + (" -- complex I removed" if "KO" in panel else "")
            rows = draw(axes[i, j], C2_NODES, C2_EDGES, vec, title, info)
            edge_rows += [("c2", key, lbl, v) for lbl, v in rows]
    fig.suptitle("Task C2: matrix acetyl-unit export under the PDH ratio floor (pFBA, per 100 glucose; "
                 "blue = drawn direction, red = reversed; KO panels: R2062 bounds 0)", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "map_c2_acetyl_routes.png", dpi=170)
    plt.close(fig)

    fit_vectors = json.loads((ART / "c2/pfba_vectors.json").read_text())
    rng = json.loads((HERE / "elastic_fit.json").read_text())["ranges_per_100_glucose"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax, tag, title in ((axes[0], "asis", "Elastic fit, R825 as supplied (total violation 6.52)"),
                           (axes[1], "rev", "Elastic fit, R825 reversible (total violation 1.43)")):
        vec = fit_vectors[f"fit|{tag}"]
        info = (f"V2 + overlay; glucose 2.414; growth 0.2663 fixed\nnet citrate export "
                f"{value(vec, C2_EDGES[9][3]):.1f} [MFA {rng['citrate_export'][0]:.1f}-{rng['citrate_export'][1]:.1f}]"
                f"\nnet malate import {value(vec, {'R807': 1, 'R902': 1, 'R1386': 1, 'R831': -1}):.1f} "
                f"[MFA {rng['malate_import'][0]:.1f}-{rng['malate_import'][1]:.1f}]\n"
                f"one optimum (pFBA on the optimal face); ranges in elastic_fit.json")
        rows = draw(ax, C1_NODES, C1_EDGES, vec, title, info)
        edge_rows += [("c1", f"fit|{tag}", lbl, v) for lbl, v in rows]
    fig.suptitle("Task C1: citrate export and malate import at the elastic-fit optimum (per 100 glucose)", fontsize=10)
    fig.tight_layout()
    fig.savefig(HERE / "map_c1_citrate_malate_fit.png", dpi=170)
    plt.close(fig)
    with (HERE / "map_edges.csv").open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["map", "panel", "edge", "value_per_100_glucose"])
        for r in edge_rows:
            w.writerow([r[0], r[1], r[2], None if r[3] is None else round(r[3], 3)])
    print("maps written")


if __name__ == "__main__":
    main()
