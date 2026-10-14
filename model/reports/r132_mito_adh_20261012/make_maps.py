"""Draw the Task B maps (PNG) from stored flux vectors; no solver calls.

Reads run.json (FVA ranges) from this folder and the pFBA vectors in
$IYALI26_RESEARCH_ROOT/artifacts/tasks/r132_mito_adh_20261012/runs/. Writes:
  map_r132_route.png        WT and complex I KO, baseline (before) and R132 removed (after)
  map_step5_glucose.png     glucose, as supplied vs R132 removed + 13C caps (two glycine-cleavage forms)
Values are flux per 100 mol glucose taken up, from one pFBA optimum; bracketed values are the FVA
range at the growth optimum where it was computed. One optimum does not prove a route is unique.

Run from the repository root with a Python that has matplotlib:
    python -B model/reports/r132_mito_adh_20261012/make_maps.py
"""
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
RUNS = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve() / "artifacts/tasks" / HERE.name / "runs"

# Node positions: cytosol y > 0.5, matrix y < 0.5.
NODES = {
    "Glc": (0.16, 0.93), "Pyr_c": (0.22, 0.80), "AcAld_c": (0.50, 0.80), "EtOH_c": (0.74, 0.80),
    "EtOH_out": (0.95, 0.93), "OAA_c": (0.08, 0.66), "Mal_c": (0.36, 0.66), "NADH_c": (0.84, 0.64),
    "Pyr_m": (0.22, 0.36), "AcAld_m": (0.56, 0.36), "EtOH_m": (0.78, 0.36), "Mal_m": (0.36, 0.16),
    "OAA_m": (0.08, 0.16), "NADH_m": (0.84, 0.14), "Ac_m": (0.56, 0.05), "AcCoA_m": (0.06, 0.30),
    "Q": (0.98, 0.40), "AKG_m": (0.09, 0.03), "SucCoA_m": (0.45, 0.03),
}
LABELS = {"Glc": "glucose", "Pyr_c": "pyruvate", "AcAld_c": "acetaldehyde", "EtOH_c": "ethanol", "EtOH_out": "ethanol (ext)",
          "OAA_c": "OAA", "Mal_c": "malate", "NADH_c": "NADH (cyt)", "Pyr_m": "pyruvate", "AcAld_m": "acetaldehyde",
          "EtOH_m": "ethanol", "Mal_m": "malate", "OAA_m": "OAA", "NADH_m": "NADH (matrix)", "Ac_m": "acetate",
          "AcCoA_m": "acetyl-CoA", "Q": "CoQ pool",
          "AKG_m": "2-oxoglutarate", "SucCoA_m": "succinyl-CoA"}
# Edge: (label, source, target, flux function, fva ids or None). Positive flux = source -> target.
EDGES = [
    ("R691 PDC", "Pyr_c", "AcAld_c", lambda f: f["R691"], ("R691",)),
    ("R130 + R1377 + R_ALCDH_nadp_hi (net oxidative)", "EtOH_c", "AcAld_c",
     lambda f: f["R130"] - f["R1377"] + f.get("R_ALCDH_nadp_hi", 0), ("R130", "R1377", "R_ALCDH_nadp_hi")),
    ("R1113 export", "EtOH_c", "EtOH_out", lambda f: -f["R1113"], ("R1113",)),
    ("R1114", "EtOH_c", "EtOH_m", lambda f: f["R1114"], ("R1114",)),
    ("R132 (reductive +)", "AcAld_m", "EtOH_m", lambda f: f.get("R132"), ("R132",)),
    ("R997", "AcAld_m", "AcAld_c", lambda f: f["R997"], ("R997",)),
    ("R135 + R2070 + R136", "AcAld_m", "Ac_m", lambda f: f["R135"] + f["R2070"] + f["R136"], ("R135", "R2070", "R136")),
    ("R690 PC", "Pyr_c", "OAA_c", lambda f: f["R690"], ("R690",)),
    ("R534 + R_MDHy (rev = OAA->malate)", "Mal_c", "OAA_c", lambda f: f["R534"] + f.get("R_MDHy", 0), ("R534", "R_MDHy")),
    ("net malate import (R807+R902+R1386-R831)", "Mal_c", "Mal_m",
     lambda f: f["R807"] + f["R902"] + f["R1386"] - f["R831"], ("R807", "R831", "R902", "R1386")),
    ("R1311", "Pyr_c", "Pyr_m", lambda f: f["R1311"], ("R1311",), -0.3),
    ("R533 mMDH", "Mal_m", "OAA_m", lambda f: f["R533"], ("R533",)),
    ("R538 ME", "Mal_m", "Pyr_m", lambda f: f["R538"], ("R538",)),
    ("R693 PDH", "Pyr_m", "AcCoA_m", lambda f: f["R693"], ("R693",)),
    ("R570 NDH2", "NADH_c", "Q", lambda f: f["R570"], ("R570",), -0.25),
    ("R2062 complex I", "NADH_m", "Q", lambda f: f["R2062"], ("R2062",)),
    ("R614 2-OGDH (E3 = R361)", "AKG_m", "SucCoA_m", lambda f: f["R614"], ("R614",), 0.0),
    ("R348/R347 G3P shuttle", "NADH_c", "Q", lambda f: f["R347"], ("R348", "R347"), 0.35),
]


def per100(vector):
    uptake = -vector["R1070"]
    return {k: 100 * v / uptake for k, v in vector.items()}, uptake


def fmt(v):
    return "n/a" if v is None else f"{v:.1f}"


def fva_text(fva, ids, uptake):
    if not fva or "error" in fva or not all(i in fva for i in ids):
        return ""
    if len(ids) != 1:
        return ""
    lo, hi = fva[ids[0]]
    lo, hi = 100 * lo / uptake, 100 * hi / uptake
    digits = 1 if max(abs(lo), abs(hi)) < 10 else 0
    return f" [{lo:.{digits}f}, {hi:.{digits}f}]"


def panel(ax, title, vector, fva, removed=()):
    f, uptake = per100(vector)
    ax.set_xlim(-0.02, 1.08)
    ax.set_ylim(-0.04, 1.0)
    ax.axis("off")
    ax.axhspan(0.5, 1.0, color="#eef4fb", zorder=0)
    ax.axhspan(-0.04, 0.5, color="#fbf3ea", zorder=0)
    ax.text(0.0, 0.97, "cytosol", fontsize=8, color="#335", va="top")
    ax.text(0.0, 0.47, "mitochondrial matrix", fontsize=8, color="#533", va="top")
    ax.axhline(0.5, color="#888", lw=1.2, ls="-")
    for key, (x, y) in NODES.items():
        ax.add_patch(FancyBboxPatch((x - 0.055, y - 0.025), 0.11, 0.05, boxstyle="round,pad=0.005",
                                    fc="white", ec="#444", lw=0.8, zorder=3))
        ax.text(x, y, LABELS[key], ha="center", va="center", fontsize=7, zorder=4)
    for edge in EDGES:
        label, src, dst, fn, ids = edge[:5]
        rad = edge[5] if len(edge) > 5 else 0.12
        rid_present = all(i in f or i in removed for i in ids)
        is_removed = any(i in removed for i in ids)
        value = None if is_removed else (fn(f) if rid_present else None)
        (x0, y0), (x1, y1) = NODES[src], NODES[dst]
        if value is not None and value < -1e-9:
            (x0, y0), (x1, y1) = (x1, y1), (x0, y0)
        active = value is not None and abs(value) > 1e-6
        color = "#c0392b" if is_removed else ("#1f4e79" if active else "#bbbbbb")
        arrow = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9,
                                lw=0.6 + min(abs(value or 0) / 60, 3.5), color=color,
                                linestyle="--" if is_removed else "-", shrinkA=14, shrinkB=14, zorder=2,
                                connectionstyle=f"arc3,rad={rad}")
        ax.add_patch(arrow)
        text = f"{label.split(' (')[0]}: " + ("REMOVED" if is_removed else fmt(abs(value) if value is not None else None))
        text += "" if is_removed else fva_text(fva, ids, uptake)
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        mx, my = mx + rad * 0.5 * (y1 - y0), my - rad * 0.5 * (x1 - x0)  # follow the arc bulge
        ax.text(mx, my, text, fontsize=5.6, ha="center", color=color, zorder=5,
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.3))
    ax.set_title(title + f"\n(uptake {uptake:.3f} mmol/gDW/h; growth {vector['biomass_C']:.4f} h$^{{-1}}$)", fontsize=8)


def main():
    run = json.loads((HERE / "run.json").read_text())
    route = {key: json.loads((RUNS / f"pfba_route_{key}.json").read_text()) for key in ("baseline", "candidate")}
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    for row, label in enumerate(("WT", "complex_I_KO")):
        for col, key in enumerate(("baseline", "candidate")):
            title = f"{'WT' if label == 'WT' else 'complex I (R2062) KO'} - " + (
                "before (baseline 0ac2addf)" if key == "baseline" else "after (R132 removed, a56f3649)")
            panel(axes[row, col], title, route[key][label], run[key]["route"][label]["fva"],
                  removed=("R132",) if key == "candidate" else ())
    fig.suptitle("Matrix NADH routes around R132 - PO1f + CSM-Leu (SENSITIVITY: medium adds amino acids), glucose <= 2.414. Flux per 100 glucose (one pFBA "
                 "optimum); [min, max] = FVA at the growth optimum.\nArrow points in the flux direction; grey = zero; "
                 "dashed red = removed step. Label value = |pFBA flux| in the arrow direction; brackets = FVA range in the reaction's written direction (signed).", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(HERE / "map_r132_route.png", dpi=160)
    step5 = json.loads((RUNS / "pfba_step5.json").read_text())
    rows = {(r["substrate"], r["variant"]): r for r in run["step5"]}
    variants = ("baseline_as_supplied", "candidate_caps_R361", "candidate_caps_gcv_module")
    fig, axes = plt.subplots(1, 3, figsize=(21, 6.5))
    for ax, name in zip(axes, variants):
        key = f"glucose|{name}"
        if key in step5 and "error" not in step5[key]:
            panel(ax, f"glucose - {name}", step5[key], rows[("glucose", name)].get("fva"),
                  removed=("R132",) if name.startswith("candidate") else ())
    fig.suptitle("Step 5 (DIAGNOSTIC; caps are an ASSUMPTION): glucose, PO1f + CSM-Leu (SENSITIVITY: adds amino acids; no overlay). Flux per 100 glucose (one pFBA "
                 "optimum); [min, max] = FVA at the growth optimum. Label value = |pFBA flux| in the arrow direction; brackets = FVA range in the reaction's written direction (signed).", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(HERE / "map_step5_glucose.png", dpi=150)
    # W29-type defined medium, V2 (primary): uracil closed, R612 restored, amino acids closed, uptake fixed.
    w29 = json.loads((HERE / "w29_defined_run.json").read_text())
    vectors = json.loads((RUNS.parent / "runs_w29_defined" / "pfba_vectors.json").read_text())
    absolute = lambda f: {k: [v[0] * 2.414 / 100, v[1] * 2.414 / 100] for k, v in f.items()}  # stored per 100 glucose
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))
    for row, label in enumerate(("WT", "complex_I_KO")):
        for col, key in enumerate(("baseline", "candidate")):
            title = f"{'WT' if label == 'WT' else 'complex I (R2062) KO'} - " + (
                "before (baseline 0ac2addf)" if key == "baseline" else "after (R132 removed, a56f3649)")
            panel(axes[row, col], title, vectors[f"route|V2|{key}|{label}"],
                  absolute(w29["route"][f"V2|{key}|{label}"]["fva_per_100"]), removed=("R132",) if key == "candidate" else ())
    fig.suptitle("Matrix NADH routes around R132 - W29-type defined medium V2 (PRIMARY: no amino acids, uracil closed, "
                 "URA3 restored, glucose fixed 2.414).\nFlux per 100 glucose (one pFBA optimum); [min, max] = FVA at the "
                 "growth optimum. Arrow points in the flux direction; grey = zero; dashed red = removed step. Label value = |pFBA flux| in the arrow direction; brackets = FVA range in the reaction's written direction (signed).", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(HERE / "map_r132_route_W29V2.png", dpi=160)
    rows = {r["variant"]: r for r in w29["step5"] if r["medium"] == "V2" and r["substrate"] == "glucose"}
    variants = ("as_supplied", "holistic_overlay", "holistic_shuttles_R361", "holistic_shuttles_gcv_module")
    fig, axes = plt.subplots(1, 4, figsize=(28, 6.5))
    for ax, name in zip(axes, variants):
        panel(ax, f"glucose V2 - {name}", vectors[f"step5|V2|glucose|{name}"], absolute(rows[name]["fva_per_100"]),
              removed=("R132",) if name.startswith("holistic_shuttles") else ())
    fig.suptitle("Step 5 (DIAGNOSTIC; overlay and caps are ASSUMPTIONS): glucose, W29-type defined medium V2 (primary). "
                 "Flux per 100 glucose (one pFBA optimum); [min, max] = FVA where computed (R570, R2062, R614, R690). Label value = |pFBA flux| in the arrow direction; brackets = FVA range in the reaction's written direction (signed).",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(HERE / "map_step5_glucose_W29V2.png", dpi=140)
    print("maps written")


if __name__ == "__main__":
    main()
