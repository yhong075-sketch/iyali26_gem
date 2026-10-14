"""Score the a56f3649 overlay screen three ways, compare it with Task A's 0ac2addf overlay screen
and draw the two 15 % panels (benchmark rule; no solution counted as lethal).

No optimization. Inputs: screen.manifest.json in this folder, its raw gene_screen.tsv, and Task A's
raw overlay screen and NGAM diagnostic for the parent 0ac2addf (hashes checked).

Usage (from the repository root, with IYALI26_RESEARCH_ROOT set):
    python -B model/reports/screen_a56f3649_overlay_w29_20261014/summarize.py
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
PRIMARY = 0.15
TASK_A_SCREEN = ("artifacts/tasks/w29_glucose_uptake_20261012/gene_screen_holistic_overlay_diagnostic.tsv",
                 "e9f4ca6f9a3dcd74bf142cb80b1c33f9c03eecbe2aa495821611c4230a2e42ca")
# Task A is uncommitted in worktree competent-turing-bc6e9b; its WT is run.json
# holistic_overlay_diagnostic.WT (run.json SHA-256 18f44c8d...).
TASK_A_WT = 0.3117523150525423
TASK_A_NGAM = (Path("/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/"
                    "competent-turing-bc6e9b/model/reports/w29_glucose_uptake_20261012/ngam_diagnostic.tsv"),
               "2d2a245dc0cfa602437bc49a9b9756e7e1d0772c33890e787ad520030bf11ab6")
VIEWS = {"benchmark_rule": "Benchmark rule: no-solution KOs unresolved",
         "lethal_view": "Labelled view: no-solution KOs counted as lethal"}
# Reference sequential blue ramp, steps 100-700 (dataviz skill, references/palette.md).
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SURFACE, TEXT_PRIMARY, TEXT_SECONDARY = "#fcfcfb", "#0b0b0b", "#52514e"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_screen(path):
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for r in rows:
        r["experimental_positive"] = r["experimental_positive"] == "True"
        r["scorable"] = r["scorable"] == "True"
        r["raw_growth"] = float(r["raw_growth"]) if r["raw_growth"] else None
    return rows


def classify(rows, wt, cutoff, excluded=frozenset(), lethal_ids=frozenset()):
    # lethal_ids: unscorable knockouts counted as lethal in a labelled view; others stay unresolved.
    counts = dict(TP=0, FN=0, FP=0, TN=0, unresolved=0)
    for r in rows:
        if r["gene_id"] in excluded:
            continue
        if not r["scorable"] and r["gene_id"] not in lethal_ids:
            counts["unresolved"] += 1
            continue
        essential = True if not r["scorable"] else r["raw_growth"] / wt < cutoff
        positive = r["experimental_positive"]
        counts[("TP" if positive else "FP") if essential else ("FN" if positive else "TN")] += 1
    tp, fn, fp, tn = (counts[k] for k in ("TP", "FN", "FP", "TN"))
    ratio = lambda a, b: a / b if b else None
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return dict(cutoff=cutoff, **counts, excluded_out_of_scope=len(excluded),
                recall=ratio(tp, tp + fn), precision=ratio(tp, tp + fp),
                specificity=ratio(tn, tn + fp),
                mcc=(tp * tn - fp * fn) / denominator if denominator else None)


def task_a_comparison(rows, wt, diagnostics, research):
    path = research / TASK_A_SCREEN[0]
    assert sha(path) == TASK_A_SCREEN[1], "Task A raw overlay screen changed"
    old = {r["gene_id"]: r for r in csv.DictReader(path.open(), delimiter="\t")}
    new = {r["gene_id"]: r for r in rows}
    assert set(old) == set(new)
    status_changes = sorted(g for g in new if new[g]["solver_status"] != old[g]["status"])
    both = [g for g in new if new[g]["raw_growth"] is not None and old[g]["growth"]]
    ratio_changes = [dict(gene_id=g, experimental_positive=new[g]["experimental_positive"],
                          parent_ratio=float(old[g]["growth"]) / TASK_A_WT, new_ratio=new[g]["raw_growth"] / wt,
                          model_reaction_ids=new[g]["model_reaction_ids"])
                     for g in sorted(both)
                     if abs(float(old[g]["growth"]) / TASK_A_WT - new[g]["raw_growth"] / wt) > 1e-6]
    changed = {}
    for c in CUTOFFS:
        before = {g: None if not old[g]["growth"] else float(old[g]["growth"]) / TASK_A_WT < c for g in old}
        after = {g: None if not new[g]["scorable"] else new[g]["raw_growth"] / wt < c for g in new}
        changed[f"{int(c * 100)}pct"] = sorted(g for g in new if before[g] != after[g])
    ngam_path, ngam_sha = TASK_A_NGAM
    assert sha(ngam_path) == ngam_sha, "Task A NGAM diagnostic changed"
    parent = {r["target"]: r for r in csv.DictReader(ngam_path.open(), delimiter="\t")
              if r["variant"] == "holistic_overlay_diagnostic"}
    ours = {d["knockout"].removeprefix("control:"): d for d in diagnostics}
    assert set(parent) == set(ours), (sorted(set(parent) ^ set(ours)))
    ngam_diffs = {k: abs(float(parent[k]["max_maintenance_atp"]) - ours[k]["max_maintenance_atp"]) for k in ours}
    label_changes = sorted(k for k in ours if parent[k]["label"] != ours[k]["label"])
    return dict(reference=str(path), reference_sha256=TASK_A_SCREEN[1],
                wt_task_a=TASK_A_WT, wt_now=wt, wt_abs_diff=abs(wt - TASK_A_WT),
                genes_compared=len(new), solver_status_changes=status_changes,
                genes_with_growth_in_both=len(both),
                max_abs_growth_diff=max(abs(new[g]["raw_growth"] - float(old[g]["growth"])) for g in both),
                ratio_changes_above_tolerance=ratio_changes, ratio_tolerance=1e-6, call_changes=changed,
                ngam_reference=str(ngam_path), ngam_reference_sha256=ngam_sha,
                ngam_knockouts_compared=len(ours), ngam_label_changes=label_changes,
                ngam_max_abs_max_maintenance_diff=max(ngam_diffs.values()))


def text_color(fill):
    # Pick the ink with the higher WCAG contrast against the cell fill.
    def luminance(hex_color):
        channel = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in to_rgb(hex_color)]
        return 0.2126 * channel[0] + 0.7152 * channel[1] + 0.0722 * channel[2]
    fill_l = luminance(fill)
    dark = (fill_l + 0.05) / (luminance(TEXT_PRIMARY) + 0.05)
    light = 1.05 / (fill_l + 0.05)
    return TEXT_PRIMARY if dark >= light else "#ffffff"


def draw(manifest, panels, diagnostics_in_scope, path_stem):
    cmap = LinearSegmentedColormap.from_list("seq_blue", RAMP)
    fig = plt.figure(figsize=(12, 7.6), facecolor=SURFACE)
    fig.suptitle("iYali26 newest candidate (a56f3649): single-gene essentiality vs consensus list (KO/WT < 15 %)",
                 x=0.04, y=0.975, ha="left", fontsize=13, fontweight="bold", color=TEXT_PRIMARY)
    fig.text(0.04, 0.925, "PO1f + CSM-Leu, glucose uptake 2.414 (W29) · in-memory diagnostic overlay, "
             "assumed caps R1384 ≤ 8 and R860 ≤ 1", ha="left", va="top", fontsize=11, color=TEXT_SECONDARY)
    controls = manifest["controls"]
    for i, (view, m) in enumerate(panels.items()):
        ax = fig.add_axes([0.12 + i * 0.48, 0.33, 0.32, 0.47])
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
        ax.set_title(VIEWS[view], fontsize=11.5, fontweight="bold", color=TEXT_PRIMARY, pad=34)
        if view == "benchmark_rule":
            second = f"{m['unresolved']} no-solution KOs left out · {m['excluded_out_of_scope']} mtDNA excluded"
        else:
            second = (f"{diagnostics_in_scope} no-solution KOs counted as lethal; all fail only on "
                      f"maintenance ATP")
        lines = [
            f"Recall {m['recall']:.1%} · Precision {m['precision']:.1%} · "
            f"Specificity {m['specificity']:.1%} · MCC {m['mcc']:.3f}",
            second,
            f"WT {manifest['WT']['raw_growth']:.4f} h⁻¹ · {manifest['gene_count']:,} genes screened · "
            f"complex I KO/WT {controls['complex_I']['ko_wt_ratio']:.3f}",
            f"Listed positives in model {manifest['in_model_positive_count']} / "
            f"{manifest['experimental_positive_count']:,} "
            f"({manifest['out_of_scope_positive_count']:,} out of scope, not FN)",
        ]
        for k, text in enumerate(lines):
            ax.text(-0.25, -0.08 - 0.08 * k, text, transform=ax.transAxes, ha="left", va="top",
                    fontsize=9, color=TEXT_PRIMARY if k == 0 else TEXT_SECONDARY)
    footnote = ("Static FBA, po1f_csm_leu_w29uptake_20261012.csv, PO1f strain profile, objective biomass_C; essential = "
                "unrounded KO/WT strictly below 0.15.\n"
                "Overlay = Task A's diagnostic package in memory (not curated); the two caps are assumptions. "
                "The benchmark rule (left) is primary; the lethal view (right) is a labelled second view.\n"
                "Cell shade = share of the experimental row. "
                "Negatives are model genes not on the consensus list (a user-defined class, not experimental "
                "non-essential calls). The list was used during model development,\n"
                "so this is a development reference, not independent validation. Screen run 2026-10-14.")
    fig.text(0.04, 0.03, footnote, ha="left", va="bottom", fontsize=8, color=TEXT_SECONDARY)
    for suffix in (".png", ".pdf"):
        fig.savefig(path_stem.with_suffix(suffix), dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main():
    research = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
    manifest = json.loads((HERE / "screen.manifest.json").read_text())
    assert manifest["status"] == "complete" and manifest["all_inputs_and_code_unchanged"]
    path = Path(manifest["outputs"]) / "gene_screen.tsv"
    assert sha(path) == manifest["output_sha256"][str(path)], "Raw screen changed"
    rows = read_screen(path)
    wt = manifest["WT"]["raw_growth"]
    mtdna = frozenset(manifest["mtdna_genes"])
    diagnostics = manifest["ngam_diagnostic"]
    unscorable = frozenset(r["gene_id"] for r in rows if not r["scorable"])
    diagnosed = {d["knockout"]: d["label"] for d in diagnostics if d["kind"] == "gene"}
    assert set(diagnosed) == unscorable, "Every gene with no solution must be diagnosed"
    maintenance = frozenset(g for g, label in diagnosed.items() if label == "cannot meet maintenance")
    scored = dict(
        benchmark_rule=[classify(rows, wt, c, mtdna) for c in CUTOFFS],
        benchmark_rule_all_native=[classify(rows, wt, c) for c in CUTOFFS],
        lethal_view=[classify(rows, wt, c, mtdna, unscorable) for c in CUTOFFS],
        maintenance_view=[classify(rows, wt, c, mtdna, maintenance) for c in CUTOFFS])
    # Rescoring from the raw rows must reproduce what the screen recorded.
    assert scored["benchmark_rule"] == manifest["metrics_mtdna_out_of_scope"]
    assert scored["benchmark_rule_all_native"] == manifest["metrics_all_native_genes"]
    comparison = task_a_comparison(rows, wt, diagnostics, research)
    index = CUTOFFS.index(PRIMARY)
    panels = {view: scored[view][index] for view in VIEWS}
    in_scope = len(unscorable - mtdna)
    stem = HERE / "confusion_matrices_15pct"
    draw(manifest, panels, in_scope, stem)
    unresolved = [dict(gene_id=r["gene_id"], experimental_positive=r["experimental_positive"],
                       mtdna_out_of_scope=r["gene_id"] in mtdna, model_reaction_ids=r["model_reaction_ids"],
                       ngam_label=diagnosed[r["gene_id"]])
                  for r in rows if not r["scorable"]]
    summary = dict(
        model=manifest["inputs"]["model"], medium=manifest["inputs"]["medium"],
        overlay=manifest["inputs"]["overlay"], overlay_status="ASSUMPTION (diagnostic), in memory only",
        primary_cutoff=PRIMARY, cutoffs=list(CUTOFFS),
        scoring_negative_class="user_defined_native_model_genes_absent_from_positive_list",
        views=dict(benchmark_rule="primary: no solution = unresolved",
                   lethal_view="labelled: every no-solution knockout counted as lethal",
                   maintenance_view="labelled (Task A proposal, not adopted): only 'cannot meet maintenance' counted as lethal"),
        validation_status="development_reference_not_independent_validation",
        wt_growth=wt, controls={k: dict(reaction=v["reaction"], status=v["solver_status"], ko_wt_ratio=v["ko_wt_ratio"])
                                for k, v in manifest["controls"].items()},
        metrics_mtdna_out_of_scope=scored, unresolved_genes=unresolved,
        ngam_label_counts={label: sum(d["label"] == label for d in diagnostics)
                           for label in sorted({d["label"] for d in diagnostics})},
        max_maintenance_atp_range=[min(d["max_maintenance_atp"] for d in diagnostics),
                                   max(d["max_maintenance_atp"] for d in diagnostics)],
        comparison_with_task_a_0ac2addf_overlay=comparison,
        calls_used=manifest["calls_used"],
        figure={str(stem.with_suffix(s).name): sha(stem.with_suffix(s)) for s in (".png", ".pdf")})
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    for view in ("benchmark_rule", "lethal_view", "maintenance_view"):
        print(view, {k: scored[view][index][k] for k in ("TP", "FN", "FP", "TN", "unresolved", "mcc")})
    print("vs Task A 0ac2addf overlay: WT diff", comparison["wt_abs_diff"], "status changes",
          len(comparison["solver_status_changes"]), "ratio changes", len(comparison["ratio_changes_above_tolerance"]),
          "call changes", {k: len(v) for k, v in comparison["call_changes"].items()},
          "NGAM label changes", len(comparison["ngam_label_changes"]),
          "max |max-maintenance diff|", comparison["ngam_max_abs_max_maintenance_diff"])


if __name__ == "__main__":
    main()
