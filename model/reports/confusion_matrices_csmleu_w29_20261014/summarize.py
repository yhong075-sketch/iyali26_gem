"""Score the two screens, run the no-solve comparison check and draw the 15 % matrices.

No optimization. Inputs are the two run_screen.py manifests in this folder, their raw
gene_screen.tsv files and Task B's raw screen of a56f3649 on the same medium
(reports/r132_mito_adh_20261012/run.json).

Copied from the SD-Leu summary of the same two models (commit 01beb03, branch
claude/confusion-matrices-iyali21-fc6ecb). Changes: the two SD-Leu comparisons are replaced by the
Task B check; a second, labelled view counts "no solution" as lethal (as in Task B); the figure
names the medium and shows the unresolved knockouts.

Usage (from the repository root, with IYALI26_RESEARCH_ROOT set):
    python -B model/reports/confusion_matrices_csmleu_w29_20261014/summarize.py
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
LABELS = ("iyali21", "a56f3649")
TITLES = {"iyali21": "iYali21 (iyli21.xml, 6974b758)",
          "a56f3649": "iYali26 newest candidate (a56f3649)"}
TASK_B_SCREEN = ("artifacts/tasks/r132_mito_adh_20261012/runs/gene_screen_candidate.tsv",
                 "c9185ecdb26fcc99620b7522b091640b1938397ead865c2c46958d4427c30024")
TASK_B_WT = 0.3124267446456666  # reports/r132_mito_adh_20261012/run.json, candidate.WT
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


def classify(rows, wt, cutoff, excluded=frozenset(), lethal_view=False):
    counts = dict(TP=0, FN=0, FP=0, TN=0, unresolved=0)
    for r in rows:
        if r["gene_id"] in excluded:
            continue
        if not r["scorable"] and not lethal_view:
            counts["unresolved"] += 1
            continue
        # Labelled second view only: a knockout with no solution is counted as lethal.
        essential = True if not r["scorable"] else r["raw_growth"] / wt < cutoff
        positive = r["experimental_positive"]
        counts[("TP" if positive else "FP") if essential else ("FN" if positive else "TN")] += 1
    tp, fn, fp, tn = (counts[k] for k in ("TP", "FN", "FP", "TN"))
    ratio = lambda a, b: a / b if b else None
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    result = dict(cutoff=cutoff, **counts, excluded_out_of_scope=len(excluded),
                  recall=ratio(tp, tp + fn), precision=ratio(tp, tp + fp),
                  specificity=ratio(tn, tn + fp),
                  mcc=(tp * tn - fp * fn) / denominator if denominator else None)
    if lethal_view:
        result["view"] = "no solution counted as lethal (labelled; not the benchmark rule)"
    return result


def load(label):
    manifest = json.loads((HERE / f"screen_{label}.manifest.json").read_text())
    assert manifest["status"] == "complete" and manifest["all_inputs_and_code_unchanged"]
    path = Path(manifest["outputs"]) / "gene_screen.tsv"
    assert sha(path) == manifest["output_sha256"][str(path)], "Raw screen changed"
    rows = read_screen(path)
    wt = manifest["WT"]["raw_growth"]
    mtdna = frozenset(manifest["mtdna_genes"])
    scored = dict(
        mtdna_out_of_scope=[classify(rows, wt, c, mtdna) for c in CUTOFFS],
        all_native_genes=[classify(rows, wt, c) for c in CUTOFFS],
        lethal_view_mtdna_out_of_scope=[classify(rows, wt, c, mtdna, lethal_view=True) for c in CUTOFFS])
    # Rescoring from the raw rows must reproduce what the screen recorded.
    assert scored["mtdna_out_of_scope"] == manifest["metrics_mtdna_out_of_scope"]
    assert scored["all_native_genes"] == manifest["metrics_all_native_genes"]
    return manifest, rows, wt, scored


def task_b_reproduction(rows, wt, research):
    path = research / TASK_B_SCREEN[0]
    assert sha(path) == TASK_B_SCREEN[1], "Task B raw screen changed"
    old = {r["gene_id"]: r for r in csv.DictReader(path.open(), delimiter="\t")}
    new = {r["gene_id"]: r for r in rows}
    assert set(old) == set(new)
    status_changes = sorted(g for g in new if new[g]["solver_status"] != old[g]["status"])
    both = [g for g in new if new[g]["raw_growth"] is not None and old[g]["growth"]]
    diffs = [abs(new[g]["raw_growth"] - float(old[g]["growth"])) for g in both]
    changed = {}
    for c in CUTOFFS:
        before = {g: None if not old[g]["growth"] else float(old[g]["growth"]) / TASK_B_WT < c for g in old}
        after = {g: None if not new[g]["scorable"] else new[g]["raw_growth"] / wt < c for g in new}
        changed[f"{int(c * 100)}pct"] = sorted(g for g in new if before[g] != after[g])
    return dict(reference=str(path), reference_sha256=TASK_B_SCREEN[1],
                wt_task_b=TASK_B_WT, wt_now=wt, wt_abs_diff=abs(wt - TASK_B_WT),
                genes_compared=len(new), solver_status_changes=status_changes,
                genes_with_growth_in_both=len(both), max_abs_growth_diff=max(diffs),
                call_changes=changed)


def unresolved_genes(rows):
    return [dict(gene_id=r["gene_id"], solver_status=r["solver_status"],
                 experimental_positive=r["experimental_positive"],
                 mtdna_out_of_scope=r["mtdna_out_of_scope"] == "True",
                 model_reaction_ids=r["model_reaction_ids"])
            for r in rows if not r["scorable"]]


def text_color(fill):
    # Pick the ink with the higher WCAG contrast against the cell fill.
    def luminance(hex_color):
        channel = [x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in to_rgb(hex_color)]
        return 0.2126 * channel[0] + 0.7152 * channel[1] + 0.0722 * channel[2]
    fill_l = luminance(fill)
    dark = (fill_l + 0.05) / (luminance(TEXT_PRIMARY) + 0.05)
    light = 1.05 / (fill_l + 0.05)
    return TEXT_PRIMARY if dark >= light else "#ffffff"


def draw(results, path_stem):
    cmap = LinearSegmentedColormap.from_list("seq_blue", RAMP)
    fig = plt.figure(figsize=(12, 7.6), facecolor=SURFACE)
    fig.suptitle("Single-gene essentiality vs consensus essential list (KO/WT < 15 %)",
                 x=0.04, y=0.975, ha="left", fontsize=14, fontweight="bold", color=TEXT_PRIMARY)
    fig.text(0.04, 0.925, "Task A medium + PO1f, glucose uptake 2.414 mmol/gDW/h",
             ha="left", va="top", fontsize=11, color=TEXT_SECONDARY)
    for i, label in enumerate(LABELS):
        manifest, m = results[label]["manifest"], results[label]["primary"]
        lethal = results[label]["primary_lethal_view"]
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
        ax.set_title(TITLES[label], fontsize=12, fontweight="bold", color=TEXT_PRIMARY, pad=34)
        lines = [
            f"Recall {m['recall']:.1%} · Precision {m['precision']:.1%} · "
            f"Specificity {m['specificity']:.1%} · MCC {m['mcc']:.3f}",
            f"WT {manifest['WT']['raw_growth']:.4f} h⁻¹ · {manifest['gene_count']:,} genes screened · "
            f"{m['excluded_out_of_scope']} mtDNA excluded · {m['unresolved']} unresolved",
            f"If no-solution KOs = lethal: TP {lethal['TP']} · FN {lethal['FN']} · "
            f"FP {lethal['FP']} · TN {lethal['TN']} · MCC {lethal['mcc']:.3f}",
            f"Listed positives in model {manifest['in_model_positive_count']} / "
            f"{manifest['experimental_positive_count']:,} "
            f"({manifest['out_of_scope_positive_count']:,} out of scope, not FN)",
        ]
        for k, text in enumerate(lines):
            ax.text(-0.25, -0.08 - 0.08 * k, text, transform=ax.transAxes, ha="left", va="top",
                    fontsize=9, color=TEXT_PRIMARY if k == 0 else TEXT_SECONDARY)
    footnote = ("Static FBA, Task A medium (po1f_csm_leu_w29uptake_20261012.csv, glucose uptake 2.414) + "
                "PO1f strain profile, objective biomass_C; essential = unrounded KO/WT strictly below 0.15.\n"
                "Matrix cells follow the benchmark rule: a knockout with no LP solution is unresolved and left out. "
                "Cell shade = share of the experimental row.\n"
                "Negatives are model genes not on the consensus list (a user-defined class, not experimental "
                "non-essential calls). The list was used during model development,\n"
                "so this is a development reference, not independent validation. Screens run 2026-10-14.")
    fig.text(0.04, 0.03, footnote, ha="left", va="bottom", fontsize=8, color=TEXT_SECONDARY)
    for suffix in (".png", ".pdf"):
        fig.savefig(path_stem.with_suffix(suffix), dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main():
    research = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
    results = {}
    for label in LABELS:
        manifest, rows, wt, scored = load(label)
        index = CUTOFFS.index(PRIMARY)
        results[label] = dict(manifest=manifest, rows=rows, wt=wt, scored=scored,
                              primary=scored["mtdna_out_of_scope"][index],
                              primary_lethal_view=scored["lethal_view_mtdna_out_of_scope"][index])
    reproduction = task_b_reproduction(results["a56f3649"]["rows"], results["a56f3649"]["wt"], research)
    genes = {label: {r["gene_id"] for r in results[label]["rows"]} for label in LABELS}
    # Shared nuclear genes whose 15 % call differs between the two models.
    by_gene = {label: {r["gene_id"]: r for r in results[label]["rows"]} for label in LABELS}
    mtdna = set(results["a56f3649"]["manifest"]["mtdna_genes"])
    differences = []
    for g in sorted((genes["iyali21"] & genes["a56f3649"]) - mtdna):
        old, new = by_gene["iyali21"][g], by_gene["a56f3649"][g]
        if old["prediction_15pct"] != new["prediction_15pct"]:
            differences.append(dict(gene_id=g, experimental_positive=new["experimental_positive"],
                iyali21_call=old["prediction_15pct"], iyali21_ratio=old["ko_wt_ratio"],
                a56f3649_call=new["prediction_15pct"], a56f3649_ratio=new["ko_wt_ratio"],
                iyali21_reactions=old["model_reaction_ids"], a56f3649_reactions=new["model_reaction_ids"]))
    with (HERE / "shared_gene_call_differences_15pct.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(differences[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(differences)
    stem = HERE / "confusion_matrices_15pct"
    draw(results, stem)
    summary = dict(
        primary_cutoff=PRIMARY, cutoffs=list(CUTOFFS),
        medium="model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv",
        scoring_negative_class="user_defined_native_model_genes_absent_from_positive_list",
        unresolved_rule="benchmark rule: no solution = unresolved, left out of the matrix; "
                        "lethal_view_* is a labelled second view only",
        validation_status="development_reference_not_independent_validation",
        models={label: dict(
            model=results[label]["manifest"]["inputs"]["model"],
            screen_manifest=f"screen_{label}.manifest.json",
            wt_growth=results[label]["wt"], gene_count=results[label]["manifest"]["gene_count"],
            in_model_positive_count=results[label]["manifest"]["in_model_positive_count"],
            out_of_scope_positive_count=results[label]["manifest"]["out_of_scope_positive_count"],
            metrics_mtdna_out_of_scope=results[label]["scored"]["mtdna_out_of_scope"],
            metrics_all_native_genes=results[label]["scored"]["all_native_genes"],
            metrics_lethal_view_mtdna_out_of_scope=results[label]["scored"]["lethal_view_mtdna_out_of_scope"],
            unresolved_genes=unresolved_genes(results[label]["rows"]))
            for label in LABELS},
        shared_gene_count=len(genes["iyali21"] & genes["a56f3649"]),
        only_in_iyali21=sorted(genes["iyali21"] - genes["a56f3649"]),
        only_in_a56f3649=sorted(genes["a56f3649"] - genes["iyali21"]),
        shared_nuclear_call_differences_15pct=len(differences),
        a56f3649_reproduction_vs_task_b=reproduction,
        figure={str(stem.with_suffix(s).name): sha(stem.with_suffix(s)) for s in (".png", ".pdf")})
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    for label in LABELS:
        print(label, {k: results[label]["primary"][k] for k in ("TP", "FN", "FP", "TN", "unresolved")},
              "lethal view", {k: results[label]["primary_lethal_view"][k] for k in ("TP", "FN", "FP", "TN")})
    print("a56f3649 vs Task B: WT diff", reproduction["wt_abs_diff"], "max |growth diff|",
          reproduction["max_abs_growth_diff"], "status changes", len(reproduction["solver_status_changes"]),
          "call changes:", {k: len(v) for k, v in reproduction["call_changes"].items()})
    print("shared nuclear call differences at 15 %:", len(differences))


if __name__ == "__main__":
    main()
