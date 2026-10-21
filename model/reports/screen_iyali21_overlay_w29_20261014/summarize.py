"""Score the iYali21 overlay screen three ways and compare it with iYali21 without the overlay.

No optimization, no figure (the paired figure is drawn in figure_iyali21_vs_a56f3649_20261014).
Inputs: screen.manifest.json here, its raw gene_screen.tsv, and the no-overlay iYali21 raw screen
from confusion_matrices_csmleu_w29_20261014 (hashes checked).

Usage (from the repository root, with IYALI26_RESEARCH_ROOT set):
    python -B model/reports/screen_iyali21_overlay_w29_20261014/summarize.py
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
NO_OVERLAY = ("artifacts/tasks/confusion_matrices_csmleu_w29_20261014/screen_iyali21/gene_screen.tsv",
              "216ff295d3e1e2b863a68a4f9215919504dd7e2427601956179edb7593be9af8")


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
    assert scored["benchmark_rule"] == manifest["metrics_mtdna_out_of_scope"]
    assert scored["benchmark_rule_all_native"] == manifest["metrics_all_native_genes"]
    # Compare with iYali21 without the overlay (same medium and profile).
    base_path = research / NO_OVERLAY[0]
    assert sha(base_path) == NO_OVERLAY[1], "No-overlay screen changed"
    base = {r["gene_id"]: r for r in read_screen(base_path)}
    new = {r["gene_id"]: r for r in rows}
    assert set(base) == set(new)
    call_changes = {f"{int(c * 100)}pct": [
        dict(gene_id=g, experimental_positive=new[g]["experimental_positive"], mtdna=g in mtdna,
             no_overlay=base[g][f"prediction_{int(c * 100)}pct"], overlay=new[g][f"prediction_{int(c * 100)}pct"],
             no_overlay_ratio=base[g]["ko_wt_ratio"], overlay_ratio=new[g]["ko_wt_ratio"],
             model_reaction_ids=new[g]["model_reaction_ids"])
        for g in sorted(new) if base[g][f"prediction_{int(c * 100)}pct"] != new[g][f"prediction_{int(c * 100)}pct"]]
        for c in CUTOFFS}
    both = [g for g in new if new[g]["raw_growth"] is not None and base[g]["raw_growth"] is not None]
    ratio_moves = sorted(g for g in both if abs(float(new[g]["ko_wt_ratio"]) - float(base[g]["ko_wt_ratio"])) > 1e-6)
    summary = dict(
        model=manifest["inputs"]["model"], medium=manifest["inputs"]["medium"],
        overlay=manifest["inputs"]["overlay"], overlay_status="ASSUMPTION (diagnostic), in memory only; adapted to iYali21",
        overlay_record=manifest["diagnostic_overlay"]["record"],
        scoring_negative_class="user_defined_native_model_genes_absent_from_positive_list",
        validation_status="development_reference_not_independent_validation",
        wt_growth=wt, controls={k: dict(reaction=v["reaction"], status=v["solver_status"], ko_wt_ratio=v["ko_wt_ratio"])
                                for k, v in manifest["controls"].items()},
        metrics_mtdna_out_of_scope=scored,
        unresolved_genes=[dict(gene_id=r["gene_id"], experimental_positive=r["experimental_positive"],
                               mtdna_out_of_scope=r["gene_id"] in mtdna, model_reaction_ids=r["model_reaction_ids"],
                               ngam_label=diagnosed[r["gene_id"]]) for r in rows if not r["scorable"]],
        ngam_label_counts={label: sum(d["label"] == label for d in diagnostics)
                           for label in sorted({d["label"] for d in diagnostics})},
        max_maintenance_atp={d["knockout"]: d["max_maintenance_atp"] for d in diagnostics if d["kind"] == "reaction"},
        comparison_with_no_overlay=dict(reference=str(base_path), reference_sha256=NO_OVERLAY[1],
                                        wt_no_overlay=float(json.loads((HERE.parent / "confusion_matrices_csmleu_w29_20261014/screen_iyali21.manifest.json").read_text())["WT"]["raw_growth"]),
                                        call_changes=call_changes, ratio_moves_above_1e6=ratio_moves),
        calls_used=manifest["calls_used"])
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    i = CUTOFFS.index(0.15)
    for view in ("benchmark_rule", "lethal_view", "maintenance_view"):
        print(view, {k: scored[view][i][k] for k in ("TP", "FN", "FP", "TN", "unresolved", "mcc")})
    print("vs no overlay: call changes", {k: len(v) for k, v in call_changes.items()}, "ratio moves", len(ratio_moves))


if __name__ == "__main__":
    main()
