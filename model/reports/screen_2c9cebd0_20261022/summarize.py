"""Score the 2c9cebd0 screen and compare it, gene by gene, with the build's own sweep and with the
b3b9fe22 screen (same medium, no overlay). No optimization, no figure.

Inputs (hashes checked): screen_2c9cebd0.manifest.json here and its raw gene_screen.tsv; the build's
sweep table (`git show origin/main:model/reports/combined_tpp_heme_r2193_20261022/ko_compare.tsv`,
exported to the research root); the b3b9fe22 manifest in ../screen_b3b9fe22_w29_20261021/.

Usage (from the repository root, with IYALI26_RESEARCH_ROOT set):
    python -B model/reports/screen_2c9cebd0_20261022/summarize.py
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
BUILD_SWEEP = ("artifacts/tasks/screen_2c9cebd0_20261022/input/ko_compare_build_sweep.tsv",
               "4aee38ac4e10a59672863e40fa338093649ec698a674fb5de8db8e16d7f2dadb")
B3B9 = HERE.parent / "screen_b3b9fe22_w29_20261021/screen_b3b9fe22.manifest.json"
CALL = {"essential": "essential", "non-essential": "nonessential", "nonessential": "nonessential",
        "unresolved": "unresolved"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(manifest_path):
    manifest = json.loads(Path(manifest_path).read_text())
    assert manifest["status"] == "complete" and manifest["all_inputs_and_code_unchanged"]
    path = Path(manifest["outputs"]) / "gene_screen.tsv"
    assert sha(path) == manifest["output_sha256"][str(path)], f"Raw screen changed: {path}"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for r in rows:
        r["experimental_positive"] = r["experimental_positive"] == "True"
        r["scorable"] = r["scorable"] == "True"
        r["raw_growth"] = float(r["raw_growth"]) if r["raw_growth"] else None
    return manifest, rows


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
    manifest, rows = load(HERE / "screen_2c9cebd0.manifest.json")
    wt = manifest["WT"]["raw_growth"]
    mtdna = frozenset(manifest["mtdna_genes"])
    unscorable = frozenset(r["gene_id"] for r in rows if not r["scorable"])
    scored = dict(benchmark_rule=[classify(rows, wt, c, mtdna) for c in CUTOFFS],
                  benchmark_rule_all_native=[classify(rows, wt, c) for c in CUTOFFS],
                  lethal_view=[classify(rows, wt, c, mtdna, unscorable) for c in CUTOFFS])
    assert scored["benchmark_rule"] == manifest["metrics_mtdna_out_of_scope"]
    assert scored["benchmark_rule_all_native"] == manifest["metrics_all_native_genes"]
    new = {r["gene_id"]: r for r in rows}
    optimal = [r for r in rows if r["solver_status"] == "optimal"]
    numerics = dict(optimal_rows=len(optimal),
                    max_mass_residual=max(float(r["max_mass_residual"]) for r in optimal),
                    max_bound_violation=max(float(r["max_bound_violation"]) for r in optimal),
                    max_knocked_abs_flux=max(float(r["knocked_max_abs_flux"]) for r in optimal))
    # 1. The build's own sweep (combined_* columns), same model and medium.
    path = research / BUILD_SWEEP[0]
    assert sha(path) == BUILD_SWEEP[1], "Build sweep table changed"
    build = {r["gene"]: r for r in csv.DictReader(path.open(), delimiter="\t")}
    assert set(build) == set(new)
    both = [g for g in new if new[g]["raw_growth"] is not None and build[g]["combined_ko"]]
    reproduction = dict(
        reference=str(path), reference_sha256=BUILD_SWEEP[1],
        status_changes=sorted(g for g in new if new[g]["solver_status"] != build[g]["combined_status"]),
        genes_with_growth_in_both=len(both),
        max_abs_growth_diff=max(abs(new[g]["raw_growth"] - float(build[g]["combined_ko"])) for g in both),
        class_changes={f"{int(c * 100)}pct": sorted(
            g for g in new if new[g][f"prediction_{int(c * 100)}pct"] != CALL[build[g][f"combined_class_{c}"]])
            for c in CUTOFFS})
    # 2. b3b9fe22 (no overlay, same medium; Gurobi at 1e-7 there).
    old_manifest, old_rows = load(B3B9)
    old = {r["gene_id"]: r for r in old_rows}
    shared = sorted(set(new) & set(old))
    vs_b3b9 = dict(
        wt_b3b9fe22=old_manifest["WT"]["raw_growth"], wt_2c9cebd0=wt, shared_genes=len(shared),
        only_in_2c9cebd0=[dict(gene_id=g, experimental_positive=new[g]["experimental_positive"],
                               prediction_15pct=new[g]["prediction_15pct"], ko_wt_ratio=new[g]["ko_wt_ratio"],
                               model_reaction_ids=new[g]["model_reaction_ids"]) for g in sorted(set(new) - set(old))],
        only_in_b3b9fe22=sorted(set(old) - set(new)),
        call_changes={f"{int(c * 100)}pct": [
            dict(gene_id=g, experimental_positive=new[g]["experimental_positive"],
                 b3b9fe22=old[g][f"prediction_{int(c * 100)}pct"], c2c9cebd0=new[g][f"prediction_{int(c * 100)}pct"],
                 b3b9fe22_ratio=old[g]["ko_wt_ratio"], c2c9cebd0_ratio=new[g]["ko_wt_ratio"],
                 model_reaction_ids=new[g]["model_reaction_ids"])
            for g in shared if old[g][f"prediction_{int(c * 100)}pct"] != new[g][f"prediction_{int(c * 100)}pct"]]
            for c in CUTOFFS})
    summary = dict(model=manifest["inputs"]["model"], medium=manifest["inputs"]["medium"],
                   solver=manifest["solver"], wt_growth=wt,
                   scoring_negative_class="user_defined_native_model_genes_absent_from_positive_list",
                   validation_status="development_reference_not_independent_validation",
                   metrics_mtdna_out_of_scope=scored, numerics=numerics,
                   unresolved_genes=manifest["unresolved_genes"],
                   reproduction_of_build_sweep=reproduction, comparison_with_b3b9fe22=vs_b3b9,
                   calls_used=manifest["calls_used"])
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    i = CUTOFFS.index(0.15)
    for view in ("benchmark_rule", "lethal_view"):
        print(view, {k: scored[view][i][k] for k in ("TP", "FN", "FP", "TN", "unresolved", "mcc")})
    print("numerics", numerics)
    print("vs build sweep: status changes", len(reproduction["status_changes"]), "max |growth diff|",
          reproduction["max_abs_growth_diff"], "class changes", {k: len(v) for k, v in reproduction["class_changes"].items()})
    print("vs b3b9fe22: only new", [g["gene_id"] for g in vs_b3b9["only_in_2c9cebd0"]], "only old", vs_b3b9["only_in_b3b9fe22"],
          "call changes", {k: len(v) for k, v in vs_b3b9["call_changes"].items()})


if __name__ == "__main__":
    main()
