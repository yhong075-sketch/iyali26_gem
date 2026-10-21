"""Score the two b3b9fe22 screens (no overlay; diagnostic overlay) and compare each, gene by gene,
with the matching a56f3649 screen. No optimization, no figure.

Inputs (hashes checked against their manifests): the b3b9fe22 manifests in this folder and in
../screen_b3b9fe22_w29_20261021/, and the a56f3649 manifests in
../confusion_matrices_csmleu_w29_20261014/ (no overlay) and ../screen_a56f3649_overlay_w29_20261014/.

Usage (from the repository root):
    python -B model/reports/screen_b3b9fe22_overlay_w29_20261021/summarize.py
"""
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPORTS = HERE.parent
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
SCREENS = {  # condition -> (b3b9fe22 manifest, a56f3649 manifest)
    "no_overlay": (REPORTS / "screen_b3b9fe22_w29_20261021/screen_b3b9fe22.manifest.json",
                   REPORTS / "confusion_matrices_csmleu_w29_20261014/screen_a56f3649.manifest.json"),
    "overlay": (HERE / "screen.manifest.json",
                REPORTS / "screen_a56f3649_overlay_w29_20261014/screen.manifest.json"),
}


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


def score(manifest, rows):
    wt = manifest["WT"]["raw_growth"]
    mtdna = frozenset(manifest["mtdna_genes"])
    unscorable = frozenset(r["gene_id"] for r in rows if not r["scorable"])
    views = dict(benchmark_rule=[classify(rows, wt, c, mtdna) for c in CUTOFFS],
                 lethal_view=[classify(rows, wt, c, mtdna, unscorable) for c in CUTOFFS])
    assert views["benchmark_rule"] == manifest["metrics_mtdna_out_of_scope"]
    if "ngam_diagnostic" in manifest:
        diagnosed = {d["knockout"]: d["label"] for d in manifest["ngam_diagnostic"] if d["kind"] == "gene"}
        assert set(diagnosed) == unscorable, "Every gene with no solution must be diagnosed"
        maintenance = frozenset(g for g, label in diagnosed.items() if label == "cannot meet maintenance")
        views["maintenance_view"] = [classify(rows, wt, c, mtdna, maintenance) for c in CUTOFFS]
    return views


def compare(new_manifest, new_rows, old_manifest, old_rows):
    new_wt, old_wt = new_manifest["WT"]["raw_growth"], old_manifest["WT"]["raw_growth"]
    new = {r["gene_id"]: r for r in new_rows}
    old = {r["gene_id"]: r for r in old_rows}
    shared = sorted(set(new) & set(old))
    changes = {f"{int(c * 100)}pct": [
        dict(gene_id=g, experimental_positive=new[g]["experimental_positive"],
             a56f3649=old[g][f"prediction_{int(c * 100)}pct"], b3b9fe22=new[g][f"prediction_{int(c * 100)}pct"],
             a56f3649_ratio=old[g]["ko_wt_ratio"], b3b9fe22_ratio=new[g]["ko_wt_ratio"],
             a56f3649_reactions=old[g]["model_reaction_ids"], b3b9fe22_reactions=new[g]["model_reaction_ids"])
        for g in shared if old[g][f"prediction_{int(c * 100)}pct"] != new[g][f"prediction_{int(c * 100)}pct"]]
        for c in CUTOFFS}
    only_new = sorted(set(new) - set(old))
    return dict(wt_a56f3649=old_wt, wt_b3b9fe22=new_wt, shared_genes=len(shared),
                only_in_b3b9fe22=[dict(gene_id=g, experimental_positive=new[g]["experimental_positive"],
                                       prediction_15pct=new[g]["prediction_15pct"],
                                       model_reaction_ids=new[g]["model_reaction_ids"]) for g in only_new],
                only_in_a56f3649=sorted(set(old) - set(new)), call_changes=changes)


def main():
    summary = dict(model="b3b9fe22cba4c8ae4408afe61740bfe75e661c383dd47a224965f06469618543",
                   scoring_negative_class="user_defined_native_model_genes_absent_from_positive_list",
                   validation_status="development_reference_not_independent_validation", conditions={})
    for condition, (new_path, old_path) in SCREENS.items():
        new_manifest, new_rows = load(new_path)
        old_manifest, old_rows = load(old_path)
        entry = dict(screen_manifest=str(Path(new_path).relative_to(REPORTS.parent)),
                     wt_growth=new_manifest["WT"]["raw_growth"],
                     metrics_mtdna_out_of_scope=score(new_manifest, new_rows),
                     unresolved_genes=new_manifest["unresolved_genes"],
                     comparison_with_a56f3649=compare(new_manifest, new_rows, old_manifest, old_rows))
        if "controls" in new_manifest:
            entry["controls"] = {k: dict(reaction=v["reaction"], status=v["solver_status"], ko_wt_ratio=v["ko_wt_ratio"])
                                 for k, v in new_manifest["controls"].items()}
            entry["ngam_label_counts"] = {label: sum(d["label"] == label for d in new_manifest["ngam_diagnostic"])
                                          for label in sorted({d["label"] for d in new_manifest["ngam_diagnostic"]})}
        entry["calls_used"] = new_manifest["calls_used"]
        summary["conditions"][condition] = entry
        i = CUTOFFS.index(0.15)
        for view, metrics in entry["metrics_mtdna_out_of_scope"].items():
            print(condition, view, {k: metrics[i][k] for k in ("TP", "FN", "FP", "TN", "unresolved", "mcc")})
        cmp_ = entry["comparison_with_a56f3649"]
        print(condition, "vs a56f3649: only new genes", len(cmp_["only_in_b3b9fe22"]), "only old", cmp_["only_in_a56f3649"],
              "call changes", {k: len(v) for k, v in cmp_["call_changes"].items()})
    (HERE / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
