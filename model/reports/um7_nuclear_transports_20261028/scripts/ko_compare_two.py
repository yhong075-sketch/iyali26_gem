"""[Logic reused unchanged for UM7 2026-10-28: BASELINE = 9d655c2b, CANDIDATE = 5f54eb22; the text below describes the run it
was first written for.]
Checked single-gene knockouts on the baseline 2c9cebd0 and the new candidate, compared and scored (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/lipoylation_acylacp_on_2c9cebd0_20261023/scripts/ko_compare_two.py \
    BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.tsv OUT_SUMMARY.json
Two-model copy of ../../combined_tpp_heme_r2193_20261022/scripts/ko_compare_three.py (same checked_solve.py:
Gurobi at 1e-9, a model reset before each knockout, the targeted check); only the model set and comparison differ. Flagged points (an optimal point failing the targeted check, or a status other than
optimal/infeasible) get a fresh re-solve, at most FRESH_LIMIT in total; infeasible knockouts stay unresolved.
Class: essential when KO/WT is strictly below the cutoff. Negatives are the user-defined class; mtDNA genes
are out of scope; the list was used in development, so scores are a development reference.
"""
import csv
import json
import math
import sys
import time
import warnings
from pathlib import Path

from checked_solve import (GENERAL_TOLERANCE, KNOCKED_TOLERANCE, LP_TIMEOUT_SECONDS, SOLVER_TOLERANCE,
                           configure, matrix, solve)
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL
from tools.dfba_new_fn_essentiality import load_experimental_reference

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
CHANGE_TOLERANCE = 1e-6
FRESH_LIMIT = 20
MAGNITUDES = ("max_knocked_flux", "max_bound_violation", "max_residual")

baseline_path, candidate_path, medium, profile, experimental, out_tsv, out_json = sys.argv[1:8]
paths = {"baseline": baseline_path, "candidate": candidate_path}
positives = set(load_experimental_reference(Path(experimental))[0])
mtdna = set(json.loads((MODEL.root / "expected/mtdna_genes_out_of_scope.json").read_text())["genes"])
fresh_used = []


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


def closed(model, gene):
    """Reactions of the gene that the knockout leaves closed (call inside the knockout context)."""
    return {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}


def fresh(path, gene):
    ctx = context(path)
    with ctx.model:
        ctx.model.genes.get_by_id(gene).knock_out()
        return solve(ctx.model, matrix(ctx.model), closed(ctx.model, gene))


def knockouts(label, path):
    ctx = context(path)
    model = ctx.model
    stoichiometry = matrix(model)
    wt = solve(model, stoichiometry, ())
    rows = {}
    for index, gene in enumerate(sorted(g.id for g in model.genes if g.id not in ctx.excluded_runtime_genes), 1):
        started = time.monotonic()
        with model:
            model.genes.get_by_id(gene).knock_out()
            row = solve(model, stoichiometry, closed(model, gene))
        row["fresh_resolve"] = None
        if not row["passed"] and row["status"] != "infeasible":
            if len(fresh_used) >= FRESH_LIMIT:
                raise RuntimeError(f"More than {FRESH_LIMIT} knockouts need a fresh re-solve; stopping as declared")
            fresh_used.append([label, gene])
            row["fresh_resolve"] = fresh(path, gene)
        rows[gene] = row
        print(f"{ctx.canonical_model_sha256[:8]} {index} {gene} {row['status']} passed={row['passed']}"
              f" fresh={None if row['fresh_resolve'] is None else row['fresh_resolve']['passed']}"
              f" {time.monotonic() - started:.2f}s", file=sys.stderr, flush=True)
    return ctx, wt, rows


def final(row):
    for point in (row, row.get("fresh_resolve") or {}):
        if point.get("passed"):
            return point["growth"]
    return None


def call(wt, growth, cutoff):
    if growth is None or not wt:
        return "unresolved"
    return "essential" if growth / wt < cutoff else "non-essential"


def score(wt, rows):
    result = {}
    for cutoff in CUTOFFS:
        counts = {"TP": 0, "FN": 0, "FP": 0, "TN": 0, "unresolved": 0}
        for gene, row in rows.items():
            if gene in mtdna:
                continue
            c = call(wt, final(row), cutoff)
            key = "unresolved" if c == "unresolved" else (
                ("TP" if gene in positives else "FP") if c == "essential" else ("FN" if gene in positives else "TN"))
            counts[key] += 1
        tp, fn, fp, tn = (counts[k] for k in ("TP", "FN", "FP", "TN"))
        d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        result[str(cutoff)] = {**counts, "recall": tp / (tp + fn) if tp + fn else None,
                               "specificity": tn / (tn + fp) if tn + fp else None,
                               "precision": tp / (tp + fp) if tp + fp else None,
                               "mcc": (tp * tn - fp * fn) / d if d else None}
    return result


def largest(rows, key):
    values = [r[key] for r in rows.values() if r.get(key) is not None]
    return {"max": max(values), "count_over_1e-9": sum(v > 1e-9 for v in values),
            "count_over_1e-7": sum(v > 1e-7 for v in values), "count_over_1e-6": sum(v > 1e-6 for v in values)}


runs = {label: knockouts(label, path) for label, path in paths.items()}
assert all(wt["passed"] for _, wt, _ in runs.values()), "A WT point failed the check"
wts = {label: wt["growth"] for label, (_, wt, _) in runs.items()}
rows = {label: r for label, (_, _, r) in runs.items()}
genes = sorted(set().union(*rows.values()))


def compare(a, b):
    changed_growth, changed_class = [], {str(c): [] for c in CUTOFFS}
    for gene in genes:
        x, y = rows[a].get(gene), rows[b].get(gene)
        gx, gy = (final(x) if x else None), (final(y) if y else None)
        cx = [call(wts[a], gx, c) if x else "absent" for c in CUTOFFS]
        cy = [call(wts[b], gy, c) if y else "absent" for c in CUTOFFS]
        if (x is None) != (y is None) or (gx is None) != (gy is None) or (
                gx is not None and gy is not None and abs(gx - gy) > CHANGE_TOLERANCE):
            changed_growth.append({"gene": gene, a: gx if x else "absent", b: gy if y else "absent"})
        for cutoff, p, q in zip(CUTOFFS, cx, cy):
            if p != q:
                changed_class[str(cutoff)].append({"gene": gene, "experimental_positive": gene in positives, a: p, b: q})
    return {"growth_changed_over_1e-6": changed_growth, "class_changed": changed_class}


with open(out_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    header = ["gene", "experimental_positive", "mtdna_out_of_scope"]
    for label in paths:
        header += [f"{label}_status", f"{label}_passed", f"{label}_fresh_passed", f"{label}_ko"]
        header += [f"{label}_{m}" for m in MAGNITUDES] + [f"{label}_class_{c}" for c in CUTOFFS]
    writer.writerow(header)
    for gene in genes:
        line = [gene, gene in positives, gene in mtdna]
        for label in paths:
            r = rows[label].get(gene)
            g = final(r) if r else None
            line += [r["status"] if r else "absent", r["passed"] if r else None,
                     (r["fresh_resolve"] or {}).get("passed") if r else None, g]
            line += [r.get(m) if r else None for m in MAGNITUDES]
            line += [call(wts[label], g, c) if r else "absent" for c in CUTOFFS]
        writer.writerow(line)
summary = {
    "condition": f"Task A medium + PO1f, Gurobi FeasibilityTol/OptimalityTol {SOLVER_TOLERANCE}, reset per knockout, targeted check",
    "check": {"knocked_reaction_flux_tolerance": KNOCKED_TOLERANCE, "general_tolerance": GENERAL_TOLERANCE},
    "lp_timeout_seconds": LP_TIMEOUT_SECONDS,
    "models": {label: {"sha256": ctx.canonical_model_sha256, "genes": len(rows[label])} for label, (ctx, _, _) in runs.items()},
    "medium_sha256": runs["candidate"][0].medium_sha256, "strain_profile_sha256": runs["candidate"][0].strain_profile_sha256,
    "negative_class": "user-defined: model genes not on the consensus list",
    "validation_status": "development reference, not independent validation",
    "wt": {label: wt for label, (_, wt, _) in runs.items()},
    "in_model_positives": {label: len(positives & set(rows[label])) for label in paths},
    "magnitudes_in_sweep": {label: {key: largest(rows[label], key) for key in MAGNITUDES} for label in paths},
    "sweep_points_failing_check": {label: sorted(g for g, r in rows[label].items() if r["status"] == "optimal" and not r["passed"]) for label in paths},
    "non_optimal_in_sweep": {label: sorted(g for g, r in rows[label].items() if r["status"] != "optimal") for label in paths},
    "fresh_resolves": fresh_used,
    "unresolved_after_fresh": {label: sorted(g for g, r in rows[label].items() if final(r) is None) for label in paths},
    "scores": {label: score(wts[label], rows[label]) for label in paths},
    "candidate_vs_baseline": compare("baseline", "candidate"),
    "lp_count": {"sweep": sum(len(r) + 1 for r in rows.values()), "fresh": len(fresh_used)},
}
with open(out_json, "w") as handle:
    json.dump(summary, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: summary[k] for k in ("wt", "in_model_positives", "sweep_points_failing_check", "fresh_resolves", "lp_count")}, indent=1))
