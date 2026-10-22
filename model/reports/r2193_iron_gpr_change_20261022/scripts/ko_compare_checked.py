"""Checked single-gene knockouts on two models, compared gene by gene and scored (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/r2193_iron_gpr_change_20261022/scripts python \
    ../model/reports/r2193_iron_gpr_change_20261022/scripts/ko_compare_checked.py \
    BASELINE.xml NEW.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.tsv OUT_SUMMARY.json
Adapted from ../heme_iron_biomass_20261022/scripts/ko_compare.py. Changes:
- every knockout starts from a reset GLPK basis, and every point gets the targeted check of
  checked_solve.py: it fails if a reaction the knockout closed carries |v| > 1e-9, or if any bound
  violation or |S v| exceeds 1e-6; all magnitudes are recorded;
- an optimal point that fails the check, or any status other than optimal or infeasible (for example
  time_limit), is solved again in a fresh context (at most FRESH_LIMIT in total); if that also fails,
  the gene is unresolved, never essential. Infeasible knockouts stay unresolved without a re-solve;
- the gene sets may differ: a gene in one model only is reported as absent from the other.
Task A medium + PO1f, GLPK, 60 s per LP. Class: essential when KO/WT is strictly below the cutoff.
Negatives are the user-defined class; mitochondrially encoded genes are out of scope; the list was
used in development, so scores are a development reference, not independent validation.
"""
import csv
import json
import math
import sys
import time
import warnings
from pathlib import Path

from checked_solve import GENERAL_TOLERANCE, KNOCKED_TOLERANCE, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL
from tools.dfba_new_fn_essentiality import load_experimental_reference

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
CHANGE_TOLERANCE = 1e-6
LP_TIMEOUT_SECONDS = 60
FRESH_LIMIT = 20

baseline_path, new_path, medium, profile, experimental, out_tsv, out_json = sys.argv[1:8]
positives = set(load_experimental_reference(Path(experimental))[0])
mtdna = set(json.loads((MODEL.root / "expected/mtdna_genes_out_of_scope.json").read_text())["genes"])
fresh_used = []


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    ctx.model.solver = "glpk"
    ctx.model.solver.configuration.timeout = LP_TIMEOUT_SECONDS
    return ctx


def closed(model, gene):
    """Reactions of the gene that the knockout leaves closed (call inside the knockout context)."""
    return {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}


def fresh(path, gene):
    ctx = context(path)
    with ctx.model:
        ctx.model.genes.get_by_id(gene).knock_out()
        return solve(ctx.model, matrix(ctx.model), closed(ctx.model, gene))


def knockouts(path):
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
        # An infeasible verdict is the solver's answer, not a stale point; it stays unresolved.
        if not row["passed"] and row["status"] != "infeasible":
            if len(fresh_used) >= FRESH_LIMIT:
                raise RuntimeError(f"More than {FRESH_LIMIT} knockouts need a fresh re-solve; stopping as declared")
            fresh_used.append([ctx.canonical_model_sha256[:8], gene])
            row["fresh_resolve"] = fresh(path, gene)
        rows[gene] = row
        print(f"{ctx.canonical_model_sha256[:8]} {index} {gene} {row['status']} passed={row['passed']}"
              f" fresh={None if row['fresh_resolve'] is None else row['fresh_resolve']['passed']}"
              f" {time.monotonic() - started:.2f}s", file=sys.stderr, flush=True)
    return ctx, wt, rows


def final(row):
    """The value used for classification: the checked sweep point, else the checked fresh point."""
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


base_ctx, base_wt, base = knockouts(baseline_path)
new_ctx, new_wt, new = knockouts(new_path)
assert base_wt["passed"] and new_wt["passed"], "WT point failed the check"
bw, nw = base_wt["growth"], new_wt["growth"]
genes = sorted(set(base) | set(new))
changed_growth, changed_class = [], {str(c): [] for c in CUTOFFS}
with open(out_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    magnitudes = ("max_knocked_flux", "max_bound_violation", "max_residual")
    writer.writerow(["gene", "experimental_positive", "mtdna_out_of_scope",
                     "base_status", "base_passed", "base_fresh_passed", "base_ko"]
                    + [f"base_{m}" for m in magnitudes]
                    + ["new_status", "new_passed", "new_fresh_passed", "new_ko"]
                    + [f"new_{m}" for m in magnitudes]
                    + [f"base_class_{c}" for c in CUTOFFS] + [f"new_class_{c}" for c in CUTOFFS])
    for gene in genes:
        b, n = base.get(gene), new.get(gene)
        bg, ng = (final(b) if b else None), (final(n) if n else None)
        bc = [call(bw, bg, c) if b else "absent" for c in CUTOFFS]
        nc = [call(nw, ng, c) if n else "absent" for c in CUTOFFS]
        writer.writerow([gene, gene in positives, gene in mtdna,
                         b["status"] if b else "absent", b["passed"] if b else None,
                         (b["fresh_resolve"] or {}).get("passed") if b else None, bg]
                        + [b.get(m) if b else None for m in magnitudes]
                        + [n["status"] if n else "absent", n["passed"] if n else None,
                           (n["fresh_resolve"] or {}).get("passed") if n else None, ng]
                        + [n.get(m) if n else None for m in magnitudes] + bc + nc)
        if (b is None) != (n is None) or (bg is None) != (ng is None) or (
                bg is not None and ng is not None and abs(bg - ng) > CHANGE_TOLERANCE):
            changed_growth.append({"gene": gene, "base": bg if b else "absent", "new": ng if n else "absent"})
        for cutoff, x, y in zip(CUTOFFS, bc, nc):
            if x != y:
                changed_class[str(cutoff)].append({"gene": gene, "experimental_positive": gene in positives,
                                                   "base": x, "new": y})
def largest(rows, key):
    values = [r[key] for r in rows.values() if r.get(key) is not None]
    return {"max": max(values), "count_over_1e-9": sum(v > 1e-9 for v in values),
            "count_over_1e-7": sum(v > 1e-7 for v in values), "count_over_1e-6": sum(v > 1e-6 for v in values)}


summary = {
    "condition": "Task A medium + PO1f, GLPK, basis reset per knockout, targeted check",
    "check": {"knocked_reaction_flux_tolerance": KNOCKED_TOLERANCE, "general_tolerance": GENERAL_TOLERANCE},
    "lp_timeout_seconds": LP_TIMEOUT_SECONDS,
    "magnitudes_in_sweep": {name: {key: largest(rows, key) for key in magnitudes}
                            for name, rows in (("baseline", base), ("new", new))},
    "baseline_sha256": base_ctx.canonical_model_sha256, "new_sha256": new_ctx.canonical_model_sha256,
    "medium_sha256": base_ctx.medium_sha256, "strain_profile_sha256": base_ctx.strain_profile_sha256,
    "negative_class": "user-defined: model genes not on the consensus list",
    "validation_status": "development reference, not independent validation",
    "mtdna_out_of_scope": sorted(mtdna),
    "wt": {"baseline": base_wt, "new": new_wt},
    "genes": {"baseline": len(base), "new": len(new), "only_in_new": sorted(set(new) - set(base)),
              "only_in_baseline": sorted(set(base) - set(new))},
    "in_model_positives": {"baseline": len(positives & set(base)), "new": len(positives & set(new))},
    "sweep_points_failing_check": {
        "baseline": sorted(g for g, r in base.items() if r["status"] == "optimal" and not r["passed"]),
        "new": sorted(g for g, r in new.items() if r["status"] == "optimal" and not r["passed"])},
    "non_optimal_in_sweep": {"baseline": sorted(g for g, r in base.items() if r["status"] != "optimal"),
                             "new": sorted(g for g, r in new.items() if r["status"] != "optimal")},
    "fresh_resolves": fresh_used,
    "unresolved_after_fresh": {"baseline": sorted(g for g, r in base.items() if final(r) is None),
                               "new": sorted(g for g, r in new.items() if final(r) is None)},
    "scores": {"baseline": score(bw, base), "new": score(nw, new)},
    "growth_changed_over_1e-6": changed_growth,
    "class_changed": changed_class,
    "lp_count": {"sweep": len(base) + len(new) + 2, "fresh": len(fresh_used)},
}
with open(out_json, "w") as handle:
    json.dump(summary, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: v for k, v in summary.items() if k not in ("mtdna_out_of_scope",)}, indent=1))
