"""Checked single-gene knockouts on one model, scored under the earlier rule (A) and the certified-infeasible rule (B).

One-model copy of ../../aconitase_aco2_split_20261029/scripts/ko_compare_two.py (bb3118a6, itself unchanged from the
dolichol round): same checked_solve.py (Gurobi at 1e-9, a model reset before each knockout, the targeted check),
same sweep loop, FRESH_LIMIT and timeout. Added (TASK.md, approved 2026-10-30): every knockout whose sweep status is
infeasible gets a certificate in a fresh context:
  1. the re-solve must again be infeasible;
  2. xMAINTENANCE must be the only reaction with lb > 0 or ub < 0;
  3. with the xMAINTENANCE lower bound set to 0, maximizing xMAINTENANCE must be optimal, pass the targeted check,
     and stay below the original lower bound by more than GENERAL_TOLERANCE;
  4. WT must be optimal (asserted for the whole run).
Rule A: every non-optimal knockout is unresolved. Rule B: a certified knockout scores as growth 0.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ko_screen_certified.py MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv UNTARGETED.json OUT.tsv OUT_SUMMARY.json
Class: essential when KO/WT is strictly below the cutoff. Negatives are the user-defined class; mtDNA genes are out
of scope; scores are given with and without the untargeted-gene scope; the list was used in development, so scores
are a development reference.
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
FRESH_LIMIT = 20
MAGNITUDES = ("max_knocked_flux", "max_bound_violation", "max_residual")
MAINTENANCE = "xMAINTENANCE"

model_path, medium, profile, experimental, untargeted_path, out_tsv, out_json = sys.argv[1:8]
positives = set(load_experimental_reference(Path(experimental))[0])
mtdna = set(json.loads((MODEL.root / "expected/mtdna_genes_out_of_scope.json").read_text())["genes"])
untargeted = set(json.loads(Path(untargeted_path).read_text())["genes"])
SCOPES = {"mtdna": mtdna, "mtdna_untargeted": mtdna | untargeted}
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


def forced(model):
    return sorted(r.id for r in model.reactions if r.lower_bound > 0 or r.upper_bound < 0)


def ceiling(model, stoichiometry, knocked):
    """Maximum maintenance flux with its lower bound set to 0 (the 'growth' field holds the objective value)."""
    with model:
        model.reactions.get_by_id(MAINTENANCE).lower_bound = 0
        model.objective = MAINTENANCE
        return solve(model, stoichiometry, knocked)


def certify(path, gene):
    ctx = context(path)
    model = ctx.model
    stoichiometry = matrix(model)
    with model:
        model.genes.get_by_id(gene).knock_out()
        knocked = closed(model, gene)
        resolve = solve(model, stoichiometry, knocked)
        forced_now = forced(model)
        bound = model.reactions.get_by_id(MAINTENANCE).lower_bound
        top = ceiling(model, stoichiometry, knocked)
    margin = bound - top["growth"] if top["status"] == "optimal" and top["growth"] is not None else None
    conditions = {
        "1_resolve_infeasible": resolve["status"] == "infeasible",
        "2_only_maintenance_forced": forced_now == [MAINTENANCE],
        "3_ceiling_optimal_passed_below_bound": (top["passed"] and margin is not None and margin > GENERAL_TOLERANCE),
    }
    return {"resolve_status": resolve["status"], "forced": forced_now, "maintenance_lower_bound": bound,
            "ceiling": top, "max_maintenance": top["growth"], "margin": margin, "conditions": conditions,
            "certified": all(conditions.values())}


def run(path):
    ctx = context(path)
    model = ctx.model
    stoichiometry = matrix(model)
    wt = solve(model, stoichiometry, ())
    if not wt["passed"]:
        raise RuntimeError(f"The WT point failed the check ({wt['status']}); stopping before the sweep as declared")
    wt_ceiling = ceiling(model, stoichiometry, ())
    rows = {}
    for index, gene in enumerate(sorted(g.id for g in model.genes if g.id not in ctx.excluded_runtime_genes), 1):
        started = time.monotonic()
        with model:
            model.genes.get_by_id(gene).knock_out()
            row = solve(model, stoichiometry, closed(model, gene))
        row["fresh_resolve"] = None
        row["certificate"] = None
        if not row["passed"] and row["status"] != "infeasible":
            if len(fresh_used) >= FRESH_LIMIT:
                raise RuntimeError(f"More than {FRESH_LIMIT} knockouts need a fresh re-solve; stopping as declared")
            fresh_used.append(gene)
            row["fresh_resolve"] = fresh(path, gene)
        if row["status"] == "infeasible":
            row["certificate"] = certify(path, gene)
        rows[gene] = row
        print(f"{ctx.canonical_model_sha256[:8]} {index} {gene} {row['status']} passed={row['passed']}"
              f" fresh={None if row['fresh_resolve'] is None else row['fresh_resolve']['passed']}"
              f" certified={None if row['certificate'] is None else row['certificate']['certified']}"
              f" {time.monotonic() - started:.2f}s", file=sys.stderr, flush=True)
    return ctx, wt, wt_ceiling, rows


def final(row, rule):
    for point in (row, row.get("fresh_resolve") or {}):
        if point.get("passed"):
            return point["growth"]
    if rule == "B" and row["certificate"] and row["certificate"]["certified"]:
        return 0.0
    return None


def call(wt, growth, cutoff):
    if growth is None or not wt:
        return "unresolved"
    return "essential" if growth / wt < cutoff else "non-essential"


def score(wt, rows, rule, excluded):
    result = {}
    for cutoff in CUTOFFS:
        counts = {"TP": 0, "FN": 0, "FP": 0, "TN": 0, "unresolved": 0}
        for gene, row in rows.items():
            if gene in excluded:
                continue
            c = call(wt, final(row, rule), cutoff)
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


ctx, wt, wt_ceiling, rows = run(model_path)
wt_growth = wt["growth"]
certified = sorted(g for g, r in rows.items() if r["certificate"] and r["certificate"]["certified"])
changed = sorted(g for g, r in rows.items()
                 if any(call(wt_growth, final(r, "A"), c) != call(wt_growth, final(r, "B"), c) for c in CUTOFFS))

with open(out_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["gene", "experimental_positive", "mtdna_out_of_scope", "untargeted_out_of_scope", "status",
                     "passed", "fresh_passed", "ko", *MAGNITUDES, "cert_resolve_status", "cert_forced",
                     "cert_maintenance_lower_bound", "cert_max_maintenance", "cert_margin", "cert_ceiling_passed",
                     "cert_certified", *[f"class_A_{c}" for c in CUTOFFS], *[f"class_B_{c}" for c in CUTOFFS]])
    for gene, r in sorted(rows.items()):
        cert = r["certificate"] or {}
        writer.writerow([gene, gene in positives, gene in mtdna, gene in untargeted, r["status"], r["passed"],
                         (r["fresh_resolve"] or {}).get("passed"), final(r, "A"), *[r.get(m) for m in MAGNITUDES],
                         cert.get("resolve_status"), ";".join(cert.get("forced", [])) or None,
                         cert.get("maintenance_lower_bound"), cert.get("max_maintenance"), cert.get("margin"),
                         (cert.get("ceiling") or {}).get("passed"), cert.get("certified"),
                         *[call(wt_growth, final(r, "A"), c) for c in CUTOFFS],
                         *[call(wt_growth, final(r, "B"), c) for c in CUTOFFS]])
summary = {
    "condition": f"Task A medium + PO1f, Gurobi FeasibilityTol/OptimalityTol {SOLVER_TOLERANCE}, reset per knockout, targeted check",
    "check": {"knocked_reaction_flux_tolerance": KNOCKED_TOLERANCE, "general_tolerance": GENERAL_TOLERANCE},
    "lp_timeout_seconds": LP_TIMEOUT_SECONDS,
    "model": {"sha256": ctx.canonical_model_sha256, "genes": len(rows)},
    "medium_sha256": ctx.medium_sha256, "strain_profile_sha256": ctx.strain_profile_sha256,
    "negative_class": "user-defined: model genes not on the consensus list",
    "validation_status": "development reference, not independent validation",
    "rules": {"A": "every non-optimal knockout unresolved (earlier rule)",
              "B": "certified maintenance-shortfall infeasible knockout scored as growth 0 (benchmark.md, decided 2026-10-30)"},
    "scopes": {"mtdna": "mtDNA genes out of scope", "mtdna_untargeted": "mtDNA and untargeted genes out of scope",
               "untargeted_genes": sorted(untargeted)},
    "wt": wt,
    "wt_maintenance_ceiling": wt_ceiling,
    "in_model_positives": len(positives & set(rows)),
    "magnitudes_in_sweep": {key: largest(rows, key) for key in MAGNITUDES},
    "sweep_points_failing_check": sorted(g for g, r in rows.items() if r["status"] == "optimal" and not r["passed"]),
    "non_optimal_in_sweep": sorted(g for g, r in rows.items() if r["status"] != "optimal"),
    "fresh_resolves": fresh_used,
    "certificates": {g: r["certificate"] for g, r in sorted(rows.items()) if r["certificate"]},
    "certified": certified,
    "not_certified_infeasible": sorted(g for g, r in rows.items() if r["certificate"] and not r["certificate"]["certified"]),
    "class_changed_A_to_B": changed,
    "scores": {f"{rule}_{scope}": score(wt_growth, rows, rule, excluded)
               for rule in ("A", "B") for scope, excluded in SCOPES.items()},
    "lp_count": {"sweep": len(rows) + 1, "wt_ceiling": 1, "fresh": len(fresh_used),
                 "certificate": 2 * sum(1 for r in rows.values() if r["certificate"])},
}
with open(out_json, "w") as handle:
    json.dump(summary, handle, indent=1)
    handle.write("\n")
print(json.dumps({k: summary[k] for k in ("wt", "certified", "not_certified_infeasible", "class_changed_A_to_B",
                                          "fresh_resolves", "lp_count")}, indent=1))
