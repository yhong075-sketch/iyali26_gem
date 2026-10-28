"""Knockout sweep of the R2135 candidate with R1726 opened in memory (the paused Kennedy change), read-only.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/r2135_acyl_coa_pool_20261027/scripts/ko_sweep_r1726_stacked.py CANDIDATE.xml MEDIUM.csv PROFILE.json \
    EXPERIMENTAL.csv CANDIDATE_KO_COMPARE.tsv OUT.tsv OUT_SUMMARY.json
Same solve and check as ko_compare_two.py (checked_solve.py: Gurobi 1e-9, reset per knockout, knocked-flux
check); no fresh re-solves (a failing point is recorded as unresolved). Classes at 1 / 5 / 10 / 15 % against the
stacked WT, compared gene by gene with the candidate's own classes from the two-model sweep
(CANDIDATE_KO_COMPARE.tsv, candidate_class_* columns). This answers whether opening R1726 on top of the R2135
fix still changes any class; nothing is written to a model file.
"""
import csv
import json
import sys
import warnings
from pathlib import Path

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL
from tools.dfba_new_fn_essentiality import load_experimental_reference

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
model_path, medium, profile, experimental, previous, out_tsv, out_json = sys.argv[1:8]
positives = set(load_experimental_reference(Path(experimental))[0])
mtdna = set(json.loads((MODEL.root / "expected/mtdna_genes_out_of_scope.json").read_text())["genes"])
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
configure(model)
model.reactions.R1726.bounds = (0.0, 1000.0)
s = matrix(model)
wt = solve(model, s, set())
assert wt["passed"], wt
before = {row["gene"]: row for row in csv.DictReader(open(previous), delimiter="\t")}
rows, changed, failing = [], {str(c): [] for c in CUTOFFS}, []
for gene in sorted(g.id for g in model.genes if g.id not in ctx.excluded_runtime_genes):
    with model:
        model.genes.get_by_id(gene).knock_out()
        knocked = {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}
        point = solve(model, s, knocked)
    growth = point["growth"] if point["passed"] else None
    if not point["passed"] and point["status"] != "infeasible":
        failing.append(gene)
    classes = ["unresolved" if growth is None else ("essential" if growth / wt["growth"] < c else "non-essential")
               for c in CUTOFFS]
    rows.append([gene, gene in positives, gene in mtdna, point["status"], point["passed"], growth, *classes])
    if gene in mtdna or gene not in before:
        continue
    for cutoff, now in zip(CUTOFFS, classes):
        was = before[gene][f"candidate_class_{cutoff}"]
        if was != now:
            changed[str(cutoff)].append({"gene": gene, "experimental_positive": gene in positives,
                                         "candidate": was, "candidate_plus_r1726": now})
with open(out_tsv, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["gene", "experimental_positive", "mtdna_out_of_scope", "status", "passed", "ko_growth",
                     *[f"class_{c}" for c in CUTOFFS]])
    writer.writerows(rows)
summary = {"model_sha256": ctx.canonical_model_sha256, "r1726_bounds_in_memory": [0.0, 1000.0], "wt": wt,
           "lp_count": len(rows) + 1, "points_failing_check_not_infeasible": failing,
           "class_changed_vs_candidate": changed}
with open(out_json, "w") as handle:
    json.dump(summary, handle, indent=1)
    handle.write("\n")
print("wt", wt["growth"], "failing", failing, {k: len(v) for k, v in changed.items()}, file=sys.stderr)
