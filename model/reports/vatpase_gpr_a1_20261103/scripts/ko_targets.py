"""Verify-1: WT and the 15 V-ATPase rule-gene knockouts on one model (16 LPs, no sweep).

Same checked_solve.py as the earlier rounds (Gurobi at 1e-9, a model reset before each solve, the targeted check) and
the same simulation context loader as ../../ceramide_pool_biomass_20261101/scripts/ko_screen_certified.py. No
certificate step: a non-optimal knockout is a stop condition here (TASK.md section 4), not a scored point.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ko_targets.py MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.json
"""
import json
import sys
import warnings
from pathlib import Path

from checked_solve import GENERAL_TOLERANCE, KNOCKED_TOLERANCE, SOLVER_TOLERANCE, configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.vatpase_a1_gpr import A_SUBUNITS, CORE, NEW_GENE, REACTIONS
from tools.dfba_new_fn_essentiality import load_experimental_reference

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
model_path, medium, profile, experimental, out_json = sys.argv[1:6]
positives = set(load_experimental_reference(Path(experimental))[0])
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
configure(model)
stoichiometry = matrix(model)
wt = solve(model, stoichiometry, ())
wt["flux"] = {rid: model.solver.primal_values[model.reactions.get_by_id(rid).forward_variable.name]
              - model.solver.primal_values[model.reactions.get_by_id(rid).reverse_variable.name] for rid in REACTIONS}
if not wt["passed"]:
    raise RuntimeError(f"The WT point failed the check ({wt['status']}); stopping as declared")
rows = {}
for gene in CORE + (NEW_GENE,) + A_SUBUNITS:
    with model:
        model.genes.get_by_id(gene).knock_out()
        knocked = sorted(r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0))
        row = solve(model, stoichiometry, set(knocked))
    row["closed"] = knocked
    row["ratio"] = row["growth"] / wt["growth"] if row["growth"] is not None else None
    row["experimental_positive"] = gene in positives
    row["class"] = {str(c): (None if row["ratio"] is None else "essential" if row["ratio"] < c else "non-essential")
                    for c in CUTOFFS}
    rows[gene] = row
    print(gene, row["status"], row["passed"], row["growth"], row["closed"], file=sys.stderr, flush=True)
summary = {
    "condition": f"Task A medium + PO1f, Gurobi FeasibilityTol/OptimalityTol {SOLVER_TOLERANCE}, reset per solve, targeted check",
    "check": {"knocked_reaction_flux_tolerance": KNOCKED_TOLERANCE, "general_tolerance": GENERAL_TOLERANCE},
    "model_sha256": ctx.canonical_model_sha256, "medium_sha256": ctx.medium_sha256,
    "strain_profile_sha256": ctx.strain_profile_sha256,
    "excluded_runtime_genes": sorted(set(rows) & set(ctx.excluded_runtime_genes)),
    "negative_class": "user-defined: model genes not on the consensus list",
    "validation_status": "development reference, not independent validation",
    "wt": wt, "knockouts": rows, "lp_count": 1 + len(rows),
    "all_passed": all(r["passed"] for r in rows.values()),
    "closed_r794_exactly_core_and_e": sorted(g for g, r in rows.items() if "R794" in r["closed"]) == sorted(CORE + (NEW_GENE,)),
}
Path(out_json).write_text(json.dumps(summary, indent=1) + "\n")
print(json.dumps({k: summary[k] for k in ("model_sha256", "lp_count", "all_passed", "closed_r794_exactly_core_and_e")}))
