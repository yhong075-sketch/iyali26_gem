"""In-memory R24 gene-rule variants and targeted knockouts (read-only; no model file is written).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python ../model/reports/r24_homocitrate_review_20261029/scripts/r24_variants.py \
    MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.json

Each variant changes only R24's gene rule (R395 stays YALI1E03387g). WT and the knockouts below use the
same model, medium, PO1f profile and objective; Gurobi at 1e-9 with a reset before every solve and the
targeted check (checked_solve.py, unchanged). Only genes in R24's old or new rule can change class, because
a variant changes no other rule, bound or stoichiometry.
"""
import json
import sys
import warnings
from pathlib import Path

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.dfba_new_fn_essentiality import load_experimental_reference

LYS4, ACO1, ACO2 = "YALI1E03387g", "YALI1D11984g", "YALI1E17935g"
VARIANTS = {
    "A_status_quo": None,
    "B_aco1_only": ACO1,
    "C_aco1_or_aco2": f"{ACO1} or {ACO2}",
    "D_aco2_only": ACO2,
    "E_lys4_or_aco1_or_aco2": f"{LYS4} or {ACO1} or {ACO2}",
}
KNOCKOUTS = ((LYS4,), (ACO1,), (ACO2,), (ACO1, ACO2))
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
REPORTED = ("R24", "R395", "R1170", "R398", "R238", "R216", "R1220", "R1202", "R1204")

model_path, medium, profile, experimental, out_json = sys.argv[1:6]
positives = set(load_experimental_reference(Path(experimental))[0])
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=Path(model_path), media_path=medium,
                                            strain_profile_path=profile)
model = ctx.model
configure(model)
stoichiometry = matrix(model)


def fluxes():
    primal = model.solver.primal_values
    return {rid: primal[model.reactions.get_by_id(rid).forward_variable.name]
            - primal[model.reactions.get_by_id(rid).reverse_variable.name] for rid in REPORTED}


def closed(genes):
    return {r.id for g in genes for r in model.genes.get_by_id(g).reactions if r.bounds == (0, 0)}


result = {"model_sha256": ctx.canonical_model_sha256, "variants": {}}
original_rule = model.reactions.R24.gene_reaction_rule
for name, rule in VARIANTS.items():
    with model:
        if rule is not None:
            model.reactions.R24.gene_reaction_rule = rule
        wt = solve(model, stoichiometry, ())
        wt["fluxes"] = fluxes() if wt["passed"] else None
        rows = {}
        for genes in KNOCKOUTS:
            with model:
                for g in genes:
                    model.genes.get_by_id(g).knock_out()
                knocked = closed(genes)
                row = solve(model, stoichiometry, knocked)
                row["closed_reactions"] = sorted(knocked)
                row["fluxes"] = fluxes() if row["passed"] else None
            ratio = row["growth"] / wt["growth"] if row["passed"] and wt["passed"] else None
            row["ko_wt"] = ratio
            label = "+".join(genes)
            if len(genes) == 1:
                positive = genes[0] in positives
                row["experimental_positive"] = positive
                row["class"] = {str(c): None if ratio is None else (
                    ("TP" if positive else "FP") if ratio < c else ("FN" if positive else "TN")) for c in CUTOFFS}
            rows[label] = row
        result["variants"][name] = {"R24_rule": model.reactions.R24.gene_reaction_rule, "wt": wt, "knockouts": rows}
    assert model.reactions.R24.gene_reaction_rule == original_rule
Path(out_json).write_text(json.dumps(result, indent=1))
for name, v in result["variants"].items():
    print(name, "|", v["R24_rule"], "| WT", v["wt"]["growth"])
    for label, row in v["knockouts"].items():
        f = row["fluxes"] or {}
        print(f"   {label:27s} {row['status']:8s} passed={row['passed']} growth={row['growth']} ko/wt={row['ko_wt']}"
              f" class={row.get('class')} R24={f.get('R24')} R1220_lys={f.get('R1220')}")
