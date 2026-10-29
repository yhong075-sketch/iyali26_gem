"""Whole-model blocked sets (FVA at 0 % of optimum, |v| <= 1e-9) for a model with in-memory um11_apply variants.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts:SCRIPTDIR \
    /opt/miniconda3/bin/python blocked_sets.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json VARIANT[,VARIANT...]
Conditions as in stage1_um11.py: 'open' (raw model, every exchange (-1000, 1000)) and 'taskA' (Task A + PO1f).
"""
import json, sys, warnings
import cobra
from cobra.flux_analysis import flux_variability_analysis
from checked_solve import configure
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from um11_apply import apply

model_path, medium, profile, out, variants = sys.argv[1:6]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    OPEN = cobra.io.read_sbml_model(model_path)
    TASKA = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile).model
for m in (OPEN, TASKA):
    configure(m)
for r in OPEN.exchanges:
    r.bounds = (-1000, 1000)
result = {}
for label in variants.split(","):
    result[label] = {}
    for cond, model in (("open", OPEN), ("taskA", TASKA)):
        with model as m:
            apply(m, label)
            m.solver.problem.reset()
            df = flux_variability_analysis(m, fraction_of_optimum=0.0, processes=1)
            blocked = sorted(r for r in df.index if max(abs(df.loc[r, "minimum"]), abs(df.loc[r, "maximum"])) <= 1e-9)
            result[label][cond] = {"n_reactions": len(df), "blocked": blocked}
            print(label, cond, "reactions", len(df), "blocked", len(blocked), flush=True)
json.dump(result, open(out, "w"), indent=1)
