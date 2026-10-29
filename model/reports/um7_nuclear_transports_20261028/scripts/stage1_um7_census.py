"""Whole-model blocked sets (FVA at 0 % of the optimum, every reaction) for the UM7 variants, open condition
(every exchange (-1000, 1000)) and Task A + PO1f; Gurobi 1e-9. Read-only; writes OUT.json.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts:<scratchpad> \
    /opt/miniconda3/bin/python stage1_um7_census.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json VARIANTS
Blocked = max(|min|, |max|) <= 1e-9 (the census tolerance at Gurobi). Variants as in stage1_um7.py.
"""
import json, sys, warnings
import cobra
from cobra.flux_analysis import flux_variability_analysis
from checked_solve import configure
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
import importlib.util

spec = importlib.util.spec_from_file_location("s1", sys.argv[0].replace("stage1_um7_census.py", "stage1_um7_apply.py"))
s1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(s1)
model_path, medium, profile, out, variants = sys.argv[1:6]
TOL = 1e-9


def blocked(model):
    df = flux_variability_analysis(model, fraction_of_optimum=0.0, processes=1)
    return sorted(r for r in df.index if max(abs(df.loc[r, "minimum"]), abs(df.loc[r, "maximum"])) <= TOL)


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    OPEN = cobra.io.read_sbml_model(model_path)
    TASKA = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile).model
configure(OPEN); configure(TASKA)
for r in OPEN.exchanges:
    r.bounds = (-1000, 1000)
res = {}
for label in variants.split(","):
    row = {}
    for cond, M in (("open", OPEN), ("taskA", TASKA)):
        with M as m:
            s1.apply(m, label)
            row[cond] = blocked(m)
    res[label] = row
    print(label, "open", len(row["open"]), "taskA", len(row["taskA"]), flush=True)
base = res.get("as_built")
if base:
    for label, row in res.items():
        for cond in ("open", "taskA"):
            row[f"{cond}_freed_vs_as_built"] = sorted(set(base[cond]) - set(row[cond]))
            row[f"{cond}_newly_blocked_vs_as_built"] = sorted(set(row[cond]) - set(base[cond]))
        print(label, {k: v for k, v in row.items() if "vs" in k}, flush=True)
json.dump(res, open(out, "w"), indent=1)
