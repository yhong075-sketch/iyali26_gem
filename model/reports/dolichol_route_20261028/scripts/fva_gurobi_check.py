"""Gurobi 1e-9 FVA on named reactions under the census open condition (every exchange at (-1000, 1000)); read-only.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    fva_gurobi_check.py OUT.json RXN_IDS(comma) MODEL.xml...
Settles GLPK-only blocked calls near its 1e-7 cutoff: a reaction counts as blocked at 1e-9 when |min| and |max| < 1e-9.
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis
from cobra.io import read_sbml_model

from checked_solve import configure

out, rxns, *paths = sys.argv[1:]
rxns = rxns.split(",")
result = {}
for path in paths:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = read_sbml_model(path)
    configure(model)
    for r in model.exchanges:
        r.bounds = (-1000.0, 1000.0)
    fva = flux_variability_analysis(model, reaction_list=[model.reactions.get_by_id(r) for r in rxns], fraction_of_optimum=0.0, processes=1)
    rows = {r: {"min": float(fva.loc[r, "minimum"]), "max": float(fva.loc[r, "maximum"]),
                "blocked_1e-9": bool(abs(fva.loc[r, "minimum"]) < 1e-9 and abs(fva.loc[r, "maximum"]) < 1e-9)} for r in rxns}
    result[path.split("/")[-1][-60:]] = rows
    print(path.split("/")[-1][-60:], "blocked at 1e-9:", sorted(r for r, v in rows.items() if v["blocked_1e-9"]), "| of", len(rxns))
    for r in rxns[:3]:
        print("   ", r, rows[r])
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
