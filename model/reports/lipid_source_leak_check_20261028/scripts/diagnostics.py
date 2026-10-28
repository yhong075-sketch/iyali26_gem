"""Follow-up diagnostics after source_check.py (read-only; at most 10 LPs, inside the 300-LP cap).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/lipid_source_leak_check_20261028/scripts/diagnostics.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
D1 C1 with the 11 P exchanges closed: max Ins1P demand (does R2171 make phosphate from nothing?).
D2 same: max Pi demand, repeated, then with R2171 closed.
D3 C5 WT with the 15 formula-less exchanges left open; D4 C5 WT as run (formula-less closed), to see which one
   kills growth only if D3 grows.
"""
import json
import sys
import warnings

from cobra import Reaction

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, medium, profile, out = sys.argv[1:5]
P_EX = ["R950", "R951", "R1005", "R1014", "R1102", "R1126", "R1269", "R1296", "R1322", "R1339", "R1341"]
FORMULA_LESS = ["R947", "R1067", "R1189", "R1294", "R1984", "R1986", "R1990", "R1992", "R1998", "R2018", "R2026",
                "R2031", "R2036", "R2214", "R2275"]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile).model
configure(model)
res, n = {}, 0


def demand(mid, closed=()):
    global n
    with model:
        for rid in closed:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        dm = Reaction(f"TESTDM_{mid}", lower_bound=0.0, upper_bound=1000.0)
        dm.add_metabolites({model.metabolites.get_by_id(mid): -1.0})
        model.add_reactions([dm])
        model.objective = dm.id
        n += 1
        r = solve(model, matrix(model), set(closed))
        if r["growth"] is not None:
            pv = model.solver.primal_values
            r["fluxes"] = {rid: pv[model.reactions.get_by_id(rid).forward_variable.name] - pv[model.reactions.get_by_id(rid).reverse_variable.name]
                           for rid in ("R2171", "R558", "R1263", "R1296", "R914")}
        return r


with model:
    for rid in P_EX:
        model.reactions.get_by_id(rid).lower_bound = max(0.0, model.reactions.get_by_id(rid).lower_bound)
    res["D1_ins1p_no_P"] = demand("m813[C_cy]")
    res["D2_pi_no_P"] = demand("m35[C_cy]")
    res["D2_pi_no_P_R2171_closed"] = demand("m35[C_cy]", ["R2171"])
with model:
    for r in model.exchanges:
        r.bounds = (-1000.0, 1000.0)
    n += 1
    res["D3_C5_formula_less_open"] = solve(model, matrix(model), set())
    for rid in FORMULA_LESS:
        model.reactions.get_by_id(rid).lower_bound = 0.0
    n += 1
    res["D4_C5_formula_less_closed"] = solve(model, matrix(model), set())
res["lp_count"] = n
json.dump(res, open(out, "w"), indent=1)
for k, v in res.items():
    print(k, v if not isinstance(v, dict) else (v.get("status"), v.get("growth"), v.get("passed"), v.get("fluxes")), file=sys.stderr)
