"""Step 6 LPs: flux for the schematic, plus two checks (Task A + PO1f, in memory).

Usage: python step6_figure_flux.py  ->  out/step6_figure_flux.json
1. pFBA at the WT optimum (Gurobi 1e-9) for: before = R1172 reopened at (-1000, 1000); after = 7ccbb568 as is;
   option = O1 (matrix 2-oxobutanoate dehydrogenase, PDH rule, inferred) + R2075 closed. Full flux vectors kept.
2. GLPK cross-check of the option: FVA (no growth floor) of the nine and the new reaction.
3. O4 (placeholder matrix CoA dilution 1.9e-5 in biomass_C): with YALI1F20788g (R834) knocked out, pFBA, and the
   CoA-moiety channels that replace R834.
Route fluxes for the figure are taken from out/step2_gurobi.json (before, R25 fixed at 1, xMAINTENANCE 0) and
out/step5b_routes.json (option, same definition); no LP is repeated for them.
"""
import cobra
from cobra.flux_analysis import pfba

from common import LP_COUNT, NINE, R1172_OLD_BOUNDS, configure, dump, fva, load

DELTA = {"R834": 1, "R1172": 1, "R1291": 1, "R2025": 1, "R2077": -1, "R2083": -1, "R2084": -1, "R2113": -1,
         "R2118": -1}


def add_obdh(model):
    r = cobra.Reaction("T_OBDHm", lower_bound=0, upper_bound=1000)
    model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in
                       {"m42[C_mi]": -1, "m73[C_mi]": -1, "m27[C_mi]": -1, "m78[C_mi]": 1, "m44[C_mi]": 1,
                        "m30[C_mi]": 1, "m28[C_mi]": 4}.items()})
    r.gene_reaction_rule = model.reactions.R693.gene_reaction_rule


def run_pfba(model):
    LP_COUNT["n"] += 2
    model.solver.problem.reset()
    sol = pfba(model)
    return {"status": sol.status, "growth": float(sol.fluxes["biomass_C"]),
            "flux": {k: float(v) for k, v in sol.fluxes.items() if abs(v) > 1e-9}}


model = load("taskA", "gurobi")
out = {"pfba": {}}
for state in ["before", "after", "option"]:
    with model:
        if state == "before":
            model.reactions.R1172.bounds = R1172_OLD_BOUNDS
        elif state == "option":
            add_obdh(model)
            model.reactions.R2075.bounds = (0, 0)
        out["pfba"][state] = run_pfba(model)
        print(state, out["pfba"][state]["status"], out["pfba"][state]["growth"], flush=True)
with model:
    model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id("m73[C_mi]"): -1.9e-5})
    model.genes.get_by_id("YALI1F20788g").knock_out()
    res = run_pfba(model)
    res["moiety_channels"] = {k: DELTA[k] * res["flux"].get(k, 0.0) for k in DELTA if abs(res["flux"].get(k, 0.0)) > 1e-12}
    res["flux"] = {k: res["flux"][k] for k in res["flux"] if k in DELTA or k in NINE or k == "biomass_C"}
    out["O4_R834_knockout"] = res
    print("O4 KO", res["status"], res["growth"], res["moiety_channels"], flush=True)
configure(model, "glpk")
with model:
    add_obdh(model)
    model.reactions.R2075.bounds = (0, 0)
    out["glpk_option_fva"] = fva(model, NINE + ["T_OBDHm"])
    print("GLPK option", {k: (v["min"], v["max"], v["min_status"], v["max_status"]) for k, v in out["glpk_option_fva"].items()})
out["lp_count"] = LP_COUNT["n"]
print("LPs", LP_COUNT["n"], "->", dump(out, "step6_figure_flux.json"))
