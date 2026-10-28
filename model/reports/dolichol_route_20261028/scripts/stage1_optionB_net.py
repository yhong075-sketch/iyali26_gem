"""Recommended variant with the three-step-detour net stoichiometry (TDA5 NAD+, DFG10 + ENV9 2 NADPH); read-only."""
import json, math, sys, warnings
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
model_path, medium, profile, nad, nadh, out = sys.argv[1:7]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
m = ctx.model; configure(m)
m.reactions.get_by_id("biomass_C").add_metabolites({m.metabolites.get_by_id("m503[C_er]"): -1e-4})
m.reactions.get_by_id("R2178").bounds = (0, 0)
r = Reaction("R_TEST_dedol_to_dol_net", lower_bound=0, upper_bound=1000); m.add_reactions([r])
# polyprenyl-PP + H2O -> polyprenol + PPi ; polyprenol + NAD+ -> polyprenal + NADH + H+ ; polyprenal + NADPH + H+ -> dolichal + NADP+ ; dolichal + NADPH + H+ -> dolichol + NADP+
# net: dedoldp + H2O + NAD+ + 2 NADPH + H+ -> dolichol + PPi + NADH + 2 NADP+
r.add_metabolites({m.metabolites.get_by_id(k): v for k, v in {"m2019[C_cy]": -1, "m32[C_cy]": -1, nad: -1, "m40[C_cy]": -2, "m10[C_cy]": -1, "m499[C_cy]": 1, "m203[C_cy]": 1, nadh: 1, "m41[C_cy]": 2}.items()})
st = matrix(m); wt = solve(m, st, set())
def no_uptake(obj):
    with m:
        for x in m.exchanges: x.bounds = (0.0, 1000.0)
        for x in m.reactions:
            if x.lower_bound > 0: x.lower_bound = 0.0
        m.objective = obj; m.solver.problem.reset(); v = m.slim_optimize(error_value=float("nan"))
        return {"status": m.solver.status, "value": None if math.isnan(v) else v}
res = {"wt": wt, "no_uptake_biomass": no_uptake("biomass_C"), "no_uptake_ntp1": no_uptake("R_NTP1")}
ids = ["R2179", "R2180", "R287", "R1101", "R2177", "R2181", "R2182", "R_TEST_dedol_to_dol_net"]
f = flux_variability_analysis(m, reaction_list=ids, fraction_of_optimum=0.999, processes=1)
res["fva_0.999"] = {i: [float(f.loc[i, "minimum"]), float(f.loc[i, "maximum"])] for i in ids}
# reversible test reaction: does a reversible lumped step open an energy or transhydrogenase cycle?
r.lower_bound = -1000
res["if_reversible"] = {"no_uptake_ntp1": no_uptake("R_NTP1"), "no_uptake_biomass": no_uptake("biomass_C")}
r.lower_bound = 0
kos = {}
for g in ("YALI1E22718g", "YALI1D26639g", "YALI1C26036g", "YALI1D21543g"):
    with m:
        m.genes.get_by_id(g).knock_out(); knocked = {x.id for x in m.genes.get_by_id(g).reactions if x.bounds == (0, 0)}
        k = solve(m, st, knocked)
    kos[g] = None if k["growth"] is None else k["growth"] / wt["growth"]
res["ko_ratio"] = kos
json.dump(res, open(out, "w"), indent=1)
print(json.dumps(res, indent=1))
