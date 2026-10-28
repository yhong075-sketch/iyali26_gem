"""Stage 1 follow-up: the recommended variant at the proposed coefficient, plus controls (read-only)."""
import json, math, sys, warnings
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
model_path, medium, profile, out = sys.argv[1:5]
DOLP_ER = "m503[C_er]"
def context():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model); return ctx
def add_lumped(model, gpr=""):
    r = Reaction("R_TEST_dedol_to_dol", lower_bound=0, upper_bound=1000); model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in {"m2019[C_cy]": -1, "m32[C_cy]": -1, "m40[C_cy]": -1, "m10[C_cy]": -1, "m499[C_cy]": 1, "m203[C_cy]": 1, "m41[C_cy]": 1}.items()})
    r.gene_reaction_rule = gpr
def build(label, coeff):
    ctx = context(); m = ctx.model
    m.reactions.get_by_id("biomass_C").add_metabolites({m.metabolites.get_by_id(DOLP_ER): -coeff})
    m.reactions.get_by_id("R2178").bounds = (0, 0)
    if label != "no_new_step": add_lumped(m)
    return m
def no_uptake(model, objective):
    with model:
        for r in model.exchanges: r.bounds = (0.0, 1000.0)
        for r in model.reactions:
            if r.lower_bound > 0: r.lower_bound = 0.0
        model.objective = objective; model.solver.problem.reset()
        v = model.slim_optimize(error_value=float("nan"))
        return {"status": model.solver.status, "value": None if math.isnan(v) else v}
res = {}
for label, coeff in (("optionB_1e-4", 1e-4), ("optionB_4e-5", 4e-5), ("optionB_1.6e-4", 1.6e-4), ("no_new_step", 1e-4)):
    m = build(label, coeff); st = matrix(m)
    wt = solve(m, st, set())
    row = {"coefficient": coeff, "wt": wt, "no_uptake_biomass": no_uptake(m, "biomass_C"), "no_uptake_ntp1": no_uptake(m, "R_NTP1")}
    if wt["growth"]:
        ids = ["R2179", "R2180", "R287", "R1101", "R2177", "R2181", "R2182", "R2102", "R289", "R_TEST_dedol_to_dol"]
        ids = [i for i in ids if i in m.reactions]
        f = flux_variability_analysis(m, reaction_list=ids, fraction_of_optimum=0.999, processes=1)
        row["fva_0.999"] = {i: [float(f.loc[i, "minimum"]), float(f.loc[i, "maximum"])] for i in ids}
        kos = {}
        for g in ("YALI1E22718g", "YALI1D26639g", "YALI1C26036g", "YALI1D21543g", "YALI1E06759g", "YALI1D21022g"):
            with m:
                m.genes.get_by_id(g).knock_out()
                knocked = {r.id for r in m.genes.get_by_id(g).reactions if r.bounds == (0, 0)}
                k = solve(m, st, knocked)
            kos[g] = None if k["growth"] is None else k["growth"] / wt["growth"]
        row["ko_ratio"] = kos
    res[label] = row
    print(label, wt["status"], wt["growth"], "no-uptake biomass", row["no_uptake_biomass"], "ntp1", row["no_uptake_ntp1"], row.get("ko_ratio"), flush=True)
json.dump(res, open(out, "w"), indent=1)
