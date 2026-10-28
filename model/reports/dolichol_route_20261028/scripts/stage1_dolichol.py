"""Stage 1 of the dolichol-chain task: in-memory LP checks (read-only, no file written except OUT.json).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    stage1_dolichol.py MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.json
Task A medium + PO1f, Gurobi at 1e-9 (checked_solve.py). Variants (applied in memory):
  as_built            : nothing
  dolp_demand         : irreversible demand (0, 1000) on dolichyl phosphate m503[C_er]
  dolp_biomass_4e-5   : m503[C_er] in biomass_C at 4e-5 mmol/gDW
  dolp_biomass_1.6e-4 : same at 1.6e-4
  dfg10_route_4e-5    : biomass 4e-5 + R2178 closed + test reaction dedoldp -> dolichol (topology test only:
                        m2019 + H2O + NADPH + H+ -> m499 + diphosphate + NADP+; chemistry not claimed)
  h5_um6_4e-5         : biomass 4e-5 + census test H5 (IPP and FPP cy<->lp links, sink on docosaprenyl-PP m432)
Per variant: WT growth (checked); FVA of the route reactions at 0 and 99.9 % of the optimum; pFBA flux at the
optimum; EGC LP (every exchange (0, 1000), positive lower bounds relaxed, max R_NTP1, status); targeted knockouts.
"""
import json, math, sys, warnings
from pathlib import Path
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis, pfba
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.dfba_new_fn_essentiality import load_experimental_reference

ROUTE = [f"R{i}" for i in range(217, 236)] + ["R2177", "R2178", "R2179", "R2180", "R2181", "R2182", "R287", "R1101",
         "R288", "R2102", "R2106", "R2105", "R2138", "R289", "R372", "R374", "R375", "R2101", "R173", "R141", "R322", "R295"]
GENES = ["YALI1C26036g", "YALI1D21543g", "YALI1D26639g", "YALI1E22718g", "YALI1D21022g", "YALI1E06759g", "YALI1F06018g",
         "YALI1F08363g", "YALI1B21618g", "YALI1B15180g", "YALI1E39395g", "YALI1D07696g", "YALI1D08050g", "YALI1B07538g",
         "YALI1A19526g", "YALI1C32184g", "YALI1E07121g", "YALI1E18121g", "YALI1D13201g"]
DOLP_ER = "m503[C_er]"
model_path, medium, profile, experimental, out = sys.argv[1:6]
positives = set(load_experimental_reference(Path(experimental))[0])
lp_count = 0


def context():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


def add_rxn(model, rid, name, mets, bounds, gpr=""):
    r = Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in mets.items()})
    r.gene_reaction_rule = gpr
    return r


def biomass_term(model, coeff):
    model.reactions.get_by_id("biomass_C").add_metabolites({model.metabolites.get_by_id(DOLP_ER): -coeff})


def apply(model, label):
    if label == "as_built":
        return
    if label == "dolp_demand":
        add_rxn(model, "DM_TEST_dolp_er", "test demand dolichyl phosphate [er]", {DOLP_ER: -1}, (0, 1000))
        return
    if label == "dolp_biomass_4e-5":
        biomass_term(model, 4e-5); return
    if label == "dolp_biomass_1.6e-4":
        biomass_term(model, 1.6e-4); return
    if label == "dfg10_route_4e-5":
        biomass_term(model, 4e-5)
        model.reactions.get_by_id("R2178").bounds = (0, 0)
        add_rxn(model, "R_TEST_dedol_to_dol", "test: dehydrodolichol-PP -> dolichol (lumped phosphatase + reductase)",
                {"m2019[C_cy]": -1, "m32[C_cy]": -1, "m40[C_cy]": -1, "m10[C_cy]": -1, "m499[C_cy]": 1, "m203[C_cy]": 1, "m41[C_cy]": 1},
                (0, 1000), "YALI1E09357g")
        return
    if label == "h5_um6_4e-5":
        biomass_term(model, 4e-5)
        add_rxn(model, "R_TEST_ipp_lp", "test IPP cy<->lp", {"m491[C_cy]": -1, "m412[C_lp]": 1}, (-1000, 1000))
        add_rxn(model, "R_TEST_fpp_lp", "test FPP cy<->lp", {"m512[C_cy]": -1, "m411[C_lp]": 1}, (-1000, 1000))
        add_rxn(model, "R_TEST_sink_m432", "test sink docosaprenyl-PP lp", {"m432[C_lp]": -1}, (-1000, 1000))
        return
    raise ValueError(label)


def checked(model, stoich, knocked=()):
    global lp_count
    lp_count += 1
    return solve(model, stoich, set(knocked))


def egc(model):
    global lp_count
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = []
        for r in model.reactions:
            if r.lower_bound > 0:
                relaxed.append([r.id, r.lower_bound]); r.lower_bound = 0.0
        model.objective = "R_NTP1"
        model.solver.problem.reset(); lp_count += 1
        v = model.slim_optimize(error_value=float("nan"))
        return {"status": model.solver.status, "max_R_NTP1_no_uptake": None if math.isnan(v) else v, "relaxed": relaxed}


def fva(model, ids):
    global lp_count
    rows = {}
    for fraction in (0.0, 0.999):
        lp_count += 2 * len(ids) + 1
        frame = flux_variability_analysis(model, reaction_list=ids, fraction_of_optimum=fraction, processes=1)
        rows[str(fraction)] = {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in ids}
    return rows


result = {"model": model_path, "variants": {}}
VARIANTS = ["as_built", "dolp_demand", "dolp_biomass_4e-5", "dolp_biomass_1.6e-4", "dfg10_route_4e-5", "h5_um6_4e-5"]
for label in VARIANTS:
    ctx = context(); model = ctx.model
    result.setdefault("model_sha256", ctx.canonical_model_sha256)
    apply(model, label)
    ids = [r for r in ROUTE if r in model.reactions] + [r.id for r in model.reactions if r.id.endswith("_TEST") or "_TEST_" in r.id]
    stoich = matrix(model)
    wt = checked(model, stoich)
    row = {"wt": wt, "egc": egc(model)}
    row["fva"] = fva(model, ids)
    row["blocked_route_reactions_at_0"] = sorted(r for r, (lo, hi) in row["fva"]["0.0"].items() if abs(lo) < 1e-9 and abs(hi) < 1e-9)
    model.solver.problem.reset(); lp_count += 1
    sol = pfba(model)
    row["pfba_route_flux"] = {r: float(sol.fluxes[r]) for r in ids if abs(sol.fluxes[r]) > 1e-12}
    row["pfba_dolp_consumers"] = {r.id: float(sol.fluxes[r.id]) for r in model.metabolites.get_by_id(DOLP_ER).reactions if abs(sol.fluxes[r.id]) > 1e-12}
    kos = {}
    for g in GENES:
        if g not in model.genes:
            kos[g] = "not in model"; continue
        with model:
            model.genes.get_by_id(g).knock_out()
            knocked = {r.id for r in model.genes.get_by_id(g).reactions if r.bounds == (0, 0)}
            k = checked(model, stoich, knocked)
        k["ratio"] = None if k["growth"] is None or not wt["growth"] else k["growth"] / wt["growth"]
        k["on_consensus_list"] = g in positives
        kos[g] = k
    row["knockouts"] = kos
    result["variants"][label] = row
    print(label, wt["status"], wt["growth"], "egc", row["egc"]["status"], row["egc"]["max_R_NTP1_no_uptake"],
          "blocked:", len(row["blocked_route_reactions_at_0"]), "/", len(ids), flush=True)
    for g, k in kos.items():
        if isinstance(k, dict) and k["ratio"] is not None and k["ratio"] < 0.15:
            print("   lethal:", g, k["ratio"], "consensus" if k["on_consensus_list"] else "", flush=True)
result["lp_count"] = lp_count
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp_count)
