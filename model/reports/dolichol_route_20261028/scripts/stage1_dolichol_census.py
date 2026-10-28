"""Open-condition blocked scan (GLPK, census method) for the dolichol Stage 1 variants. Read-only.
Usage: python stage1_dolichol_census.py MODEL.xml OUT.json   (PYTHONPATH = census scripts dir)"""
import json, sys
from common import load_model, apply_condition, blocked
from cobra import Reaction
model_path, out = sys.argv[1:3]
DOLP_ER = "m503[C_er]"
def add_rxn(M, rid, mets, bounds):
    r = Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1]); M.add_reactions([r])
    r.add_metabolites({M.metabolites.get_by_id(k): v for k, v in mets.items()})
def apply(M, label):
    if label == "dolp_demand": add_rxn(M, "DM_TEST_dolp_er", {DOLP_ER: -1}, (0, 1000))
    elif label == "dolp_biomass_4e-5": M.reactions.get_by_id("biomass_C").add_metabolites({M.metabolites.get_by_id(DOLP_ER): -4e-5})
    elif label == "dfg10_route_4e-5":
        M.reactions.get_by_id("biomass_C").add_metabolites({M.metabolites.get_by_id(DOLP_ER): -4e-5})
        M.reactions.get_by_id("R2178").bounds = (0, 0)
        add_rxn(M, "R_TEST_dedol_to_dol", {"m2019[C_cy]": -1, "m32[C_cy]": -1, "m40[C_cy]": -1, "m10[C_cy]": -1, "m499[C_cy]": 1, "m203[C_cy]": 1, "m41[C_cy]": 1}, (0, 1000))
    elif label == "h5_um6_4e-5":
        M.reactions.get_by_id("biomass_C").add_metabolites({M.metabolites.get_by_id(DOLP_ER): -4e-5})
        add_rxn(M, "R_TEST_ipp_lp", {"m491[C_cy]": -1, "m412[C_lp]": 1}, (-1000, 1000))
        add_rxn(M, "R_TEST_fpp_lp", {"m512[C_cy]": -1, "m411[C_lp]": 1}, (-1000, 1000))
        add_rxn(M, "R_TEST_sink_m432", {"m432[C_lp]": -1}, (-1000, 1000))
res = {}
base = None
for label in ["as_built", "dolp_demand", "dolp_biomass_4e-5", "dfg10_route_4e-5", "h5_um6_4e-5"]:
    M = load_model(model_path, "glpk"); apply_condition(M, None); apply(M, label)
    b, note = blocked(M)
    b = {x for x in b if not x.startswith("R_TEST") and not x.startswith("DM_TEST")}
    if base is None: base = b
    res[label] = {"n_blocked": len(b), "note": note, "freed_vs_as_built": sorted(base - b), "newly_blocked_vs_as_built": sorted(b - base)}
    print(label, len(b), "freed", len(base - b), "newly blocked", len(b - base), note, flush=True)
json.dump(res, open(out, "w"), indent=1)
