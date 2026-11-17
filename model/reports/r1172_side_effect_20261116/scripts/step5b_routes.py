"""Step 5b: which route carries each of the nine under each option (Task A + PO1f, Gurobi 1e-9, in memory).

Usage: python step5b_routes.py  ->  out/step5b_routes.json
Variants: O1, O1 + R2075 closed (the package the user authorized on the other line), O3, O4, with the same
definitions as step5_options.py (the variant code is repeated here because that script runs on import).
For O1 + R2075 closed also WT, FVA of the nine and the energy check. Route = xMAINTENANCE lower bound 0, target
fixed at min(1, half its maximum), minimum total absolute flux (as step 2 part C). Each route records the
CoA-moiety channels across the inner membrane, the gene-less reactions it uses, and the direction of R2112.
"""
import cobra
from optlang.symbolics import Zero

from common import LP_COUNT, NINE, TOL, dump, extreme, flux_vector, fva, load, lp

DELTA = {"R834": 1, "R1172": 1, "R1291": 1, "R2025": 1, "R2077": -1, "R2083": -1, "R2084": -1, "R2113": -1,
         "R2118": -1}


def add_rxn(model, rid, stoich, bounds, rule):
    r = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoich.items()})
    r.gene_reaction_rule = rule
    return r


def apply(model, variant):
    if variant.startswith("O1"):
        add_rxn(model, "T_OBDHm", {"m42[C_mi]": -1, "m73[C_mi]": -1, "m27[C_mi]": -1, "m78[C_mi]": 1,
                                    "m44[C_mi]": 1, "m30[C_mi]": 1, "m28[C_mi]": 4}, (0, 1000),
                model.reactions.R693.gene_reaction_rule)
        if variant == "O1_R2075closed":
            model.reactions.R2075.bounds = (0, 0)
    elif variant == "O3":
        for comp in ("C_cy", "C_mi"):
            model.add_metabolites([cobra.Metabolite(f"T_pcarn[{comp}]", formula="C10H19NO4", charge=0,
                                                    name="O-propanoylcarnitine", compartment=comp)])
        add_rxn(model, "T_CRPTc", {"m1970[C_cy]": -1, "m369[C_cy]": -1, "m10[C_cy]": -4, "T_pcarn[C_cy]": 1,
                                    "m69[C_cy]": 1}, (-1000, 1000), model.reactions.R191.gene_reaction_rule)
        add_rxn(model, "T_PCARt", {"T_pcarn[C_cy]": -1, "m374[C_mi]": -1, "m369[C_cy]": 1, "T_pcarn[C_mi]": 1},
                (0, 1000), model.reactions.R825.gene_reaction_rule)
        add_rxn(model, "T_CRPTm", {"T_pcarn[C_mi]": -1, "m73[C_mi]": -1, "m78[C_mi]": 1, "m374[C_mi]": 1,
                                    "m28[C_mi]": 4}, (-1000, 1000), model.reactions.R193.gene_reaction_rule)
    elif variant == "O4":
        model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id("m73[C_mi]"): -1.9e-5})


def min_total_flux(model, target, value):
    with model:
        model.reactions.get_by_id(target).bounds = (value, value)
        coefficients = {}
        for r in model.reactions:
            coefficients[r.forward_variable] = 1.0
            coefficients[r.reverse_variable] = 1.0
        model.objective = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective.set_linear_coefficients(coefficients)
        total, status = lp(model)
        flux = flux_vector(model) if total is not None else None
    return total, status, flux


def energy_check(model):
    out = {}
    with model:
        for r in model.boundary:
            r.lower_bound = max(r.lower_bound, 0.0)
        for r in model.reactions:
            if r.lower_bound > 0:
                r.lower_bound = 0.0
        for obj in ("R_NTP1", "biomass_C"):
            with model:
                model.objective = obj
                out[obj] = dict(zip(("value", "status"), lp(model)))
    return out


model = load("taskA", "gurobi")
maxima = {}  # maxima of the nine per variant from step5_options.json are not reused: recomputed only for the new one
import json
from common import ART
prior = json.load(open(ART / "out" / "step5_options.json"))["variants"]
result = {}
for variant in ["O1", "O1_R2075closed", "O3", "O4"]:
    with model:
        apply(model, variant)
        res = {}
        if variant == "O1_R2075closed":
            res["wt"] = dict(zip(("value", "status"), lp(model)))
            res["fva"] = fva(model, NINE + ["T_OBDHm"])
            res["energy"] = energy_check(model)
            ranges = res["fva"]
        else:
            ranges = prior[variant]["fva"]
        model.reactions.xMAINTENANCE.lower_bound = 0.0
        routes = {}
        for target in NINE:
            hi, lo = ranges[target]["max"], ranges[target]["min"]
            value = min(1.0, 0.5 * hi) if hi > TOL else (max(-1.0, 0.5 * lo) if lo < -TOL else 0.0)
            if value == 0.0:
                routes[target] = {"note": "blocked in this variant"}
                continue
            total, status, flux = min_total_flux(model, target, value)
            entry = {"fixed_at": value, "status": status, "total_abs_flux": total}
            if flux is not None:
                nz = {k: v for k, v in flux.items() if abs(v) > 1e-9}
                entry["nonzero"] = dict(sorted(nz.items()))
                entry["geneless_used"] = sorted(k for k in nz if not model.reactions.get_by_id(k).gene_reaction_rule
                                                and not model.reactions.get_by_id(k).boundary)
                entry["moiety_channels"] = {k: DELTA[k] * flux[k] for k in DELTA if abs(flux[k]) > 1e-9}
                entry["R2112"] = flux["R2112"]
            routes[target] = entry
        res["routes"] = routes
        result[variant] = res
        print(variant, res.get("wt"), res.get("energy"),
              {t: (r.get("status"), r.get("moiety_channels"), round(r.get("R2112", 0), 4)) for t, r in routes.items()},
              flush=True)
result["lp_count"] = LP_COUNT["n"]
print("LPs", LP_COUNT["n"], "->", dump(result, "step5b_routes.json"))
