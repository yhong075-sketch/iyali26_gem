"""Step 5: in-memory tests of the options on 7ccbb568 (Task A + PO1f, Gurobi 1e-9). Nothing is written to SBML.

Usage: python step5_options.py   ->  out/step5_options.json

Variants (each applied inside `with model:` on the as-is model, R1172 stays closed):
- B   : as is (control).
- O1  : matrix 2-oxobutanoate dehydrogenase, m42 + CoA + NAD -> propionyl-CoA + CO2 + NADH + 4 H+ (C_mi),
        bounds (0, 1000), gene rule = R693 (pyruvate dehydrogenase), INFERRED. Same chemistry as R_OBDHm of the
        methylcitrate round on the other line.
- O1b : as O1 with the branched-chain 2-oxo acid dehydrogenase rule of R1969 instead (INFERRED).
- O3  : propionyl-carnitine shuttle: R191-like cytosolic transfer (m1970 + carnitine + 4 H+ <=> propionylcarnitine
        + CoA, rule of R191), R825-like carrier (propionylcarnitine[cy] + carnitine[mi] -> carnitine[cy] +
        propionylcarnitine[mi], rule of R825), R193-like matrix transfer (propionylcarnitine + CoA <=>
        propionyl-CoA + carnitine + 4 H+, rule of R193). New species O-propanoylcarnitine C10H19NO4, charge 0, in
        C_cy and C_mi. Rules INFERRED (acetyl enzymes assumed to take propionyl).
- O4  : matrix CoA growth dilution: biomass_C also consumes 1.9e-5 m73[C_mi] (10 % of its cytosolic CoA
        coefficient). PLACEHOLDER value, no evidence for it; tests only what such a term would do.
Per variant: added reactions' mass balance; WT growth; FVA (no growth floor) of the nine and the added reactions;
energy check (every boundary lower bound and the xMAINTENANCE lower bound set to 0, maximise R_NTP1 and biomass_C,
status kept); targeted knockouts (gene.knock_out in a context, growth, KO/WT, class at 15 %).
"""
import cobra

from common import LP_COUNT, NINE, dump, fva, load, lp

GENES = sorted(set(
    ["YALI1E01019g", "YALI1E03300g", "YALI1C24124g", "YALI1F39620g", "YALI1F03803g", "YALI1B02194g",
     "YALI1C04093g", "YALI1E00588g", "YALI1F31215g", "YALI1C02550g", "YALI1F20788g"]  # genes of the nine
    + ["YALI1E31947g", "YALI1F27556g", "YALI1D26360g", "YALI1B13179g", "YALI1D30618g"]  # PDH rule (O1)
    + ["YALI1D11159g", "YALI1F07630g", "YALI1D30866g"]  # BCKDH rule (O1b), E3 shared with PDH
    + ["YALI1F28064g", "YALI1B13948g", "YALI1C03415g"]  # carnitine shuttle rules (O3)
    + ["YALI1D03047g", "YALI0E18568g", "YALI1F13202g", "YALI1F08845g"]))  # Thr deaminase, R2285, R265, R2024


def add_rxn(model, rid, stoich, bounds, rule):
    r = cobra.Reaction(rid, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in stoich.items()})
    r.gene_reaction_rule = rule
    return r


def apply(model, variant):
    added = []
    if variant in ("O1", "O1b"):
        rule = model.reactions.R693.gene_reaction_rule if variant == "O1" else model.reactions.R1969.gene_reaction_rule
        added.append(add_rxn(model, "T_OBDHm", {"m42[C_mi]": -1, "m73[C_mi]": -1, "m27[C_mi]": -1, "m78[C_mi]": 1,
                                                 "m44[C_mi]": 1, "m30[C_mi]": 1, "m28[C_mi]": 4}, (0, 1000), rule))
    elif variant == "O3":
        pc = {}
        for comp in ("C_cy", "C_mi"):
            mt = cobra.Metabolite(f"T_pcarn[{comp}]", formula="C10H19NO4", charge=0, name="O-propanoylcarnitine",
                                  compartment=comp)
            model.add_metabolites([mt])
            pc[comp] = mt.id
        added.append(add_rxn(model, "T_CRPTc", {"m1970[C_cy]": -1, "m369[C_cy]": -1, "m10[C_cy]": -4,
                                                 pc["C_cy"]: 1, "m69[C_cy]": 1}, (-1000, 1000),
                             model.reactions.R191.gene_reaction_rule))
        added.append(add_rxn(model, "T_PCARt", {pc["C_cy"]: -1, "m374[C_mi]": -1, "m369[C_cy]": 1, pc["C_mi"]: 1},
                             (0, 1000), model.reactions.R825.gene_reaction_rule))
        added.append(add_rxn(model, "T_CRPTm", {pc["C_mi"]: -1, "m73[C_mi]": -1, "m78[C_mi]": 1, "m374[C_mi]": 1,
                                                 "m28[C_mi]": 4}, (-1000, 1000),
                             model.reactions.R193.gene_reaction_rule))
    elif variant == "O4":
        model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id("m73[C_mi]"): -1.9e-5})
    return added


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
                v, s = lp(model)
            out[obj] = {"value": v, "status": s}
    return out


model = load("taskA", "gurobi")
result = {"genes": GENES, "variants": {}}
for variant in ["B", "O1", "O1b", "O3", "O4"]:
    with model:
        added = apply(model, variant)
        res = {"added": {r.id: {"equation": r.build_reaction_string(), "rule": r.gene_reaction_rule,
                                "bounds": list(r.bounds), "mass_balance": r.check_mass_balance()} for r in added}}
        wt, s = lp(model)
        res["wt"], res["wt_status"] = wt, s
        targets = NINE + [r.id for r in added]
        if variant != "O1b":  # O1b is the same LP as O1; only the gene rule differs
            res["fva"] = fva(model, targets)
            res["energy"] = energy_check(model)
        kos = {}
        genes = GENES if variant != "O1b" else ["YALI1E31947g", "YALI1F27556g", "YALI1D26360g", "YALI1B13179g",
                                                 "YALI1D30618g", "YALI1D11159g", "YALI1F07630g", "YALI1D30866g"]
        for gid in genes:
            with model:
                model.genes.get_by_id(gid).knock_out()
                g, gs = lp(model)
            ratio = (g / wt) if (g is not None and wt) else None
            kos[gid] = {"growth": g, "status": gs, "ko_over_wt": ratio,
                        "class_15": None if ratio is None else ("essential" if ratio < 0.15 else "non-essential")}
        res["ko"] = kos
        result["variants"][variant] = res
        print(variant, "wt", wt, s, "| blocked of nine:",
              [r for r in NINE if res.get("fva", {}).get(r, {}).get("blocked")] if "fva" in res else "n/a",
              "| energy", res.get("energy"), "| added", {k: v["mass_balance"] for k, v in res["added"].items()},
              flush=True)
result["lp_count"] = LP_COUNT["n"]
print("LPs", LP_COUNT["n"], "->", dump(result, "step5_options.json"))
