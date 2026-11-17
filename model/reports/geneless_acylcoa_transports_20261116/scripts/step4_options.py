"""Step 4: options for the eight reactions, tested in memory on 7ccbb568 (Task A + PO1f). Nothing is written to SBML.

Usage: python step4_options.py   ->  out/step4_options.json

Bases
- B0: 7ccbb568 as is.
- B1: B0 + T_OBDHm, a matrix 2-oxobutanoate dehydrogenase (m42 + CoA + NAD -> propionyl-CoA + CO2 + NADH + 4 H+,
      bounds (0, 1000), rule of R693, INFERRED), with R2075 closed: the package of the methylcitrate line (candidate
      b24a5e78 on cbf2ef95, uncommitted) and option "O1 + R2075 closed" of the R1172 round.
Variants (applied on a base inside `with model:`)
- K  : keep (control).
- M  : merge the cytosolic propionyl-CoA m1875 (no formula) into m1970 (C24H36N7O17P3S, -4). R2024 then makes
       m1970; it closes in H and charge only with 3 H+ on the product side, which is added (both close, O unaffected).
       R2025 becomes an exact duplicate of R2113 (m1970[cy] <=> m78[mi]) and is set to (0, 0) as the duplicate.
- MC : M + R2113 closed: no propionyl-CoA crosses the inner membrane.
- C6 : the six blocked reactions closed: R2077, R2083, R2084, R2118, R1291, R2117.
- V  : the valine steps R2117, R2116, R2115, R2114 and R2112 moved into the matrix (same gene rules and bounds),
       with R2118 closed. New matrix species (formulas INFERRED from KEGG and the model's matrix acyl-CoA form):
       methacrylyl-CoA C25H40N7O17P3S 0, (S)-3-hydroxyisobutyryl-CoA C25H42N7O18P3S 0 (as m1973),
       (S)-3-hydroxyisobutyrate C4H7O3 -1, (S)-methylmalonate semialdehyde C4H5O3 -1. H+ on each moved reaction is
       set so that H and charge both close; the script stops if a moved reaction does not balance.
- ALL: M + MC + V + R2077, R2083, R2084, R1291 closed.
C6 on B0 is not solved: the six are blocked on B0 (step 1, both solvers), so closing them leaves the feasible set
unchanged (provable), and every LP result equals B0's.
Per variant (Gurobi): WT; FVA without growth floor of the TARGETS present; energy check; targeted knockouts on the
variants in KO_VARIANTS. GLPK repeats WT and the energy check for every variant.
"""
import cobra

from common import (CUTOFFS, KO_GENES, LP_COUNT, TOL, classify, dump, energy_check, fva, load, lp)

TARGETS = ["R2025", "R2113", "R2077", "R2083", "R2084", "R2118", "R1291", "R2117", "R2117_mi", "R2114",
           "R2114_mi", "R2112", "R2112_mi", "R25", "R95", "R552", "R834", "R2074", "R2075", "R2024", "R2022",
           "R2285", "R1969", "R1292", "R2082", "T_OBDHm"]
VARIANTS = [("B0", "K"), ("B0", "M"), ("B0", "MC"), ("B0", "V"), ("B0", "ALL"),
            ("B1", "K"), ("B1", "M"), ("B1", "MC"), ("B1", "C6"), ("B1", "V"), ("B1", "ALL")]
KO_VARIANTS = {("B0", "M"), ("B0", "V"), ("B0", "ALL"), ("B1", "K"), ("B1", "MC"), ("B1", "ALL")}
NO_FVA = {("B0", "K")}  # step 1 has the no-floor FVA of the as-is model

NEW_MI = {  # id: (name, formula, charge)
    "V_macoa[C_mi]": ("methacrylyl-CoA (in memory)", "C25H40N7O17P3S", 0),
    "V_hibcoa[C_mi]": ("(S)-3-hydroxyisobutyryl-CoA (in memory)", "C25H42N7O18P3S", 0),
    "V_hib[C_mi]": ("(S)-3-hydroxyisobutyrate (in memory)", "C4H7O3", -1),
    "V_mmsa[C_mi]": ("(S)-methylmalonate semialdehyde (in memory)", "C4H5O3", -1),
}
# cytosolic species -> matrix species for the moved valine steps; H+ (m10 -> m28) is recomputed
MAP = {"m1971[C_cy]": "m1829[C_mi]", "m528[C_cy]": "m563[C_mi]", "m1972[C_cy]": "V_macoa[C_mi]",
       "m32[C_cy]": "m26[C_mi]", "m1973[C_cy]": "V_hibcoa[C_mi]", "m1974[C_cy]": "V_hib[C_mi]",
       "m69[C_cy]": "m73[C_mi]", "m122[C_cy]": "m27[C_mi]", "m123[C_cy]": "m30[C_mi]",
       "m1975[C_cy]": "V_mmsa[C_mi]", "m1970[C_cy]": "m78[C_mi]", "m82[C_cy]": "m44[C_mi]",
       "m566[C_mi]": "m566[C_mi]"}
MOVED = ["R2117", "R2116", "R2115", "R2114", "R2112"]


def met(model, mid):
    return model.metabolites.get_by_id(mid)


def add_rxn(model, rid, stoich, bounds, rule, name=""):
    r = cobra.Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([r])
    r.add_metabolites({met(model, k): v for k, v in stoich.items()})
    r.gene_reaction_rule = rule
    return r


def close_h(rxn, proton_id):
    """Add H+ so that H and charge both close; stop if anything else is off or the two disagree."""
    bal = rxn.check_mass_balance()
    if not bal:
        return 0
    extra = {k: v for k, v in bal.items() if k not in ("H", "charge")}
    h, q = bal.get("H", 0), bal.get("charge", 0)
    if extra or abs(h - q) > 1e-9:
        raise SystemExit(f"{rxn.id} does not close with H+ alone: {bal}")
    rxn.add_metabolites({met(rxn.model, proton_id): -h})
    if rxn.check_mass_balance():
        raise SystemExit(f"{rxn.id} still unbalanced: {rxn.check_mass_balance()}")
    return -h


def apply_base(model, base):
    log = {}
    if base == "B1":
        r = add_rxn(model, "T_OBDHm", {"m42[C_mi]": -1, "m73[C_mi]": -1, "m27[C_mi]": -1, "m78[C_mi]": 1,
                                       "m44[C_mi]": 1, "m30[C_mi]": 1, "m28[C_mi]": 4}, (0, 1000),
                    model.reactions.R693.gene_reaction_rule, "2-oxobutanoate dehydrogenase (in memory, inferred)")
        model.reactions.R2075.bounds = (0, 0)
        log["T_OBDHm"] = {"equation": r.build_reaction_string(), "balance": r.check_mass_balance(),
                          "rule": r.gene_reaction_rule}
    return log


def merge(model, log):
    old, new = met(model, "m1875[C_cy]"), met(model, "m1970[C_cy]")
    for rid in ["R2024", "R2025"]:
        r = model.reactions.get_by_id(rid)
        c = r.metabolites[old]
        r.add_metabolites({old: -c, new: c})
    log["R2024_H_added"] = close_h(model.reactions.R2024, "m10[C_cy]")
    log["R2024"] = {"equation": model.reactions.R2024.build_reaction_string(),
                    "balance": model.reactions.R2024.check_mass_balance()}
    log["R2025_after_merge"] = model.reactions.R2025.build_reaction_string()
    model.reactions.R2025.bounds = (0, 0)  # exact duplicate of R2113 after the merge
    if len(old.reactions) != 0:
        raise SystemExit(f"m1875 still used by {[r.id for r in old.reactions]}")


def valine_to_matrix(model, log):
    for mid, (name, formula, charge) in NEW_MI.items():
        model.add_metabolites([cobra.Metabolite(mid, name=name, formula=formula, charge=charge,
                                                compartment="C_mi")])
    moved = {}
    for rid in MOVED:
        r = model.reactions.get_by_id(rid)
        stoich = {}
        for m, c in r.metabolites.items():
            if m.id in ("m10[C_cy]",):
                continue  # H+ recomputed below
            if m.id not in MAP:
                raise SystemExit(f"{rid}: no matrix map for {m.id}")
            stoich[MAP[m.id]] = stoich.get(MAP[m.id], 0) + c
        nr = add_rxn(model, f"{rid}_mi", stoich, r.bounds, r.gene_reaction_rule, f"{r.name} (moved to C_mi)")
        h = close_h(nr, "m28[C_mi]")
        r.bounds = (0, 0)
        moved[rid] = {"equation": nr.build_reaction_string(), "H_added": h, "rule": nr.gene_reaction_rule,
                      "bounds": list(nr.bounds), "balance": nr.check_mass_balance()}
    model.reactions.R2118.bounds = (0, 0)
    log["moved"] = moved


def apply_variant(model, variant):
    log = {}
    if variant in ("M", "MC", "ALL"):
        merge(model, log)
    if variant in ("MC", "ALL"):
        model.reactions.R2113.bounds = (0, 0)
    if variant == "C6":
        for rid in ["R2077", "R2083", "R2084", "R2118", "R1291", "R2117"]:
            model.reactions.get_by_id(rid).bounds = (0, 0)
    if variant in ("V", "ALL"):
        valine_to_matrix(model, log)
    if variant == "ALL":
        for rid in ["R2077", "R2083", "R2084", "R1291"]:
            model.reactions.get_by_id(rid).bounds = (0, 0)
    return log


def run_variants(model, result):
    for base, variant in VARIANTS:
        key = f"{base}/{variant}"
        with model:
            res = {"edits": apply_base(model, base)}
            res["edits"].update(apply_variant(model, variant))
            wt, s = lp(model)
            res["wt"], res["wt_status"] = wt, s
            if (base, variant) not in NO_FVA:
                res["fva"] = fva(model, [t for t in TARGETS if t in model.reactions])
            res["energy"] = energy_check(model)
            if (base, variant) in KO_VARIANTS:
                kos = {}
                for gid in KO_GENES:
                    with model:
                        model.genes.get_by_id(gid).knock_out()
                        g, gs = lp(model)
                        if g is None:  # one fresh re-solve; never read a failed solve as 0
                            g, gs = lp(model)
                    kos[gid] = {"growth": g, "status": gs, "ratio": None if (g is None or not wt) else g / wt,
                                "class": {str(c): x for c, x in classify(g, wt).items()}}
                res["ko"] = kos
            result["variants"][key] = res
        blocked = sorted(k for k, v in res.get("fva", {}).items() if v.get("blocked"))
        print(key, "wt", wt, s, "| energy", {k: v["value"] for k, v in res["energy"].items()},
              "| blocked", blocked, "| LPs", LP_COUNT["n"], flush=True)


def run_glpk(result):
    """GLPK cross-check of WT and the energy check."""
    gmodel = load("taskA", "glpk")
    for base, variant in VARIANTS:
        with gmodel:
            apply_base(gmodel, base)
            apply_variant(gmodel, variant)
            wt, s = lp(gmodel)
            result["glpk"][f"{base}/{variant}"] = {"wt": wt, "wt_status": s, "energy": energy_check(gmodel)}
        print("glpk", base, variant, wt, s, result["glpk"][f"{base}/{variant}"]["energy"], flush=True)


def main():
    result = {"variants": {}, "glpk": {}, "ko_genes": KO_GENES, "cutoffs": list(CUTOFFS)}
    model = load("taskA", "gurobi")
    run_variants(model, result)
    run_glpk(result)
    result["lp_count"] = LP_COUNT["n"]
    print("LPs", LP_COUNT["n"], "->", dump(result, "step4_options.json"))


if __name__ == "__main__":  # guard added after the run (logic unchanged; as-run copy in scripts_as_run/)
    main()
