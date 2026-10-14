"""Task C, Amendment 2: checks required before using the R825-reversible ASSUMPTION (diagnostic, in memory).

1. Loop check, R825 reversible (lb -1000), at the V2 elastic-fit state (overlay, 22 ranges) and at WT
   (V2 with and without overlay; PO1f + CSM-Leu):
   - FVA of R825, R1275, R1205 (standard and cobra loopless=True, i.e. CycleFreeFlux per extremum);
   - cycle flux in the representative pFBA vector: min(-v825, v1275, v1205) with v825 < 0;
   - the fit total / WT growth with the loop blocked (R1275 = 0, R1205 = 0, both);
   - full-model loopless MILP (cobra add_loopless) of the fit LP (reported as supplementary: it removes
     every internal cycle in the model, not only this one).
   Plus complex I KO/WT under the PDH floor (V2, overlay, 13C caps), loop open vs blocked, because
   Task D will use the assumption there.
2. Carnitine: maximum demand of each carnitine species (no uptake allowed) on V2, V1 and CSM, with R825
   as supplied and reversible; every non-transport reaction producing a carnitine species; name scan for
   carnitine-biosynthesis intermediates; carnitine exchange bounds in each medium.

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/r825_loop_carnitine_check.py
Writes r825_loop_carnitine_check.json here. Never overwrites.
"""
import json
import time
from datetime import datetime, timezone

from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.flux_analysis.loopless import add_loopless

from common import (HERE, MU, UPTAKE, Counter_, add_caps, apply_overlay, check_pins, growth, load, load_caps,
                    pdh_floor, r825_reversible)
from elastic_fit import MAPPING, ranges

LOOP = ("R825", "R1275", "R1205")
CARNITINE = ("m369[C_cy]", "m374[C_mi]", "m371[C_pe]", "m1069[C_ex]")
CARNITINE_FORMS = ("m369[C_cy]", "m374[C_mi]", "m371[C_pe]", "m370[C_cy]", "m372[C_pe]", "m373[C_mi]", "m1069[C_ex]")
BLOCKS = {"none": (), "R1275": ("R1275",), "R1205": ("R1205",), "both": ("R1275", "R1205")}


def per100(x):
    return None if x is None else 100.0 * x / UPTAKE


def add_fit(model):
    """Elastic-fit slacks (as elastic_fit.fit); returns the total-violation expression."""
    rng = ranges()
    prob = model.problem
    model.reactions.biomass_C.bounds = (MU, MU)
    extra, terms = [], []
    for key, _, _, coeffs in MAPPING:
        lb, ub, _ = rng[key]
        e = (100.0 / UPTAKE) * sum(c * model.reactions.get_by_id(r).flux_expression for r, c in coeffs.items())
        sp, sm = prob.Variable(f"sp_{key}", lb=0), prob.Variable(f"sm_{key}", lb=0)
        extra += [sp, sm, prob.Constraint(e - sp, ub=ub), prob.Constraint(e + sm, lb=lb)]
        terms.append((sp + sm) / (ub - lb))
    model.add_cons_vars(extra)
    return sum(terms)


def solve_fit(model):
    total = add_fit(model)
    model.objective = model.problem.Objective(total, direction="min")
    Counter_.lps += 1
    model.slim_optimize()
    return (model.solver.objective.value if model.solver.status == "optimal" else None), model.solver.status, total


def loop_fva(model):
    out = {}
    for loopless in (False, True):
        Counter_.lps += 1 + 2 * len(LOOP)
        try:
            f = flux_variability_analysis(model, reaction_list=list(LOOP), fraction_of_optimum=1.0, loopless=loopless,
                                          processes=1)
            out["loopless" if loopless else "standard"] = {r: [per100(f.loc[r, "minimum"]), per100(f.loc[r, "maximum"])]
                                                            for r in LOOP}
        except Exception as exc:
            out["loopless" if loopless else "standard"] = {"error": f"{type(exc).__name__}: {exc}"}
    return out


def cycle_flux(vec):
    v825, v1275, v1205 = vec["R825"], vec["R1275"], vec["R1205"]
    return per100(max(0.0, min(-v825, v1275, v1205)))


def main():
    inputs = check_pins()
    caps = load_caps()
    started = time.monotonic()
    record = {"started_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
              "command": "python -B model/reports/w29_13c_consistency_20261012/r825_loop_carnitine_check.py"}
    models = {m: load(m) for m in ("V2", "V1", "CSM")}

    # 1a. Fit state, V2, overlay, R825 reversible.
    sim, model = models["V2"]
    fit = {}
    for name, blocked in BLOCKS.items():
        with model:
            apply_overlay(model)
            r825_reversible(model)
            for rid in blocked:
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            opt, status, total = solve_fit(model)
            entry = {"blocked": list(blocked), "total_violation": opt, "status": status}
            if name == "none" and opt is not None:
                model.add_cons_vars([model.problem.Constraint(total, ub=opt + 1e-6, name="face")])
                model.objective = "biomass_C"
                entry["fva_per_100"] = loop_fva(model)
                Counter_.lps += 2
                vec = pfba(model).fluxes
                entry["pfba_point_per_100"] = {r: per100(vec[r]) for r in LOOP}
                entry["cycle_flux_in_pfba_point_per_100"] = cycle_flux(vec)
            fit[name] = entry
        print("fit", name, fit[name]["total_violation"], flush=True)
    # Full-model loopless MILP of the fit (supplementary).
    with model:
        apply_overlay(model)
        r825_reversible(model)
        model.solver.problem.Params.TimeLimit = 600.0
        add_loopless(model)
        t0 = time.monotonic()
        opt, status, _ = solve_fit(model)
        fit["full_loopless_milp"] = {"total_violation": opt, "status": status, "seconds": time.monotonic() - t0,
                                     "note": "add_loopless removes every internal cycle in the model"}
        model.solver.problem.Params.TimeLimit = 60.0
    print("fit loopless MILP", fit["full_loopless_milp"], flush=True)
    record["fit_state_V2"] = fit

    # 1b. WT states.
    wt = {}
    for medium, overlay in (("V2", True), ("V2", False), ("CSM", False)):
        sim, model = models[medium]
        key = f"{medium}|overlay={overlay}"
        wt[key] = {}
        for name, blocked in BLOCKS.items():
            with model:
                if overlay:
                    apply_overlay(model)
                r825_reversible(model)
                for rid in blocked:
                    model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                g, s = growth(model)
                entry = {"blocked": list(blocked), "growth": g, "status": s}
                if name == "none" and g:
                    entry["fva_per_100"] = loop_fva(model)
                    Counter_.lps += 2
                    vec = pfba(model).fluxes
                    entry["pfba_point_per_100"] = {r: per100(vec[r]) for r in LOOP}
                    entry["cycle_flux_in_pfba_point_per_100"] = cycle_flux(vec)
                wt[key][name] = entry
            print("wt", key, name, wt[key][name]["growth"], flush=True)
    record["wt_states"] = wt

    # 1c. Complex I KO/WT under the PDH floor (V2, overlay, 13C caps), loop open vs blocked.
    sim, model = models["V2"]
    ko = {}
    for name, blocked in BLOCKS.items():
        res = {}
        for label, rid_ko in (("WT", None), ("KO", "R2062")):
            with model:
                apply_overlay(model)
                add_caps(model, caps, "module")
                pdh_floor(model)
                r825_reversible(model)
                for rid in blocked:
                    model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                if rid_ko:
                    model.reactions.get_by_id(rid_ko).bounds = (0.0, 0.0)
                res[label] = growth(model)
        ko[name] = {"blocked": list(blocked), "WT": res["WT"][0], "KO": res["KO"][0], "KO_status": res["KO"][1],
                    "KO_over_WT": None if not res["WT"][0] or res["KO"][0] is None else res["KO"][0] / res["WT"][0]}
        print("ko", name, ko[name], flush=True)
    record["complex_I_floor_caps_V2"] = ko

    # 2. Carnitine.
    carn = {"exchange_bounds": {}, "max_demand_per_gDW_h": {}, "producers": {}, "name_scan": []}
    for medium, (sim, model) in models.items():
        carn["exchange_bounds"][medium] = list(model.reactions.R943.bounds)
        for rev in (False, True):
            for mid in CARNITINE:
                with model:
                    if rev:
                        r825_reversible(model)
                    dm = model.add_boundary(model.metabolites.get_by_id(mid), type="demand", reaction_id=f"DM_test_{mid}")
                    model.objective = dm
                    Counter_.lps += 1
                    val = model.slim_optimize()
                    carn["max_demand_per_gDW_h"][f"{medium}|R825={'rev' if rev else 'asis'}|{mid}"] = \
                        val if model.solver.status == "optimal" else model.solver.status
    sim, model = models["V2"]
    for mid in CARNITINE_FORMS:
        met = model.metabolites.get_by_id(mid)
        carn["producers"][mid] = sorted(
            {r.id: {"name": r.name, "gpr": r.gene_reaction_rule or None, "bounds": list(r.bounds),
                    "compartments": sorted({x.compartment for x in r.metabolites})}
             for r in met.reactions}.items())
    for x in model.metabolites:
        n = x.name.lower()
        if any(k in n for k in ("trimethyl", "butyrobetaine", "carnitine")):
            carn["name_scan"].append([x.id, x.name])
    carn["in_biomass"] = [x.id for x in model.reactions.biomass_C.metabolites if "carnitine" in x.name.lower()]
    record["carnitine"] = carn
    print("carnitine max demand", carn["max_demand_per_gDW_h"], flush=True)

    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started)
    (HERE / "r825_loop_carnitine_check.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
