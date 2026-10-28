"""Checks on the built files: baseline 00e8e397 vs the R1726 candidate (read-only, about 40 LPs).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/kennedy_closures_20261027/scripts/built_checks.py BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Task A medium + PO1f, Gurobi 1e-9 (checked_solve.py). Per model: WT (checked); FVA of R1726 / R1725 at 100 %,
99.9 % and growth >= 0 (processes=1); no-uptake LPs with every exchange (0, 1000) and every positive lower
bound relaxed to 0 (max R_NTP1 and max biomass_C, status recorded); the PSD1 + PSD2 double knockout; the PLD
gene YALI1E22633g knockout; the mechanism control with choline R1041 and ethanolamine R1115 uptake 1 and the
PE-methylation genes knocked out one at a time; the targeted R1726-route fluxes of a pFBA point.
"""
import hashlib
import json
import math
import sys
import warnings
from pathlib import Path

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

ROUTE = ["R1726", "R212", "R213", "R1603", "R1042", "R1806", "R1842", "R1743", "R1776", "R1725", "R629", "R1841",
         "R1835", "R1845", "R1844"]
baseline_path, candidate_path, medium, profile, out = sys.argv[1:6]


def context(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


def no_uptake(model, objective):
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = [[r.id, r.lower_bound] for r in model.reactions if r.lower_bound > 0]
        for rid, _ in relaxed:
            model.reactions.get_by_id(rid).lower_bound = 0.0
        model.objective = objective
        model.solver.problem.reset()
        value = model.slim_optimize(error_value=float("nan"))
        return {"objective": objective, "status": model.solver.status,
                "value": None if math.isnan(value) else value, "relaxed": relaxed}


def knocked(model, genes, reactions=()):
    with model:
        for gene in genes:
            model.genes.get_by_id(gene).knock_out()
        for rid in reactions:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        closed = {r.id for g in genes for r in model.genes.get_by_id(g).reactions if r.bounds == (0, 0)} | set(reactions)
        return solve(model, matrix(model), closed)


def exchange_of(model, mid):
    found = [r for r in model.metabolites.get_by_id(mid).reactions if r.boundary]
    assert len(found) == 1, (mid, [r.id for r in found])
    return found[0]


result = {}
for label, path in (("baseline", baseline_path), ("candidate", candidate_path)):
    ctx = context(path)
    model = ctx.model
    row = {"path": path, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
           "bounds": {rid: list(model.reactions.get_by_id(rid).bounds) for rid in ("R1726", "R1725")}}
    row["wt"] = solve(model, matrix(model), set())
    row["fva"] = {}
    for fraction in (1.0, 0.999, 0.0):
        frame = flux_variability_analysis(model, reaction_list=["R1726", "R1725"], fraction_of_optimum=fraction, processes=1)
        row["fva"][str(fraction)] = {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])]
                                     for rid in ("R1726", "R1725")}
    row["no_uptake"] = [no_uptake(model, "R_NTP1"), no_uptake(model, "biomass_C")]
    row["psd_double_ko"] = knocked(model, ("YALI1D27274g", "YALI1D04432g"))
    row["pld_gene_ko"] = knocked(model, ("YALI1E22633g",))
    fluxes = pfba(model).fluxes
    row["pfba_route"] = {rid: float(fluxes[rid]) for rid in ROUTE}
    with model:
        for mid in ("m1037[C_ex]", "m1047[C_ex]"):
            exchange_of(model, mid).lower_bound = -1.0
        row["choline_ethanolamine_control"] = {
            "wt": solve(model, matrix(model), set()),
            "YALI1E07347g_ko": knocked(model, ("YALI1E07347g",)),
            "YALI1E15432g_ko": knocked(model, ("YALI1E15432g",))}
    result[label] = row
    print(label, row["sha256"][:8], row["bounds"], row["wt"]["status"], row["wt"]["growth"], row["wt"]["passed"],
          [(r["objective"], r["status"], r["value"]) for r in row["no_uptake"]], file=sys.stderr, flush=True)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
