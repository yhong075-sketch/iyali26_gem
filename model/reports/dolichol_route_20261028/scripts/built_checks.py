"""Checks on the built files: baseline ce9e1112 vs the dolichol-route candidate (read-only, about 60 LPs).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/dolichol_route_20261028/scripts/built_checks.py BASELINE.xml CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
CANDIDATE may be the byte-identical seed-1 copy: cobra on macOS rejects absolute paths over 260 characters.
Task A medium + PO1f, Gurobi 1e-9 (checked_solve.py). Per model: WT (checked); FVA of the route at 100 %,
99.9 % and growth >= 0 (processes=1); no-uptake LPs with every exchange (0, 1000) and every positive lower bound
relaxed to 0 (max R_NTP1 and max biomass_C, status recorded); on the candidate the same no-uptake LPs with
R_DOLSYN made reversible in memory; the route fluxes of a pFBA point; targeted knockouts (SEC59, CWH8, the two
cis-prenyltransferase genes alone and together, the Dol-P-mannose synthase genes); mass balance of the touched
reactions as built.
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

ROUTE = ["R2179", "R2180", "R2178", "R_DOLSYN", "R287", "R1101", "R2102", "R289", "R2177", "R2181", "R2182", "R322"]
KNOCKOUTS = {"SEC59 YALI1E22718g": ("YALI1E22718g",), "CWH8 YALI1D26639g": ("YALI1D26639g",),
             "RER2-type YALI1C26036g": ("YALI1C26036g",), "cis-PT YALI1D21543g": ("YALI1D21543g",),
             "both cis-PT": ("YALI1C26036g", "YALI1D21543g"),
             "DPM YALI1D08050g": ("YALI1D08050g",), "DPM YALI1B07538g": ("YALI1B07538g",), "DPM YALI1A19526g": ("YALI1A19526g",)}
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


def knocked(model, genes):
    with model:
        for gene in genes:
            model.genes.get_by_id(gene).knock_out()
        closed = {r.id for g in genes for r in model.genes.get_by_id(g).reactions if r.bounds == (0, 0)}
        return solve(model, matrix(model), closed)


result = {}
for label, path in (("baseline", baseline_path), ("candidate", candidate_path)):
    ctx = context(path)
    model = ctx.model
    ids = [rid for rid in ROUTE if rid in model.reactions]
    row = {"path": path, "sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
           "canonical_sha256": ctx.canonical_model_sha256,
           "bounds": {rid: list(model.reactions.get_by_id(rid).bounds) for rid in ids},
           "biomass_m503": model.reactions.biomass_C.metabolites.get(model.metabolites.get_by_id("m503[C_er]"))}
    row["wt"] = solve(model, matrix(model), set())
    row["fva"] = {}
    for fraction in (1.0, 0.999, 0.0):
        frame = flux_variability_analysis(model, reaction_list=ids, fraction_of_optimum=fraction, processes=1)
        row["fva"][str(fraction)] = {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in ids}
    row["no_uptake"] = [no_uptake(model, "R_NTP1"), no_uptake(model, "biomass_C")]
    if "R_DOLSYN" in model.reactions:
        with model:
            model.reactions.R_DOLSYN.lower_bound = -1000.0
            row["no_uptake_dolsyn_reversible"] = [no_uptake(model, "R_NTP1"), no_uptake(model, "biomass_C")]
    fluxes = pfba(model).fluxes
    row["pfba_route"] = {rid: float(fluxes[rid]) for rid in ids}
    row["knockouts"] = {}
    for name, genes in KNOCKOUTS.items():
        point = knocked(model, genes)
        point["ratio"] = None if point["growth"] is None else point["growth"] / row["wt"]["growth"]
        row["knockouts"][name] = point
    row["mass_balance_as_built"] = {}
    for rid in ("R2178", "R_DOLSYN", "R287"):
        if rid in model.reactions:
            try:
                row["mass_balance_as_built"][rid] = model.reactions.get_by_id(rid).check_mass_balance()
            except ValueError as error:  # polymeric formula such as C20H37O4P(C5H8)n cannot be parsed
                row["mass_balance_as_built"][rid] = f"not checkable: {error}"
    result[label] = row
    print(label, row["sha256"][:8], row["wt"]["status"], row["wt"]["growth"], row["wt"]["passed"],
          [(r["objective"], r["status"], r["value"]) for r in row["no_uptake"]],
          {k: v["ratio"] for k, v in row["knockouts"].items()}, file=sys.stderr, flush=True)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1, default=str)
    handle.write("\n")
