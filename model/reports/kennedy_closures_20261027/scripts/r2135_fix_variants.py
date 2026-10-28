"""In-memory test of candidate fixes for the R2135 acetyl-CoA -> acyl-CoA leak and the 0.01-scaled R2145 (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/kennedy_closures_20261027/scripts/r2135_fix_variants.py BASELINE.xml R1726_CANDIDATE.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.json
Task A medium + PO1f, Gurobi 1e-9, checked solves. Variants, applied to both models:
  as_built; remap (R2135 = m567[C_em] <=> m1984[C_er], the ER-membrane pooled acyl-CoA instead of acetyl-CoA);
  close (R2135 (0, 0)); rescale (R2145 = m1993 + H2O -> m1994 + Pi, all coefficients 1); remap+rescale;
  close+rescale.
Per variant: WT (checked); phosphate exchange flux at a pFBA point; FVA at 99.9 % of R2135, R2168, R2145, R2134;
the no-uptake LPs (every exchange (0, 1000), positive lower bounds relaxed; max R_NTP1 and max biomass_C);
targeted knockouts of the lipid genes listed below. Diagnostic only: nothing is written to a model file.
"""
import json
import math
import sys
import warnings
from pathlib import Path

from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.dfba_new_fn_essentiality import load_experimental_reference

GENES = ("YALI1B19844g", "YALI1B25427g", "YALI1C07638g", "YALI1D22124g",  # changed class with R1726 open
         "YALI1D10264g", "YALI1E38810g", "YALI1E20049g", "YALI1D35593g",  # R2168/R2165, R1771, PDAT, R2145
         "YALI1F09542g", "YALI1D05038g", "YALI1D08000g", "YALI1C15420g", "YALI1E22633g",
         "YALI1D10930g", "YALI1E17431g", "YALI1E07347g", "YALI1E15432g", "YALI1F27002g")
FVA = ["R2135", "R2168", "R2145", "R2134"]
baseline, candidate, medium, profile, experimental, out = sys.argv[1:7]
positives = set(load_experimental_reference(Path(experimental))[0])


def remap(model):
    r = model.reactions.R2135
    r.subtract_metabolites({model.metabolites.get_by_id("m68[C_cy]"): -1.0})
    r.add_metabolites({model.metabolites.get_by_id("m567[C_em]"): -1.0})


def close(model):
    model.reactions.R2135.bounds = (0.0, 0.0)


def rescale(model):
    r = model.reactions.R2145
    r.add_metabolites({model.metabolites.get_by_id("m1993[C_cy]"): -0.99, model.metabolites.get_by_id("m1994[C_cy]"): 0.99})


VARIANTS = {"as_built": (), "remap": (remap,), "close": (close,), "rescale": (rescale,),
            "remap+rescale": (remap, rescale), "close+rescale": (close, rescale)}


def no_uptake(model, objective):
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        for r in model.reactions:
            if r.lower_bound > 0:
                r.lower_bound = 0.0
        model.objective = objective
        model.solver.problem.reset()
        value = model.slim_optimize(error_value=float("nan"))
        return [model.solver.status, None if math.isnan(value) else value]


result = {"genes": {g: g in positives for g in GENES}, "models": {}}
for label, path in (("baseline", baseline), ("r1726_candidate", candidate)):
    rows = {}
    for name, edits in VARIANTS.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile).model
        configure(model)
        for edit in edits:
            edit(model)
        s = matrix(model)
        wt = solve(model, s, set())
        row = {"equations": {rid: model.reactions.get_by_id(rid).build_reaction_string() for rid in ("R2135", "R2145")},
               "wt": wt}
        if wt["growth"]:
            fluxes = pfba(model).fluxes
            row["pfba"] = {rid: float(fluxes[rid]) for rid in FVA + ["R1296", "R1771", "R1726"]}
            frame = flux_variability_analysis(model, reaction_list=FVA, fraction_of_optimum=0.999, processes=1)
            row["fva_0.999"] = {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in FVA}
        row["no_uptake"] = {"R_NTP1": no_uptake(model, "R_NTP1"), "biomass_C": no_uptake(model, "biomass_C")}
        kos = {}
        for gene in GENES:
            with model:
                model.genes.get_by_id(gene).knock_out()
                knocked = {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}
                point = solve(model, s, knocked)
            kos[gene] = None if point["growth"] is None or not wt["growth"] else point["growth"] / wt["growth"]
        row["ko_ratio"] = kos
        rows[name] = row
        print(label, name, wt["status"], wt["growth"], row["no_uptake"],
              {g: (None if v is None else round(v, 4)) for g, v in kos.items()}, file=sys.stderr, flush=True)
    result["models"][label] = rows
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
