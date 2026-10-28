"""Pre-check for re-running the R1726 opening on the R2135 candidate e9d35857 (read-only, about 25 LPs).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/kennedy_closures_20261027/scripts/rerun_precheck.py CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Gurobi 1e-9. (1) Open condition (every exchange (-1000, 1000), as in the census), R1726 closed vs opened in
memory: FVA (fraction 0, processes=1) of the UM12 reactions. (2) Task A + PO1f with choline R1041 and
ethanolamine R1115 uptake 1, R1726 closed vs opened: WT and knockouts of the PE-methylation genes
YALI1E07347g and YALI1E15432g (checked solves).
"""
import json
import sys
import warnings

from cobra.flux_analysis import flux_variability_analysis

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

UM12 = ["R1726", "R212", "R213", "R1603", "R830", "R1041", "R1042", "R1743", "R1842", "R1809", "R354"]
path, medium, profile, out = sys.argv[1:5]


def load():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = load_effective_simulation_context(model_path=path, media_path=medium, strain_profile_path=profile).model
    configure(model)
    return model


result = {}
for opened in (False, True):
    tag = "r1726_open" if opened else "as_built"
    model = load()
    if opened:
        model.reactions.R1726.bounds = (0.0, 1000.0)
    with model:
        for r in model.exchanges:
            r.bounds = (-1000.0, 1000.0)
        frame = flux_variability_analysis(model, reaction_list=UM12, fraction_of_optimum=0.0, processes=1)
        fva = {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in UM12}
    for mid in ("m1037[C_ex]", "m1047[C_ex]"):
        [r for r in model.metabolites.get_by_id(mid).reactions if r.boundary][0].lower_bound = -1.0
    s = matrix(model)
    control = {"wt": solve(model, s, set())}
    for gene in ("YALI1E07347g", "YALI1E15432g"):
        with model:
            model.genes.get_by_id(gene).knock_out()
            knocked = {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}
            control[gene] = solve(model, s, knocked)
    result[tag] = {"open_condition_fva": fva, "choline_ethanolamine_control": control}
    print(tag, {k: [round(a, 6), round(b, 6)] for k, (a, b) in fva.items()},
          {k: (v["status"], None if v["growth"] is None else round(v["growth"], 7)) for k, v in control.items()}, file=sys.stderr)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
