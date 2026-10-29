"""Hidden-source check on 9d655c2b: can the model grow, or make acyl-CoA / phosphate / inositol, with the source removed?

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/lipid_source_leak_check_20261028/scripts/source_check.py MODEL.xml MEDIA_DIR PROFILE.json OUT.json
Read-only (TASK.md section 4). Gurobi 1e-9 through checked_solve.py (reset per solve, bound / residual /
knocked-flux checks). Demand reactions are added in memory inside `with model:` and never written out.
"""
import hashlib
import json
import math
import sys
import warnings
from pathlib import Path

from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis, pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

model_path, media_dir, profile, out = sys.argv[1:5]
LP_CAP = 300
NONZERO = 1e-7
MEDIA = {"C1": "po1f_csm_leu_w29uptake_20261012.csv", "C2": "po1f_csm_leu_w29uptake_20261012.csv",
         "C3": "po1f_csm_leu_w29glycerol_20261012.csv", "C4": "po1f_csm_leu_w29oleate_20261012.csv",
         "C5": "po1f_csm_leu_w29uptake_20261012.csv"}
FORMULA_LESS = ["R947", "R1067", "R1189", "R1294", "R1984", "R1986", "R1990", "R1992", "R1998", "R2018", "R2026",
                "R2031", "R2036", "R2214", "R2275"]
SOURCES = {
    "acyl": dict(genes=["YALI1B19844g"],
                 exchanges=["R951", "R1083", "R1236", "R1281", "R1288", "R1289", "R1329", "R1440", "R1444"],
                 demands=["m567[C_em]", "m570[C_lp]", "m790[C_mm]", "m1984[C_er]", "m1349[C_cy]"]),
    "phosphate": dict(genes=[],
                      exchanges=["R950", "R951", "R1005", "R1014", "R1102", "R1126", "R1269", "R1296", "R1322",
                                 "R1339", "R1341"],
                      demands=["m35[C_cy]"]),
    "inositol": dict(genes=["YALI1B05798g"], exchanges=["R950", "R1263"], demands=["m814[C_cy]", "m813[C_cy]"]),
}
TRACE_CLOSURES = ["R2171", "R2165", "R1726", "R1842", "R2168"]
LP = {"n": 0}


def count(n=1):
    LP["n"] += n
    if LP["n"] > LP_CAP:
        raise RuntimeError(f"LP cap {LP_CAP} reached")


def context(medium):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=str(Path(media_dir) / medium),
                                                strain_profile_path=profile)
    configure(ctx.model)
    return ctx.model


def set_condition(model, cond, source=None):
    """Apply condition-specific bounds inside an open `with model:` block."""
    if cond == "C2":
        for rid, name in (("R1041", "choline"), ("R1115", "ethanolamine")):
            r = model.reactions.get_by_id(rid)
            assert r.boundary and name in r.name, (rid, r.name)
            r.lower_bound = -1.0
    if cond == "C5":
        for r in model.exchanges:
            r.bounds = (-1000.0, 1000.0)
        for rid in FORMULA_LESS:
            model.reactions.get_by_id(rid).lower_bound = 0.0


def remove_source(model, cond, source):
    """Knock out the source genes and close the source exchanges to uptake; return the closed reaction ids."""
    spec = SOURCES[source]
    closed = set()
    for gene in spec["genes"]:
        g = model.genes.get_by_id(gene)
        g.knock_out()
        closed |= {r.id for r in g.reactions if r.bounds == (0, 0)}
    if not (source == "acyl" and cond == "C4"):  # C4 keeps oleate: the positive control for the acyl test
        for rid in spec["exchanges"]:
            model.reactions.get_by_id(rid).lower_bound = max(0.0, model.reactions.get_by_id(rid).lower_bound)
    return closed


def checked(model, knocked):
    count()
    res = solve(model, matrix(model), knocked)
    if res["growth"] is not None:
        res["R2171"] = model.solver.primal_values[model.reactions.R2171.forward_variable.name] - \
            model.solver.primal_values[model.reactions.R2171.reverse_variable.name]
    return res


def demand_max(model, mid, knocked):
    with model:
        dm = Reaction(f"TESTDM_{mid}", lower_bound=0.0, upper_bound=1000.0)
        dm.add_metabolites({model.metabolites.get_by_id(mid): -1.0})
        model.add_reactions([dm])
        model.objective = dm.id
        return checked(model, knocked)


def trace(model, mid, knocked, value):
    """pFBA route at the demand maximum, then the maximum with each candidate reaction closed."""
    with model:
        dm = Reaction(f"TESTDM_{mid}", lower_bound=0.0, upper_bound=1000.0)
        dm.add_metabolites({model.metabolites.get_by_id(mid): -1.0})
        model.add_reactions([dm])
        model.objective = dm.id
        model.solver.problem.reset()
        count(2)
        sol = pfba(model, fraction_of_optimum=0.999)
        route = {rid: float(v) for rid, v in sol.fluxes.items() if abs(v) > 1e-9}
        closures = {}
        for rid in TRACE_CLOSURES:
            with model:
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                res = checked(model, knocked | {rid})
                closures[rid] = res["growth"]
    return {"demand_max": value, "route_n": len(route), "route": route, "closures": closures}


def run_block(model, cond, variant, result):
    block = {"wt": None, "sources": {}}
    with model:
        set_condition(model, cond)
        block["wt"] = checked(model, set())
        for source, spec in SOURCES.items():
            with model:
                closed = remove_source(model, cond, source)
                row = {"knocked_reactions": sorted(closed), "growth": checked(model, closed), "demands": {}}
                for mid in spec["demands"]:
                    row["demands"][mid] = demand_max(model, mid, closed)
                block["sources"][source] = row
    result["blocks"][f"{cond}|{variant}"] = block
    print(cond, variant, "WT", block["wt"]["growth"], {s: (r["growth"]["growth"], {m: d["growth"] for m, d in r["demands"].items()})
                                                       for s, r in block["sources"].items()}, "LPs", LP["n"], file=sys.stderr, flush=True)


result = {"model": model_path, "sha256": hashlib.sha256(Path(model_path).read_bytes()).hexdigest(),
          "media": {c: m for c, m in MEDIA.items()}, "lp_cap": LP_CAP, "nonzero_threshold": NONZERO,
          "blocks": {}, "r2171_fva": {}, "traces": {}}
models = {}
try:
    for cond, medium in MEDIA.items():
        if medium not in models:
            models[medium] = context(medium)
        model = models[medium]
        run_block(model, cond, "base", result)
        with model:
            set_condition(model, cond)
            fva = {}
            for frac in (0.999, 0.0):
                count(3)
                model.solver.problem.reset()
                f = flux_variability_analysis(model, reaction_list=["R2171"], fraction_of_optimum=frac, processes=1)
                fva[str(frac)] = [float(f.loc["R2171", "minimum"]), float(f.loc["R2171", "maximum"])]
            result["r2171_fva"][cond] = fva
        if cond in ("C1", "C5"):
            with model:
                model.reactions.R1720.bounds = (0.0, 1000.0)
                run_block(model, cond, "R1720_open", result)
    # traces of every non-zero source-removed demand (not the C4 acyl positive control)
    for key, block in result["blocks"].items():
        cond, variant = key.split("|")
        for source, row in block["sources"].items():
            if source == "acyl" and cond == "C4":
                continue
            for mid, res in row["demands"].items():
                v = res["growth"]
                if v is not None and v > NONZERO:
                    model = models[MEDIA[cond]]
                    with model:
                        set_condition(model, cond)
                        if variant == "R1720_open":
                            model.reactions.R1720.bounds = (0.0, 1000.0)
                        closed = remove_source(model, cond, source)
                        result["traces"][f"{key}|{source}|{mid}"] = trace(model, mid, closed, v)
                    print("traced", key, source, mid, v, "LPs", LP["n"], file=sys.stderr, flush=True)
except RuntimeError as exc:
    result["stopped"] = str(exc)
    print("STOPPED:", exc, file=sys.stderr, flush=True)
result["lp_count"] = LP["n"]
with open(out, "w") as handle:
    json.dump(result, handle, indent=1, default=lambda o: None if isinstance(o, float) and math.isnan(o) else str(o))
    handle.write("\n")
