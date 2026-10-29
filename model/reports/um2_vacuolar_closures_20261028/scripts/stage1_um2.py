"""Stage 1 of the UM2 vacuolar-closure task: in-memory LP checks per option (read-only; writes only OUT.json).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python stage1_um2.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json [VARIANTS]
Gurobi at 1e-9 (checked_solve.py), real HOME. Conditions: 'open' = raw model, every exchange (-1000, 1000)
(the census's "any medium"); 'taskA' = Task A medium + PO1f profile. Variants are applied in memory only.
Per variant: WT growth on Task A (checked); FVA at 0 % of the optimum of the 37 UM2 reactions and context
reactions under both conditions, and at 99.9 % on Task A; pFBA flux on Task A; energy-generating-cycle LPs
(every exchange (0, 1000), positive lower bounds relaxed incl. xMAINTENANCE; max R_NTP1 and max biomass;
status reported); a no-uptake FVA of the same reactions under those relaxed bounds (any non-zero range is a
flux mode with no substrate, i.e. an internal cycle).
"""
import json, math, sys, warnings
from pathlib import Path
import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

UM2 = ("R795 R869 R870 R871 R876 R878 R879 R882 R883 R887 R888 R897 R898 R918 R925 R927 R1163 R1363 R2018 "
       "R2019 R2020 R2021 R2026 R2027 R2028 R2029 R2030 R2031 R2032 R2033 R2034 R2035 R2036 R2037 R2038 R2039 "
       "R2040").split()
CONTEXT = ["R794", "R1159", "R2055", "R2056", "R1655", "R_NTP1", "xMAINTENANCE"]
GROUPS = {  # option building blocks; bounds as (lb, ub)
    "pump": {"R795": (0, 1000)},
    "water": {"R1363": (0, 1000)},
    "hdiff": {"R1163": (0, 1000)},        # as written: H+ cyt -> vac
    "hleak": {"R1163": (-1000, 0)},       # reversed: passive leak vac -> cyt only
    "aaimp": {r: (0, 1000) for r in ("R869", "R878", "R887", "R897")},          # H+ antiport into vacuole
    "aaexp": {r: (0, 1000) for r in ("R870", "R879", "R883", "R888", "R898")},  # H+ symport out of vacuole
    "aspglu": {r: (0, 1000) for r in ("R871", "R876")},
    "pa": {r: (0, 1000) for r in ("R918", "R925", "R927")},
}
model_path, medium, profile, out = sys.argv[1:5]
VARIANTS = sys.argv[5].split(",") if len(sys.argv) > 5 else [
    "as_built", "pump", "water", "pump+water", "pump+water+hleak", "pump+water+hdiff",
    "pump+water+aaexp", "pump+water+aaimp+aaexp", "pump+water+aspglu", "pump+water+pa",
    "pump+water+aaimp+aaexp+aspglu+pa", "pump+water+hleak+aaimp+aaexp+aspglu+pa",
    "pump+water+hdiff+aaimp+aaexp+aspglu+pa"]
lp = {"n": 0}


def apply(model, label):
    for part in label.split("+"):
        if part == "as_built":
            continue
        for rid, b in GROUPS[part].items():
            model.reactions.get_by_id(rid).bounds = b


def fva(model, ids, frac):
    lp["n"] += 2 * len(ids)
    df = flux_variability_analysis(model, reaction_list=ids, fraction_of_optimum=frac, processes=1)
    return {r: [float(df.loc[r, "minimum"]), float(df.loc[r, "maximum"])] for r in ids}


def relaxed_no_uptake(model):
    relaxed = []
    for r in model.exchanges:
        r.bounds = (0.0, 1000.0)
    for r in model.reactions:
        if r.lower_bound > 0:
            relaxed.append([r.id, r.lower_bound]); r.lower_bound = 0.0
    return relaxed


def egc(model, ids):
    out = {}
    with model:
        out["relaxed"] = relaxed_no_uptake(model)
        for obj in ("R_NTP1", "biomass_C"):
            model.objective = obj
            model.solver.problem.reset(); lp["n"] += 1
            v = model.slim_optimize(error_value=float("nan"))
            out[obj] = {"status": model.solver.status, "max_no_uptake": None if math.isnan(v) else v}
        model.objective = "biomass_C"
        model.reactions.biomass_C.bounds = (0, 0)
        model.solver.problem.reset()
        f = fva(model, ids, 0.0)
        out["no_uptake_fva_nonzero"] = {r: v for r, v in f.items() if max(abs(v[0]), abs(v[1])) > 1e-9}
    return out


def open_model():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = cobra.io.read_sbml_model(Path(model_path))
    configure(m)
    for r in m.exchanges:
        r.bounds = (-1000, 1000)
    return m


def taska_model():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=Path(model_path), media_path=Path(medium),
                                                strain_profile_path=Path(profile))
    configure(ctx.model)
    return ctx


OPEN, CTX = open_model(), taska_model()
TASKA = CTX.model
ids = UM2 + CONTEXT
result = {"inputs": {"model": model_path, "model_sha256": CTX.canonical_model_sha256, "medium": medium,
                     "medium_sha256": CTX.medium_sha256, "profile": profile,
                     "profile_sha256": CTX.strain_profile_sha256}, "groups": {k: {r: list(b) for r, b in v.items()}
                                                                               for k, v in GROUPS.items()},
          "variants": {}}
blocked = lambda f: sorted(r for r in UM2 if max(abs(f[r][0]), abs(f[r][1])) <= 1e-9)
for label in VARIANTS:
    row = {}
    with OPEN as mo:
        apply(mo, label)
        mo.solver.problem.reset(); lp["n"] += 1
        row["open_growth"] = mo.slim_optimize(error_value=float("nan"))
        row["fva_open_0"] = fva(mo, ids, 0.0)
        row["egc_open_model"] = egc(mo, ids)
    with TASKA as m:
        apply(m, label)
        lp["n"] += 1
        wt = solve(m, matrix(m), set())
        row["wt"] = wt
        row["fva_taskA_0"] = fva(m, ids, 0.0)
        row["fva_taskA_999"] = fva(m, ids, 0.999)
        m.solver.problem.reset(); lp["n"] += 1
        sol = pfba(m)
        row["pfba_taskA"] = {r: float(sol.fluxes[r]) for r in ids if abs(sol.fluxes[r]) > 1e-12}
        row["egc_taskA_model"] = egc(m, ids)
    row["um2_blocked_open"] = blocked(row["fva_open_0"])
    row["um2_blocked_taskA"] = blocked(row["fva_taskA_0"])
    result["variants"][label] = row
    e = row["egc_taskA_model"]
    print(f"{label:42s} WT {wt['status']} {wt['growth']!r} passed={wt['passed']} | open growth {row['open_growth']:.10g} | "
          f"UM2 blocked open {len(row['um2_blocked_open'])}/37 taskA {len(row['um2_blocked_taskA'])}/37 | "
          f"EGC ntp1 {e['R_NTP1']['status']} {e['R_NTP1']['max_no_uptake']} bio {e['biomass_C']['status']} "
          f"{e['biomass_C']['max_no_uptake']} | no-uptake nonzero {sorted(e['no_uptake_fva_nonzero'])} "
          f"/open-model {sorted(row['egc_open_model']['no_uptake_fva_nonzero'])}", flush=True)
    print("     blocked open:", row["um2_blocked_open"], flush=True)
    print("     free taskA:", sorted(set(UM2) - set(row["um2_blocked_taskA"])), flush=True)
result["lp_count"] = lp["n"]
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp["n"])
