"""Stage 1 of the UM11 nuclear-RNR task: in-memory LP checks (read-only; writes only OUT.json).

Adapted from model/reports/um7_nuclear_transports_20261028/scripts/stage1_um7.py (same conditions and checks;
only the reaction lists and variants differ).
Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts:SCRIPTDIR \
    /opt/miniconda3/bin/python stage1_um11.py MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.json [VARIANTS]
Gurobi at 1e-9 (checked_solve.py), real HOME. Conditions: 'open' = raw model, every exchange (-1000, 1000);
'taskA' = Task A medium + PO1f profile. Per variant: WT growth on Task A (checked); FVA at 0 % of the optimum
of the UM11 and context reactions under both conditions and at 99.9 % on Task A; pFBA flux on Task A;
energy-generating-cycle LPs (every exchange (0, 1000), positive lower bounds relaxed incl. xMAINTENANCE; max
R_NTP1 and max biomass; status reported); targeted single-gene knockouts on Task A with the stale-point check.
"""
import json, math, sys, warnings
from pathlib import Path
import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.dfba_new_fn_essentiality import load_experimental_reference
from um11_apply import apply

UM11 = ["R707", "R709", "R711", "R1914", "R782", "R1033", "R1076", "R1077", "R1093", "R1136", "R1353"]
CONTEXT = ["R706", "R708", "R710", "R1913", "R702", "R703", "R704", "R705", "R755", "R781", "R789", "R1104",
           "R1313", "R1286", "R1008", "R1021", "R1361"]
model_path, medium, profile, experimental, out = sys.argv[1:6]
positives = set(load_experimental_reference(Path(experimental))[0])
lp = {"n": 0}


def fva(model, ids, frac):
    lp["n"] += 2 * len(ids)
    df = flux_variability_analysis(model, reaction_list=ids, fraction_of_optimum=frac, processes=1)
    return {r: [float(df.loc[r, "minimum"]), float(df.loc[r, "maximum"])] for r in ids}


def egc(model):
    out = {}
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = []
        for r in model.reactions:
            if r.lower_bound > 0:
                relaxed.append([r.id, r.lower_bound]); r.lower_bound = 0.0
        for obj in ("R_NTP1", "biomass_C"):
            model.objective = obj
            model.solver.problem.reset(); lp["n"] += 1
            v = model.slim_optimize(error_value=float("nan"))
            out[obj] = {"status": model.solver.status, "max_no_uptake": None if math.isnan(v) else v}
        out["relaxed"] = relaxed
    return out


def open_model():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = cobra.io.read_sbml_model(model_path)
    configure(m)
    for r in m.exchanges:
        r.bounds = (-1000, 1000)
    return m


def taska_model():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


VARIANTS = ["as_built", "trx", "trxss", "trx+trxss", "trx+trxss+h2o1914", "ntrr", "close_nu"]
if len(sys.argv) > 6:
    VARIANTS = sys.argv[6].split(",")
OPEN, CTX = open_model(), taska_model()
TASKA = CTX.model
result = {"inputs": {"model": model_path, "model_sha256": CTX.canonical_model_sha256, "medium": medium,
                     "medium_sha256": CTX.medium_sha256, "profile": profile,
                     "profile_sha256": CTX.strain_profile_sha256, "experimental": experimental}, "variants": {}}
genes = sorted({g.id for rid in UM11 + CONTEXT for g in TASKA.reactions.get_by_id(rid).genes})
result["genes"] = {g: sorted(r.id for r in TASKA.genes.get_by_id(g).reactions) for g in genes}
for label in VARIANTS:
    row = {}
    with OPEN as mo:
        apply(mo, label)
        ids = UM11 + CONTEXT + [r.id for r in mo.reactions if r.id.startswith("T_")]
        mo.solver.problem.reset(); lp["n"] += 1
        row["open_growth"] = mo.slim_optimize(error_value=float("nan"))
        row["open_status"] = mo.solver.status
        row["fva_open_0"] = fva(mo, ids, 0.0)
    with TASKA as m:
        apply(m, label)
        stoich = matrix(m)
        lp["n"] += 1
        wt = solve(m, stoich, set())
        row["wt"] = wt
        row["fva_taskA_0"] = fva(m, ids, 0.0)
        row["fva_taskA_999"] = fva(m, ids, 0.999)
        m.solver.problem.reset(); lp["n"] += 1
        sol = pfba(m)
        row["pfba_taskA"] = {r: float(sol.fluxes[r]) for r in ids if abs(sol.fluxes[r]) > 1e-12}
        row["egc"] = egc(m)
        kos = {}
        for g in genes:
            with m:
                gene = m.genes.get_by_id(g)
                before = {r.id: r.bounds for r in gene.reactions}
                gene.knock_out()
                knocked = {rid for rid, b in before.items() if m.reactions.get_by_id(rid).bounds == (0, 0) and b != (0, 0)}
                lp["n"] += 1
                k = solve(m, stoich, knocked)
            k["ratio"] = None if k["growth"] is None or not wt["growth"] else k["growth"] / wt["growth"]
            k["knocked"] = sorted(knocked)
            k["on_consensus_list"] = g in positives
            kos[g] = k
        row["knockouts"] = kos
    blocked = lambda f: sorted(r for r in UM11 if max(abs(f[r][0]), abs(f[r][1])) <= 1e-9)
    row["um11_blocked_open"] = blocked(row["fva_open_0"])
    row["um11_blocked_taskA"] = blocked(row["fva_taskA_0"])
    result["variants"][label] = row
    print(f"{label:22s} WT {wt['status']} {wt['growth']!r} passed={wt['passed']} | open {row['open_status']} "
          f"{row['open_growth']:.6g} | UM11 blocked open {len(row['um11_blocked_open'])}/11 taskA "
          f"{len(row['um11_blocked_taskA'])}/11 | EGC ntp1 {row['egc']['R_NTP1']['status']} "
          f"{row['egc']['R_NTP1']['max_no_uptake']} bio {row['egc']['biomass_C']['status']} "
          f"{row['egc']['biomass_C']['max_no_uptake']}", flush=True)
    print("     blocked open:", row["um11_blocked_open"], "| taskA:", row["um11_blocked_taskA"], flush=True)
    print("     pFBA nonzero:", {k: round(v, 9) for k, v in row["pfba_taskA"].items()}, flush=True)
    for g, k in kos.items():
        tag = "LETHAL<15%" if k["ratio"] is not None and k["ratio"] < 0.15 else ""
        if tag or not k["passed"] or (k["ratio"] is not None and k["ratio"] < 0.999999):
            print(f"     KO {g} ratio={k['ratio']} growth={k['growth']} status={k['status']} passed={k['passed']} {tag} "
                  f"{'consensus' if k['on_consensus_list'] else ''}", flush=True)
result["lp_count"] = lp["n"]
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp["n"])
