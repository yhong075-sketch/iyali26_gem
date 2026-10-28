"""Stage 1 of the UM7 nuclear-transport task: in-memory LP checks (read-only; writes only OUT.json).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts \
    /opt/miniconda3/bin/python stage1_um7.py MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.json
Gurobi at 1e-9 (checked_solve.py), real HOME. Conditions: 'open' = raw model, every exchange (-1000, 1000)
(the census's "any medium"); 'taskA' = Task A medium + PO1f profile. Variants are applied in memory only.
Per variant: WT growth on Task A (checked); FVA at 0 % of the optimum of the UM7 reactions and context
reactions under both conditions, and at 99.9 % on Task A; pFBA flux on Task A; energy-generating-cycle LPs
(every exchange (0, 1000), positive lower bounds relaxed incl. xMAINTENANCE; max R_NTP1 and max biomass;
status reported); targeted single-gene knockouts on Task A with the stale-point check.
"""
import json, math, sys, warnings
from pathlib import Path
import cobra
from cobra import Metabolite, Reaction
from cobra.flux_analysis import flux_variability_analysis, pfba
from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.dfba_new_fn_essentiality import load_experimental_reference

UM7 = ["R92", "R396", "R567", "R568", "R569", "R579", "R580", "R965", "R1000", "R1001", "R1017", "R1055",
       "R1082", "R1212", "R1214", "R1267", "R1270"]
CONTEXT = ["R1912", "R1170", "R24", "R578", "R90", "R91", "R2024", "R2268"]
EXTRA_GENES = ["YALI1F33579g"]
model_path, medium, profile, experimental, out = sys.argv[1:6]
positives = set(load_experimental_reference(Path(experimental))[0])
lp = {"n": 0}

# test reactions (gene-less unless stated); ids prefixed T_ so they cannot collide with model ids
LINKS = {
    "T_ADPRtn": ("ADP-ribose nucleus<->cytosol (pore diffusion)", {"m818[C_nu]": -1, "m1725[C_cy]": 1}),
    "T_PPItn": ("diphosphate nucleus<->cytosol (pore diffusion)", {"m210[C_nu]": -1, "m203[C_cy]": 1}),
    "T_NADtn": ("NAD nucleus<->cytosol (pore diffusion)", {"m817[C_nu]": -1, "m122[C_cy]": 1}),
}


def add(model, rid, name, mets, bounds=(-1000, 1000)):
    r = Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in mets.items()})
    return r


def homocitrate_cy(model):
    if "m_T_hcit[C_cy]" not in model.metabolites:
        model.add_metabolites([Metabolite("m_T_hcit[C_cy]", formula="C7H10O7", charge=0, name="homocitrate (test)",
                                          compartment="C_cy")])
    add(model, "T_HCITtn", "homocitrate nucleus<->cytosol (pore diffusion)", {"m628[C_nu]": -1, "m_T_hcit[C_cy]": 1})
    add(model, "T_HCITtm", "homocitrate cytosol<->mitochondrion (carrier unknown)", {"m_T_hcit[C_cy]": -1, "m75[C_mi]": 1})


def apply(model, label):
    parts = label.split("+")
    for p in parts:
        if p == "as_built":
            continue
        elif p in ("adpr", "ppi", "nad"):
            rid = {"adpr": "T_ADPRtn", "ppi": "T_PPItn", "nad": "T_NADtn"}[p]
            add(model, rid, *LINKS[rid])
        elif p == "close396":
            model.reactions.R396.bounds = (0, 0)
        elif p == "hcit":
            homocitrate_cy(model)
        elif p == "open92":
            model.reactions.R92.bounds = (0, 1000)
        elif p == "censusH11":  # census targeted test H11: homocitrate nu<->mi + sinks on PPi[nu], ADP-ribose[nu]
            add(model, "T_HCIT_nm", "homocitrate nu<->mi (census H11)", {"m628[C_nu]": -1, "m75[C_mi]": 1})
            add(model, "T_SK_ppi_nu", "sink PPi[nu]", {"m210[C_nu]": -1})
            add(model, "T_SK_adpr_nu", "sink ADP-ribose[nu]", {"m818[C_nu]": -1})
        elif p == "collapse":  # R567 rewritten onto the cytosolic species; nuclear R567/R569/R579 closed
            add(model, "T_R567c", "sirtuin NAD hydrolysis, cytosolic copy (test)",
                {"m32[C_cy]": -1, "m122[C_cy]": -1, "m10[C_cy]": 2, "m1725[C_cy]": 1, "m829[C_cy]": 1}, (0, 1000))
            for rid in ("R567", "R569", "R579"):
                model.reactions.get_by_id(rid).bounds = (0, 0)
        else:
            raise ValueError(p)


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
    return ctx.model


VARIANTS = ["as_built", "adpr", "ppi", "nad", "adpr+ppi", "adpr+ppi+nad", "adpr+ppi+nad+close396",
            "adpr+ppi+nad+hcit", "adpr+ppi+nad+hcit+open92", "censusH11", "collapse", "collapse+ppi"]
if len(sys.argv) > 6:
    VARIANTS = sys.argv[6].split(",")
result = {"inputs": {"model": model_path, "medium": medium, "profile": profile, "experimental": experimental},
          "variants": {}}
OPEN, TASKA = open_model(), taska_model()
genes = sorted({g.id for rid in UM7 + CONTEXT for g in TASKA.reactions.get_by_id(rid).genes} | set(EXTRA_GENES))
result["genes"] = {g: sorted(r.id for r in TASKA.genes.get_by_id(g).reactions) for g in genes}
for label in VARIANTS:
    row = {}
    with OPEN as mo:
        apply(mo, label)
        ids = UM7 + CONTEXT + [r.id for r in mo.reactions if r.id.startswith("T_")]
        mo.solver.problem.reset(); lp["n"] += 1
        row["open_growth"] = mo.slim_optimize(error_value=float("nan"))
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
    blocked = lambda f: sorted(r for r in UM7 if max(abs(f[r][0]), abs(f[r][1])) <= 1e-9)
    row["um7_blocked_open"] = blocked(row["fva_open_0"])
    row["um7_blocked_taskA"] = blocked(row["fva_taskA_0"])
    result["variants"][label] = row
    print(f"{label:28s} WT {wt['status']} {wt['growth']!r} passed={wt['passed']} | open growth {row['open_growth']:.6g} | "
          f"UM7 blocked open {len(row['um7_blocked_open'])}/17 taskA {len(row['um7_blocked_taskA'])}/17 | "
          f"EGC ntp1 {row['egc']['R_NTP1']['status']} {row['egc']['R_NTP1']['max_no_uptake']} "
          f"bio {row['egc']['biomass_C']['status']} {row['egc']['biomass_C']['max_no_uptake']}", flush=True)
    print("     blocked open:", row["um7_blocked_open"], flush=True)
    for g, k in kos.items():
        tag = "LETHAL<15%" if k["ratio"] is not None and k["ratio"] < 0.15 else ""
        if tag or not k["passed"] or (k["ratio"] is not None and k["ratio"] < 0.999999):
            print(f"     KO {g} ratio={k['ratio']} status={k['status']} passed={k['passed']} {tag} "
                  f"{'consensus' if k['on_consensus_list'] else ''}", flush=True)
result["lp_count"] = lp["n"]
json.dump(result, open(out, "w"), indent=1)
print("LPs", lp["n"])
