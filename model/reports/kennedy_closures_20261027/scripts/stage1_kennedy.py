"""Stage 1 of the Kennedy-closure task: in-memory LP checks of opening R1725 / R1726 (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    stage1_kennedy.py MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT.json
Task A medium + PO1f, Gurobi at 1e-9 (checked_solve.py: reset before each solve, targeted knocked-flux check).
Variants: closed (as built), R1725 open, R1726 open, both open; opened bounds (0, 1000) (both are declared
irreversible). Per variant: WT growth; pFBA fluxes of the PC / PE routes (both open); FVA of R1725 / R1726 at
100 % and 99.9 % of the optimum and at growth >= 0; the energy-generating-cycle LP (every exchange (0, 1000),
every positive lower bound relaxed to 0, max R_NTP1, status recorded); targeted knockouts of the PC / PE route
genes. Mechanism control (not a benchmark condition): choline and ethanolamine exchanges opened to uptake 1.
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

TARGETS = ("R1725", "R1726")
VARIANTS = {"closed": (), "R1725_open": ("R1725",), "R1726_open": ("R1726",), "both_open": TARGETS}
ROUTE = ["R1725", "R1726", "R212", "R213", "R1603", "R293", "R2080", "R635", "R1602", "SPHPL", "R730", "R663",
         "R1841", "R196", "R197", "R629", "R630", "R631", "R1835", "R1845", "R1844", "R2187", "R1743", "R1842",
         "R354", "R1730", "R1731", "R1729", "R1728", "R1839", "R1787", "R1786", "R1706", "R1705", "R1793", "R1799",
         "R1803", "R_PEt_em_er", "R1806", "R1807"]
CHOLINE_ETN = ("m1037[C_ex]", "m1047[C_ex]")
model_path, medium, profile, experimental, out = sys.argv[1:6]
positives = set(load_experimental_reference(Path(experimental))[0])
lp_count = 0


def context():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx


def open_targets(model, ids):
    for rid in ids:
        model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)


def checked(model, stoichiometry, knocked=()):
    global lp_count
    lp_count += 1
    return solve(model, stoichiometry, set(knocked))


def exchange_of(model, mid):
    met = model.metabolites.get_by_id(mid)
    found = [r for r in met.reactions if r.boundary]
    assert len(found) == 1, (mid, [r.id for r in found])
    return found[0]


def egc(model):
    global lp_count
    with model:
        for r in model.exchanges:
            r.bounds = (0.0, 1000.0)
        relaxed = []
        for r in model.reactions:
            if r.lower_bound > 0:
                relaxed.append([r.id, r.lower_bound])
                r.lower_bound = 0.0
        model.objective = "R_NTP1"
        model.solver.problem.reset()
        lp_count += 1
        value = model.slim_optimize(error_value=float("nan"))
        return {"status": model.solver.status, "max_R_NTP1_no_uptake": None if math.isnan(value) else value,
                "relaxed_positive_lower_bounds": relaxed}


def fva(model):
    global lp_count
    rows = {}
    for fraction in (1.0, 0.999, 0.0):
        lp_count += 2 * len(TARGETS) + 1
        frame = flux_variability_analysis(model, reaction_list=list(TARGETS), fraction_of_optimum=fraction, processes=1)
        rows[str(fraction)] = {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in TARGETS}
    return rows


def route_genes(model):
    genes = set()
    for rid in ROUTE:
        if rid in model.reactions:
            genes |= {g.id for g in model.reactions.get_by_id(rid).genes}
    return sorted(genes)


def knockouts(model, stoichiometry, genes, wt):
    rows = {}
    for gene in genes:
        with model:
            model.genes.get_by_id(gene).knock_out()
            knocked = {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}
            row = checked(model, stoichiometry, knocked)
        row["ratio"] = None if row["growth"] is None else row["growth"] / wt
        rows[gene] = row
    return rows


result = {"model": model_path, "variants": {}, "mechanism_control_choline_ethanolamine_uptake_1": {}}
ctx = context()
result["model_sha256"] = ctx.canonical_model_sha256
genes = route_genes(ctx.model)
result["route_genes"] = {g: {"on_consensus_list": g in positives,
                             "reactions": sorted(r.id for r in ctx.model.genes.get_by_id(g).reactions)} for g in genes}
for label, ids in VARIANTS.items():
    ctx = context()
    model = ctx.model
    open_targets(model, ids)
    stoichiometry = matrix(model)
    wt = checked(model, stoichiometry)
    row = {"wt": wt, "fva": fva(model), "egc": egc(model)}
    if label in ("closed", "both_open"):
        lp_count += 1
        fluxes = pfba(model).fluxes
        row["pfba_route_fluxes"] = {rid: float(fluxes[rid]) for rid in ROUTE if rid in model.reactions}
        row["knockouts"] = knockouts(model, stoichiometry, genes, wt["growth"])
    result["variants"][label] = row
    print(label, wt["status"], wt["growth"], row["egc"]["status"], row["egc"]["max_R_NTP1_no_uptake"], file=sys.stderr, flush=True)

for label, ids in (("closed", ()), ("both_open", TARGETS)):
    ctx = context()
    model = ctx.model
    open_targets(model, ids)
    for mid in CHOLINE_ETN:
        exchange_of(model, mid).lower_bound = -1.0
    stoichiometry = matrix(model)
    wt = checked(model, stoichiometry)
    lp_count += 1
    fluxes = pfba(model).fluxes
    result["mechanism_control_choline_ethanolamine_uptake_1"][label] = {
        "exchanges": [exchange_of(model, mid).id for mid in CHOLINE_ETN], "wt": wt,
        "pfba_route_fluxes": {rid: float(fluxes[rid]) for rid in ROUTE if rid in model.reactions},
        "knockouts": knockouts(model, stoichiometry, genes, wt["growth"])}

result["lp_count_approx"] = lp_count
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
print("LPs", lp_count, file=sys.stderr)
