"""Mechanism classification of the false-negative genes of the a1237821 screen.

Reads the screen in $IYALI26_RESEARCH_ROOT/artifacts/tasks/essentiality_fn_causes_20261027/screen_a1237821/
(same model, medium, profile, solver settings as run_screen.py) and, for every consensus-essential gene
predicted non-essential at 15 %, records:

- GPR structure: which reactions the gene KO closes (from the screen) and which OR partners keep the
  others open;
- control 1, "bypass closed": all reactions of the gene closed regardless of GPR (reaction-level KO);
- WT flux variability (fraction_of_optimum 1.0) of every reaction the gene is on, and whether the
  reaction is structurally blocked (find_blocked_reactions, medium as set, no growth constraint);
- control 2, uptake shift: pFBA of the reaction-level KO versus WT pFBA, listing medium uptakes that rise
  (amino acids, adenine, uracil, vitamins are capped in the Task A medium) and the reactions whose
  flux changes most (the bypass route, one optimum only).

Investigation only: nothing is written to the model. Outputs go to the task's artifact folder
(fn_classification.tsv, fn_classification.json) and the compact CSV to the report folder.

Usage (repository root, IYALI26_RESEARCH_ROOT set):
    python -B model/reports/essentiality_fn_causes_20261027/classify_fn.py
"""
import csv
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "platform"))

from cobra.flux_analysis import find_blocked_reactions, flux_variability_analysis, pfba
from cobra.util.solver import linear_reaction_coefficients
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL
from tools.diagnose_closed_energy import configure_solver

HERE = Path(__file__).resolve().parent
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0, FeasibilityTol=1e-9, OptimalityTol=1e-9, IntFeasTol=1e-9)
CUTOFF = 0.15
TOL = 1e-6          # growth / flux resolution used for the categories (FVA is unreliable below ~1e-7)
MODEL_SHA = "a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    research = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
    task = research / "artifacts/tasks" / HERE.name
    manifest = json.loads((HERE / "screen_a1237821.manifest.json").read_text())
    assert manifest["status"] == "complete"
    inputs = {k: Path(v["path"]) for k, v in manifest["inputs"].items()}
    assert {k: sha(p) for k, p in inputs.items()} == {k: v["sha256"] for k, v in manifest["inputs"].items()}
    assert sha(inputs["model"]) == MODEL_SHA
    screen = list(csv.DictReader((task / "screen_a1237821/gene_screen.tsv").open(), delimiter="\t"))
    wt_growth = manifest["WT"]["raw_growth"]
    medium = {r["exchange"]: float(r["uptake"]) for r in csv.DictReader(inputs["medium"].open())}
    fn_rows = [r for r in screen if r["experimental_positive"] == "True" and r["mtdna_out_of_scope"] != "True"
               and r["prediction_15pct"] == "nonessential"]
    assert len(fn_rows) == 212, len(fn_rows)

    sim = load_effective_simulation_context(model_path=inputs["model"], media_path=inputs["medium"],
                                            strain_profile_path=inputs["profile"])
    model = sim.model
    model.solver = "gurobi"
    configure_solver(model, {"solver": PARAMS})
    assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {"biomass_C": 1.0}
    lps = Counter()

    def growth():
        lps["lp"] += 1
        s = model.optimize()
        ok = s.status == "optimal" and s.objective_value is not None and math.isfinite(s.objective_value)
        return (float(s.objective_value) if ok else None), s.status

    g0, st = growth()
    assert st == "optimal" and abs(g0 - wt_growth) < 1e-9, (g0, wt_growth)

    # Whole-model blocked set (medium as set, no growth constraint).
    t = time.monotonic()
    blocked = set(find_blocked_reactions(model, open_exchanges=False, processes=1))
    lps["lp"] += 2 * len(model.reactions)
    print(f"blocked reactions: {len(blocked)} of {len(model.reactions)} ({time.monotonic()-t:.0f}s)", flush=True)

    # WT pFBA and uptake pattern.
    wt_pfba = pfba(model).fluxes
    lps["lp"] += 2
    exchanges = {r.id for r in model.reactions if r.boundary}

    # WT FVA at the optimum for every reaction the FN genes touch.
    touched = sorted({rid for r in fn_rows for rid in r["model_reaction_ids"].split(";") if rid})
    t = time.monotonic()
    fva = flux_variability_analysis(model, reaction_list=touched, fraction_of_optimum=1.0, processes=1)
    lps["lp"] += 1 + 2 * len(touched)
    print(f"WT FVA of {len(touched)} reactions ({time.monotonic()-t:.0f}s)", flush=True)

    capped = {rid for rid, v in medium.items() if 0 < v < 1000}  # glucose, iron, adenine, amino acids, uracil
    results = []
    for index, row in enumerate(fn_rows, 1):
        gid = row["gene_id"]
        gene = model.genes.get_by_id(gid)
        rxns = sorted(gene.reactions, key=lambda r: r.id)
        closed = sorted(json.loads(row["changed_reaction_bounds"]))
        kept = [r for r in rxns if r.id not in closed]
        partners = {r.id: dict(rule=r.gene_reaction_rule,
                               other_genes=sorted(g.id for g in r.genes if g.id != gid)) for r in kept}
        rxn_info = {}
        for r in rxns:
            lo, hi = float(fva.loc[r.id, "minimum"]), float(fva.loc[r.id, "maximum"])
            rxn_info[r.id] = dict(name=r.name, reaction=r.reaction, bounds=list(r.bounds),
                                  rule=r.gene_reaction_rule, wt_flux=float(wt_pfba[r.id]),
                                  fva_min=lo, fva_max=hi, blocked=r.id in blocked,
                                  closed_by_gene_ko=r.id in closed,
                                  used_at_optimum=max(abs(lo), abs(hi)) > TOL,
                                  required_at_optimum=(lo > TOL) or (hi < -TOL))
        # Control 1: reaction-level KO (every reaction of the gene closed).
        with model:
            for r in rxns:
                r.bounds = (0.0, 0.0)
            g_rxn, st_rxn = growth()
            uptake_shift, route = {}, {}
            if g_rxn is not None and g_rxn / wt_growth >= CUTOFF:
                try:
                    ko_pfba = pfba(model).fluxes
                    lps["lp"] += 2
                    for rid in exchanges:
                        d = float(ko_pfba[rid] - wt_pfba[rid])
                        # uptake is negative flux on boundary reactions; "rises" = more negative
                        if d < -TOL:
                            uptake_shift[rid] = dict(wt=float(wt_pfba[rid]), ko=float(ko_pfba[rid]),
                                                     cap=medium.get(rid), capped=rid in capped)
                    diff = (ko_pfba - wt_pfba)
                    diff = diff[diff.abs() > TOL].abs().sort_values(ascending=False)
                    route = {rid: dict(wt=float(wt_pfba[rid]), ko=float(ko_pfba[rid]))
                             for rid in diff.index[:12] if rid not in exchanges}
                except Exception as exc:  # recorded, never hidden
                    route = {"error": f"{type(exc).__name__}: {exc}"}
        ratio_gene = float(row["ko_wt_ratio"])
        ratio_rxn = g_rxn / wt_growth if g_rxn is not None else None
        closed_info = [rxn_info[c] for c in closed]
        # Category (mechanical rule; see REPORT.md for the definitions).
        if not rxns:
            cat, sub = "e_gene_in_model_without_reactions", "e1_no_reaction"
        elif not closed:
            cat = "a_or_partner_keeps_reaction_open"
            sub = "a1_reaction_ko_lethal" if (ratio_rxn is not None and ratio_rxn < CUTOFF) or st_rxn != "optimal" \
                else "a2_reaction_ko_still_grows"
        elif all(c["blocked"] for c in closed_info):
            cat, sub = "d_closed_reactions_blocked_in_wt", "d1_structurally_blocked"
        elif not any(c["used_at_optimum"] for c in closed_info):
            cat, sub = "b_closed_reactions_unused_at_optimum", "b1_route_optional"
        else:
            rises = [rid for rid, v in uptake_shift.items() if v["capped"] and rid != "R1070"]
            if ratio_gene < 0.999 and rises:
                cat, sub = "c_medium_supplies_product", "c1_growth_capped_by_uptake"
            elif rises:
                cat, sub = "c_medium_supplies_product", "c2_uptake_rises_no_growth_cost"
            else:
                cat, sub = "b_alternative_route_carries_flux", ("b2_partial_growth_cost" if ratio_gene < 0.999 else "b3_no_growth_cost")
        partial = ratio_gene < 0.999
        results.append(dict(gene_id=gid, ko_wt_ratio=ratio_gene, partial_growth_effect=partial,
            category=cat, subcategory=sub, reactions=[r.id for r in rxns], closed_by_gene_ko=closed,
            or_partners=partners, reaction_ko_ratio=ratio_rxn, reaction_ko_status=st_rxn,
            uptake_shift=uptake_shift, route_change=route, reaction_info=rxn_info,
            model_name_unverified=row["model_name_unverified"],
            list_function=row["source_putative_function"].strip()))
        if index % 40 == 0 or index == len(fn_rows):
            print(f"{index}/{len(fn_rows)} genes; LPs {lps['lp']}", flush=True)

    # Restore check.
    g1, _ = growth()
    assert abs(g1 - wt_growth) < 1e-9
    assert {k: sha(p) for k, p in inputs.items()} == {k: v["sha256"] for k, v in manifest["inputs"].items()}

    out = task / "fn_classification"
    out.mkdir(exist_ok=True)
    (out / "fn_classification.json").write_text(json.dumps(dict(
        model_sha256=MODEL_SHA, wt_growth=wt_growth, cutoff=CUTOFF, tolerance=TOL, solver=PARAMS,
        lps=lps["lp"], blocked_reaction_count=len(blocked), blocked_reactions=sorted(blocked),
        category_counts=dict(Counter(r["category"] for r in results)),
        subcategory_counts=dict(Counter(r["subcategory"] for r in results)),
        genes=results), indent=1, allow_nan=False) + "\n")
    fields = ["gene_id", "ko_wt_ratio", "category", "subcategory", "reactions", "closed_by_gene_ko",
              "or_partner_genes", "reaction_ko_ratio", "reaction_ko_status", "closed_wt_fva",
              "capped_uptakes_rising", "route_top_changes", "list_function"]
    with (HERE / "fn_causes_a1237821.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow(dict(gene_id=r["gene_id"], ko_wt_ratio=f"{r['ko_wt_ratio']:.4f}", category=r["category"],
                subcategory=r["subcategory"], reactions=";".join(r["reactions"]),
                closed_by_gene_ko=";".join(r["closed_by_gene_ko"]),
                or_partner_genes=";".join(sorted({g for p in r["or_partners"].values() for g in p["other_genes"]})),
                reaction_ko_ratio="" if r["reaction_ko_ratio"] is None else f"{r['reaction_ko_ratio']:.4f}",
                reaction_ko_status=r["reaction_ko_status"],
                closed_wt_fva=";".join(f"{c}:[{r['reaction_info'][c]['fva_min']:.3g},{r['reaction_info'][c]['fva_max']:.3g}]"
                                       + ("(blocked)" if r['reaction_info'][c]['blocked'] else "") for c in r["closed_by_gene_ko"]),
                capped_uptakes_rising=";".join(f"{k}:{v['wt']:.3g}->{v['ko']:.3g}" for k, v in r["uptake_shift"].items() if v["capped"]),
                route_top_changes=";".join(f"{k}:{v['wt']:.3g}->{v['ko']:.3g}" for k, v in list(r["route_change"].items())[:6] if isinstance(v, dict)),
                list_function=r["list_function"][:160]))
    print(json.dumps(dict(Counter(r["category"] for r in results)), indent=1))
    print(json.dumps(dict(Counter(r["subcategory"] for r in results)), indent=1))
    print("LPs", lps["lp"])


if __name__ == "__main__":
    main()
