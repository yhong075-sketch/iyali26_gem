"""Targeted controls for the false-negative causes (investigation only, in-memory, no model change).

Same model, medium, profile and solver as run_screen.py / classify_fn.py (SHA-checked through the manifest).

Part 1, every FN gene (212):
  - control "aa_closed": the 11 amino-acid and the adenine uptakes of the Task A medium set to 0 (uracil stays,
    PO1f is ura3-), WT and gene KO re-solved; ratio KO/WT on that medium;
  - control "vit_closed": the 9 vitamin uptakes set to 0, same;
  - dead-end check of every reaction the gene KO closes that is blocked in WT: participants with no other
    producer or consumer anywhere in the model (structural reason for the block).
Part 2, deep batch (DEEP): reaction-level KO, then FVA (fraction_of_optimum 0.99 of the KO optimum) over all
  reactions; reactions required in the KO (|flux| >= 1e-6 throughout) that are not required in WT are the
  bypass route. WT FVA at 0.99 is computed once for comparison.
Part 3: GPR text of the complex / multi-gene reactions named in REPORT.md, and WT FVA of the R659 loop.

Usage (repository root, IYALI26_RESEARCH_ROOT set):
    python -B model/reports/essentiality_fn_causes_20261027/trace_controls.py
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

from cobra.flux_analysis import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.diagnose_closed_energy import configure_solver

HERE = Path(__file__).resolve().parent
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0, FeasibilityTol=1e-9, OptimalityTol=1e-9, IntFeasTol=1e-9)
TOL = 1e-6
MODEL_SHA = "a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82"
AA_ADENINE = ["R1003", "R1202", "R1204", "R1215", "R1217", "R1220", "R1222", "R1223", "R1231", "R1232", "R1233", "R1234"]
VITAMINS = ["R1029", "R1340", "R1305", "R2270", "R946", "R1263", "R971", "R1127", "R1268", "R1314"]
DEEP = ["YALI1C08262g", "YALI1F11049g", "YALI1F38116g", "YALI1D18130g", "YALI1A03069g", "YALI1B05798g",
        "YALI1D12493g", "YALI1A21156g", "YALI1F00821g", "YALI1D11984g", "YALI1B28659g", "YALI1E30886g",
        "YALI1F27556g", "YALI1C09033g", "YALI1F18105g", "YALI1B01136g"]
GPR_REACTIONS = ["R171", "R794", "R211", "R2193", "R759", "R760", "R723", "R673", "R251", "R737", "R177", "R815",
                 "R121", "R1085", "R60", "R540", "R712", "R4", "R289", "R1815", "R2081", "R334", "R239", "R1556",
                 "R419", "R420", "R88", "R2062", "R740", "R693", "R344", "R411", "R412", "R659", "r0210", "R_CYST_II",
                 "R793", "R532", "R548", "R166", "R568", "R569", "R697", "R698", "R306", "R14", "R274", "R559"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    research = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
    task = research / "artifacts/tasks" / HERE.name
    manifest = json.loads((HERE / "screen_a1237821.manifest.json").read_text())
    inputs = {k: Path(v["path"]) for k, v in manifest["inputs"].items()}
    assert {k: sha(p) for k, p in inputs.items()} == {k: v["sha256"] for k, v in manifest["inputs"].items()}
    assert sha(inputs["model"]) == MODEL_SHA
    wt_growth = manifest["WT"]["raw_growth"]
    cls = json.loads((task / "fn_classification/fn_classification.json").read_text())
    genes = cls["genes"]
    blocked = set(cls["blocked_reactions"])
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

    g0, _ = growth()
    assert abs(g0 - wt_growth) < 1e-9
    for rid in AA_ADENINE + VITAMINS:
        assert model.reactions.get_by_id(rid).lower_bound < 0, rid  # uptake is open in the medium

    def close(ids):
        for rid in ids:
            model.reactions.get_by_id(rid).lower_bound = 0.0

    # Part 1 -------------------------------------------------------------------------------------------------
    with model:
        close(AA_ADENINE)
        wt_aa, st = growth()
        assert st == "optimal"
    with model:
        close(VITAMINS)
        wt_vit, st = growth()
        assert st == "optimal"
    print(f"WT on Task A medium {wt_growth:.7f}; amino acids+adenine closed {wt_aa:.7f}; vitamins closed {wt_vit:.7f}", flush=True)
    part1 = {}
    for g in genes:
        gid = g["gene_id"]
        gene = model.genes.get_by_id(gid)
        out = {}
        for label, ids, wt_ref in (("aa_closed", AA_ADENINE, wt_aa), ("vit_closed", VITAMINS, wt_vit)):
            with model:
                close(ids)
                gene.knock_out()
                gk, st = growth()
            if gene._model is None:
                gene.functional = True
            out[label] = dict(status=st, growth=gk, ratio=(gk / wt_ref if gk is not None else None),
                              essential_15=(gk is not None and gk / wt_ref < 0.15) or st != "optimal")
        dead = {}
        for rid in g["closed_by_gene_ko"]:
            if rid not in blocked:
                continue
            r = model.reactions.get_by_id(rid)
            ends = []
            for met in r.metabolites:
                others = [x for x in met.reactions if x.id != rid]
                if not others:
                    ends.append(f"{met.id} (no other reaction)")
                elif all(x.id in blocked for x in others):
                    ends.append(f"{met.id} (all {len(others)} other reactions blocked)")
            dead[rid] = ends
        out["blocked_reaction_dead_ends"] = dead
        part1[gid] = out
    print(f"Part 1 done; LPs {lps['lp']}", flush=True)

    # Part 2 -------------------------------------------------------------------------------------------------
    t = time.monotonic()
    wt_fva = flux_variability_analysis(model, fraction_of_optimum=0.99, processes=1)
    lps["lp"] += 1 + 2 * len(model.reactions)
    wt_required = {rid for rid in wt_fva.index if min(abs(wt_fva.loc[rid, "minimum"]), abs(wt_fva.loc[rid, "maximum"])) >= TOL
                   and wt_fva.loc[rid, "minimum"] * wt_fva.loc[rid, "maximum"] > 0}
    print(f"WT FVA(0.99): {len(wt_required)} required reactions ({time.monotonic()-t:.0f}s)", flush=True)
    part2 = {}
    by_id = {g["gene_id"]: g for g in genes}
    for gid in DEEP:
        g = by_id[gid]
        rxns = g["reactions"]
        with model:
            for rid in rxns:
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            gk, st = growth()
            entry = dict(reactions=rxns, reaction_ko_growth=gk, status=st, ratio=(gk / wt_growth if gk else None))
            if gk and gk > TOL:
                t = time.monotonic()
                ko_fva = flux_variability_analysis(model, fraction_of_optimum=0.99, processes=1)
                lps["lp"] += 1 + 2 * len(model.reactions)
                ko_required = {rid for rid in ko_fva.index
                               if min(abs(ko_fva.loc[rid, "minimum"]), abs(ko_fva.loc[rid, "maximum"])) >= TOL
                               and ko_fva.loc[rid, "minimum"] * ko_fva.loc[rid, "maximum"] > 0}
                new = sorted(ko_required - wt_required)
                lost = sorted(wt_required - ko_required)
                entry["newly_required"] = {rid: dict(name=model.reactions.get_by_id(rid).name,
                    reaction=model.reactions.get_by_id(rid).reaction, rule=model.reactions.get_by_id(rid).gene_reaction_rule,
                    ko_min=float(ko_fva.loc[rid, "minimum"]), ko_max=float(ko_fva.loc[rid, "maximum"]),
                    wt_min=float(wt_fva.loc[rid, "minimum"]), wt_max=float(wt_fva.loc[rid, "maximum"])) for rid in new}
                entry["no_longer_required"] = lost
                print(f"  {gid}: KO/WT {entry['ratio']:.4f}; newly required {len(new)}, no longer required {len(lost)} ({time.monotonic()-t:.0f}s)", flush=True)
            else:
                print(f"  {gid}: reaction KO growth {gk} ({st})", flush=True)
        part2[gid] = entry

    # Part 3 -------------------------------------------------------------------------------------------------
    gprs = {}
    for rid in GPR_REACTIONS:
        if rid in model.reactions:
            r = model.reactions.get_by_id(rid)
            gprs[rid] = dict(name=r.name, reaction=r.reaction, bounds=list(r.bounds), rule=r.gene_reaction_rule,
                             n_genes=len(r.genes), blocked=rid in blocked,
                             wt_fva_099=[float(wt_fva.loc[rid, "minimum"]), float(wt_fva.loc[rid, "maximum"])])
    loop = {}
    for rid in ("R659", "r0210", "R_CYST_II", "R793", "R239"):
        if rid in model.reactions:
            with model:
                model.objective = rid
                model.objective_direction = "max"
                mx, _ = growth()
                model.objective_direction = "min"
                mn, _ = growth()
            loop[rid] = dict(unconstrained_min=mn, unconstrained_max=mx)
    g1, _ = growth()
    assert abs(g1 - wt_growth) < 1e-9
    (task / "fn_classification/trace_controls.json").write_text(json.dumps(dict(
        model_sha256=MODEL_SHA, wt_growth=wt_growth, wt_aa_closed=wt_aa, wt_vit_closed=wt_vit,
        aa_adenine_exchanges=AA_ADENINE, vitamin_exchanges=VITAMINS, lps=lps["lp"],
        part1=part1, part2=part2, gprs=gprs, loop_ranges=loop), indent=1, allow_nan=False) + "\n")
    print(json.dumps(dict(aa_closed_essential=sum(v["aa_closed"]["essential_15"] for v in part1.values()),
                          vit_closed_essential=sum(v["vit_closed"]["essential_15"] for v in part1.values()),
                          lps=lps["lp"])))


if __name__ == "__main__":
    main()
