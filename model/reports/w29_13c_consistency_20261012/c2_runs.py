"""Task C2: in-memory effects of making R825 (Crc1 carnitine/acetylcarnitine antiport) reversible.

Diagnostic only; no model is written. Requested by the user (2026-10-13):
1. WT growth, R825 as supplied vs reversible: W29 V2, and PO1f + CSM-Leu at uptake 2.414 (with and
   without the holistic overlay).
2. Full gene benchmark before / after (PO1f + CSM-Leu screen context, model as supplied, benchmark
   rule; second labelled view counts "no solution" as lethal), listing every changed call.
3. Complex I (R2062) KO/WT with and without PDH >= 0.866629 x glucose uptake, grid over medium
   {V2, V1} x overlay {no, yes} x 13C caps {none, module} x R825 {as is, reversible}; FVA of the
   acetyl-CoA export routes for the PDH-floor runs.
4. Step-5 reference check (Amendment 1): candidate + overlay + module caps, glucose, V2 / V1, growth
   ratio vs 0.99 / 1.00 and V2 NDH2 ~204, complex I ~187 per 100 glucose (2 % stop on growth ratio).
5. pFBA vectors for the maps (acetyl-CoA routes; elastic-fit representative points for C1).

Run from the repository root:
    IYALI26_RESEARCH_ROOT=... .venv/bin/python -B model/reports/w29_13c_consistency_20261012/c2_runs.py
Writes c2_runs.json and CSV tables here; gene screens and flux vectors go to
$IYALI26_RESEARCH_ROOT/artifacts/tasks/w29_13c_consistency_20261012/c2/. Never overwrites.
"""
import csv
import json
import math
import platform
import subprocess
import time
from collections import Counter
from datetime import datetime, timezone

import cobra
import gurobipy
from cobra.flux_analysis import pfba

import common
from common import (ART, HERE, MU, PARAMS, PINS, REPO, UPTAKE, Counter_, add_caps, apply_overlay, check_pins, fva,
                    growth, load, load_caps, pdh_floor, pfba_vector, r825_reversible)
from elastic_fit import MAPPING, ranges
from tools.dfba_new_fn_essentiality import load_experimental_reference

OUT = ART / "c2"
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
ROUTES = ("R1311", "R693", "R236", "R817", "R831", "R832", "R833", "R807", "R902", "R1386", "R1894", "R825",
          "R1275", "R1205", "R193", "R191", "R83", "R412", "R1172", "R410", "R805", "R945", "R1974", "R89", "R1708",
          "R690", "R538", "R691", "R1378", "R134", "R90", "R2062", "R570", "R305")


def per100(vec, uptake):
    return {k: 100 * v / uptake for k, v in vec.items()}


def metrics(rows, cutoff, excluded, lethal_view=False):
    """Task B scoring (run_task_b.py), unchanged."""
    c = Counter(TP=0, FN=0, FP=0, TN=0, unresolved=0)
    for row in rows:
        if row["gene_id"] in excluded:
            continue
        ratio = row["ratio"]
        if ratio is None and lethal_view:
            ratio = 0.0
        if ratio is None:
            c["unresolved"] += 1
            continue
        essential = ratio < cutoff
        c[("TP" if row["positive"] else "FP") if essential else ("FN" if row["positive"] else "TN")] += 1
    tp, fn, fp, tn = (c[k] for k in ("TP", "FN", "FP", "TN"))
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return dict(cutoff=cutoff, view="no solution counted as lethal" if lethal_view else "benchmark rule", **c,
                recall=tp / (tp + fn), precision=tp / (tp + fp) if tp + fp else None, specificity=tn / (tn + fp),
                mcc=(tp * tn - fp * fn) / den if den else None)


def screen(sim, model, positives, mtdna):
    wt, status = growth(model)
    assert wt and wt > 0, status
    genes = sorted(g.id for g in model.genes if g.id not in sim.excluded_runtime_genes)
    assert len(genes) == 1073, len(genes)
    rows = []
    for gene_id in genes:
        gene = model.genes.get_by_id(gene_id)
        with model:
            gene.knock_out()
            g, s = growth(model)
        rows.append({"gene_id": gene_id, "positive": gene_id in positives, "mtdna": gene_id in mtdna,
                     "reactions": ";".join(sorted(r.id for r in gene.reactions)), "status": s, "growth": g,
                     "ratio": None if g is None else g / wt})
    return wt, rows


def fit_point(model, keys):
    """Elastic-fit optimum, then pFBA at fixed biomass on the optimal face (one representative vector)."""
    rng = ranges()
    prob = model.problem
    model.reactions.biomass_C.bounds = (MU, MU)
    extra, terms = [], []
    for key, _, _, coeffs in MAPPING:
        if key not in keys:
            continue
        lb, ub, _ = rng[key]
        e = (100.0 / UPTAKE) * sum(c * model.reactions.get_by_id(r).flux_expression for r, c in coeffs.items())
        sp, sm = prob.Variable(f"sp_{key}", lb=0), prob.Variable(f"sm_{key}", lb=0)
        extra += [sp, sm, prob.Constraint(e - sp, ub=ub), prob.Constraint(e + sm, lb=lb)]
        terms.append((sp + sm) / (ub - lb))
    model.add_cons_vars(extra)
    total = sum(terms)
    model.objective = prob.Objective(total, direction="min")
    Counter_.lps += 1
    model.slim_optimize()
    opt = model.solver.objective.value
    model.add_cons_vars([prob.Constraint(total, ub=opt + 1e-6, name="face")])
    model.objective = "biomass_C"
    Counter_.lps += 2
    return opt, pfba(model).fluxes.to_dict()


def main():
    inputs = check_pins()
    caps = load_caps()
    assert not OUT.exists(), OUT
    started = time.monotonic()
    positives, _ = load_experimental_reference(PINS["experimental"][0])
    positives = set(positives)
    mtdna = json.loads(PINS["mtdna_scope"][0].read_text())["genes"]
    record = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "command": "python -B model/reports/w29_13c_consistency_20261012/c2_runs.py",
              "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
              "inputs": inputs, "solver": PARAMS,
              "versions": {"python": platform.python_version(), "cobra": cobra.__version__,
                           "gurobi": list(gurobipy.gurobi.version())},
              "pdh_floor_fraction": common.PDH_FLOOR,
              "scoring": "mtDNA genes out of scope (benchmark.md 2026-10-08); negatives are model genes absent from "
                         "the positive list (user-defined class); benchmark rule: no solution = unresolved; second "
                         "view counts no solution as lethal (labelled); development reference, not independent "
                         "validation"}
    models = {m: load(m) for m in ("V2", "V1", "CSM")}
    for m, (sim, _) in models.items():
        record.setdefault("simulation_context", {})[m] = sim.provenance()

    # 1. WT growth.
    wt = []
    for medium in ("V2", "CSM"):
        sim, model = models[medium]
        for overlay in (False, True):
            for rev in (False, True):
                with model:
                    if overlay:
                        apply_overlay(model)
                    if rev:
                        r825_reversible(model)
                    g, s = growth(model)
                    uptake = -model.reactions.R1070.flux if g else None
                wt.append({"medium": medium, "overlay": overlay, "R825": "reversible" if rev else "as supplied",
                           "growth": g, "status": s, "glucose_uptake": uptake})
    record["wt_growth"] = wt

    # 2. Benchmark (PO1f + CSM-Leu screen context, model as supplied).
    sim, model = models["CSM"]
    screens = {}
    for rev in (False, True):
        label = "R825_reversible" if rev else "as_supplied"
        with model:
            if rev:
                r825_reversible(model)
            w, rows = screen(sim, model, positives, mtdna)
            controls = {}
            for name, rid in (("complex_I", "R2062"), ("complex_III", "R305"), ("complex_IV", "R304"),
                              ("ATP_synthase", "R171")):
                with model:
                    model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                    g, s = growth(model)
                controls[name] = {"reaction": rid, "growth": g, "status": s, "ratio": None if g is None else g / w}
        screens[label] = rows
        excluded = frozenset(mtdna)
        record.setdefault("benchmark", {})[label] = {
            "WT": w, "controls": controls,
            "scores": [metrics(rows, c, excluded, v) for v in (False, True) for c in CUTOFFS]}
        print("benchmark", label, w, {k: round(v["ratio"], 4) if v["ratio"] is not None else None
                                      for k, v in controls.items()}, flush=True)
    base = {r["gene_id"]: r for r in screens["as_supplied"]}
    diffs, calls = [], {}
    for r in screens["R825_reversible"]:
        b = base[r["gene_id"]]
        if (b["ratio"] is None) != (r["ratio"] is None) or (
                r["ratio"] is not None and abs(r["ratio"] - b["ratio"]) > 1e-6):
            diffs.append({"gene_id": r["gene_id"], "positive": r["positive"], "mtdna": r["mtdna"],
                          "reactions": r["reactions"], "ratio_as_supplied": b["ratio"],
                          "ratio_R825_reversible": r["ratio"], "status_as_supplied": b["status"],
                          "status_R825_reversible": r["status"]})
    for c in CUTOFFS:
        calls[str(c)] = [d["gene_id"] for d in diffs if d["gene_id"] not in mtdna and (
            (d["ratio_as_supplied"] is None) != (d["ratio_R825_reversible"] is None)
            or (d["ratio_as_supplied"] is not None and (d["ratio_as_supplied"] < c) != (d["ratio_R825_reversible"] < c)))]
    record["benchmark_ratio_changes"] = diffs
    record["benchmark_call_changes"] = calls

    # 3. Complex I x PDH floor grid, and 4. step-5 reference check.
    grid, vectors = [], {}
    for medium in ("V2", "V1"):
        sim, model = models[medium]
        for overlay in (False, True):
            for capped in (False, True):
                for rev in (False, True):
                    for floor in (False, True):
                        entry = {"medium": medium, "overlay": overlay, "caps": "module" if capped else None,
                                 "R825": "reversible" if rev else "as supplied", "pdh_floor": floor}
                        for label, ko in (("WT", None), ("complex_I_KO", "R2062")):
                            with model:
                                if overlay:
                                    apply_overlay(model)
                                if capped:
                                    add_caps(model, caps, "module")
                                if rev:
                                    r825_reversible(model)
                                if floor:
                                    pdh_floor(model)
                                if ko:
                                    model.reactions.get_by_id(ko).bounds = (0.0, 0.0)
                                g, s = growth(model)
                                entry[label] = {"growth": g, "status": s}
                                if g and (floor or (overlay and capped and not rev)):
                                    f = fva(model, ROUTES)
                                    entry[label]["fva_per_100"] = {k: [100 * v[0] / UPTAKE, 100 * v[1] / UPTAKE]
                                                                   for k, v in f.items()} if "error" not in f else f
                                if g and medium == "V2" and overlay and capped:
                                    v = pfba_vector(model)
                                    vectors[f"{label}|{'rev' if rev else 'asis'}|{'floor' if floor else 'nofloor'}"] = v
                        wt_g, ko_g = entry["WT"]["growth"], entry["complex_I_KO"]["growth"]
                        entry["complex_I_KO_over_WT"] = None if not wt_g or ko_g is None else ko_g / wt_g
                        grid.append(entry)
                        print("grid", medium, overlay, capped, rev, floor, wt_g, ko_g, entry["complex_I_KO_over_WT"],
                              flush=True)
    record["complex_I_grid"] = grid
    ref = {"V2": 0.99, "V1": 1.00}
    checks = []
    for medium, target in ref.items():
        e = next(x for x in grid if x["medium"] == medium and x["overlay"] and x["caps"] == "module"
                 and x["R825"] == "as supplied" and not x["pdh_floor"])
        ratio = e["WT"]["growth"] / MU
        checks.append({"medium": medium, "target": target, "observed_ratio": ratio,
                       "relative_difference": ratio / target - 1, "within_2pct": abs(ratio / target - 1) <= 0.02,
                       "NDH2_R570_per100": e["WT"]["fva_per_100"]["R570"],
                       "complex_I_R2062_per100": e["WT"]["fva_per_100"]["R2062"]})
    record["step5_reference_check"] = checks
    record["step5_stop_fired"] = not all(c["within_2pct"] for c in checks)

    # 5. Elastic-fit representative points (V2, overlay, all 22 ranges) for the C1 map.
    sim, model = models["V2"]
    for rev in (False, True):
        with model:
            apply_overlay(model)
            if rev:
                r825_reversible(model)
            with model:
                opt, vec = fit_point(model, [k for k, *_ in MAPPING])
        vectors[f"fit|{'rev' if rev else 'asis'}"] = vec
        record.setdefault("fit_points", {})["rev" if rev else "asis"] = opt
    OUT.mkdir(parents=True, exist_ok=False)
    for label, rows in screens.items():
        with (OUT / f"gene_screen_{label}.tsv").open("x") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    (OUT / "pfba_vectors.json").open("x").write(json.dumps(vectors) + "\n")

    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started,
                  inputs_unchanged=all(common.sha(p) == d for p, d in PINS.values()),
                  outputs={p.name: common.sha(p) for p in sorted(OUT.iterdir())})
    (HERE / "c2_runs.json").open("x").write(json.dumps(record, indent=1, default=str) + "\n")
    with (HERE / "c2_complex_I_grid.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["medium", "overlay", "caps", "R825", "pdh_floor", "WT_growth", "WT_status", "KO_growth",
                    "KO_status", "complex_I_KO_over_WT"])
        for e in grid:
            w.writerow([e["medium"], e["overlay"], e["caps"], e["R825"], e["pdh_floor"], e["WT"]["growth"],
                        e["WT"]["status"], e["complex_I_KO"]["growth"], e["complex_I_KO"]["status"],
                        e["complex_I_KO_over_WT"]])
    with (HERE / "c2_benchmark_scores.csv").open("x", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["model", "view", "cutoff", "TP", "FP", "FN", "TN", "unresolved", "recall", "precision",
                    "specificity", "mcc"])
        for label, b in record["benchmark"].items():
            for s in b["scores"]:
                w.writerow([label, s["view"], s["cutoff"], s["TP"], s["FP"], s["FN"], s["TN"], s["unresolved"],
                            s["recall"], s["precision"], s["specificity"], s["mcc"]])
    with (HERE / "c2_benchmark_ratio_changes.csv").open("x", newline="") as f:
        fields = ["gene_id", "positive", "mtdna", "reactions", "ratio_as_supplied", "ratio_R825_reversible",
                  "status_as_supplied", "status_R825_reversible"]
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(diffs)
    print("call changes", calls, "| ratio changes", len(diffs), "| step5 checks", checks, "| LPs", Counter_.lps,
          "| elapsed", round(record["elapsed_seconds"]), "s")


if __name__ == "__main__":
    main()
