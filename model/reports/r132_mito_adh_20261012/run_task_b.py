"""Task B runs: R132-removal candidate vs baseline.

1. Controls: WT and KO/WT of R2062, R305, R304, R171 (Task A medium + PO1f profile).
2. Full gene benchmark for both models, scored as in Task A (benchmark rule and "lethal" view).
3. Route analysis: FVA at the optimum plus one pFBA vector, WT and complex I KO, both models.
4. Step 5 (diagnostic): 13C caps from model/conditions/assumptions/w29_13c_upper_caps_20261012.json
   on glucose, glycerol and oleate, plus a labelled reconciliation on an amino-acid-free medium.

Run from the repository root:
    .venv/bin/python -B model/reports/r132_mito_adh_20261012/run_task_b.py
Writes run.json and small CSV tables here; per-gene screens and flux vectors go to
$IYALI26_RESEARCH_ROOT/artifacts/tasks/r132_mito_adh_20261012/runs/. Never overwrites.
"""
import csv
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO / "platform"))

import cobra
import gurobipy
from cobra.flux_analysis import flux_variability_analysis, pfba
from cobra.util.solver import linear_reaction_coefficients
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from tools.dfba_new_fn_essentiality import load_experimental_reference
from tools.diagnose_closed_energy import configure_solver

RESEARCH = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
OUT = RESEARCH / "artifacts/tasks" / HERE.name / "runs"
MODELS = {
    "baseline": (REPO / "model/candidates/E5_coq9_respiratory_package_20261008.xml",
                 "0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba"),
    "candidate": (REPO / "model/candidates/E5_coq9_resp_pkg_r132_20261012.xml",
                  "a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515"),
}
MEDIA = {
    "glucose": (REPO / "model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv",
                "6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45"),
    "glycerol": (REPO / "model/conditions/media/po1f_csm_leu_w29glycerol_20261012.csv",
                 "23371f0aad1e3b464bc6eef9099dbde53aa134d86f057f7a7d98fda63dc98c4d"),
    "oleic_acid": (REPO / "model/conditions/media/po1f_csm_leu_w29oleate_20261012.csv",
                   "707f15d80b7d1bbe8634bc3fcbfbde4337b8a509ba11d20a75022b4f4c34b7e4"),
}
INPUTS = {
    "profile": (RESEARCH / "state/strain_profiles/po1f_sd_leu.json",
                "35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383"),
    "experimental": (RESEARCH / "state/essentiality/repository/consensus_essential_genes.csv",
                     "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b"),
    "mtdna_scope": (REPO / "model/expected/mtdna_genes_out_of_scope.json",
                    "16b3e9b4418f24821b8775936fd377ae21f77afa8e8af2b8949a246f62c44996"),
    "caps": (REPO / "model/conditions/assumptions/w29_13c_upper_caps_20261012.json", None),
    "w29_fluxes": (RESEARCH / "reference/w29_worland2024/glucose_central_fluxes.csv",
                   "844615d0edfd8ffb033c1bce15ed746e9bb8ea02cbbe75016771b5c5fd08389f"),
}
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0, FeasibilityTol=1e-7, OptimalityTol=1e-7, IntFeasTol=1e-7)
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
COMPLEXES = {"complex_I": "R2062", "complex_III": "R305", "complex_IV": "R304", "ATP_synthase": "R171"}
AMINO_ACID_ROWS = ("R1003", "R1202", "R1204", "R1215", "R1217", "R1220", "R1222", "R1223", "R1231", "R1232",
                   "R1233", "R1234")
# Shuttle and replacement-route reactions (step 4 maps and FVA).
ROUTE = ("R130", "R1377", "R_ALCDH_nadp_hi", "R1114", "R132", "R997", "R135", "R2070", "R136", "R1113", "R691",
         "R690", "R534", "R_MDHy", "R807", "R831", "R902", "R1386", "R533", "R538", "R693", "R1311", "R570",
         "R2062", "R348", "R347", "R614", "R361", "R357", "R360", "R362", "R305", "R304", "R171")
STEP5_FVA = ("R570", "R2062", "R690", "R807", "R831", "R902", "R1386", "R693", "R614", "R361", "R357", "R360",
             "R362", "R132")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class Counter_:
    lps = 0


def context(model_key, medium_key):
    sim = load_effective_simulation_context(model_path=MODELS[model_key][0], media_path=MEDIA[medium_key][0],
                                            strain_profile_path=INPUTS["profile"][0])
    model = sim.model
    model.solver = "gurobi"
    configure_solver(model, {"solver": PARAMS})
    assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {"biomass_C": 1.0}
    return sim, model


def growth(model):
    Counter_.lps += 1
    solution = model.optimize()
    value = solution.objective_value
    ok = solution.status == "optimal" and value is not None and math.isfinite(value) and value >= 0
    return (float(value) if ok else None), solution.status


def fva(model, reactions):
    present = [r for r in reactions if r in model.reactions]
    Counter_.lps += 1 + 2 * len(present)
    try:
        frame = flux_variability_analysis(model, reaction_list=present, fraction_of_optimum=1.0, processes=1)
        return {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in present}
    except Exception as exc:  # recorded, not hidden
        return {"error": f"{type(exc).__name__}: {exc}"}


def pfba_vector(model):
    Counter_.lps += 2
    try:
        return pfba(model).fluxes.to_dict()
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def metrics(rows, cutoff, excluded, lethal_view=False):
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


def screen(model_key, positives, mtdna):
    sim, model = context(model_key, "glucose")
    wt, status = growth(model)
    assert wt and wt > 0, f"WT infeasible or zero: {model_key} {status}"
    result = {"simulation_context": sim.provenance(), "WT": wt, "reaction_ko": {}}
    for label, rid in COMPLEXES.items():
        with model:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            g, s = growth(model)
        result["reaction_ko"][label] = {"reaction": rid, "growth": g, "status": s, "ratio": None if g is None else g / wt}
    genes = sorted(g.id for g in model.genes if g.id not in sim.excluded_runtime_genes)
    assert len(genes) == 1073, len(genes)
    rows = []
    for gene_id in genes:
        gene = model.genes.get_by_id(gene_id)
        original = gene.functional
        with model:
            gene.knock_out()
            g, s = growth(model)
        if gene._model is None:
            gene.functional = original
        rows.append({"gene_id": gene_id, "positive": gene_id in positives, "mtdna": gene_id in mtdna,
                     "reactions": ";".join(sorted(r.id for r in gene.reactions)), "status": s,
                     "growth": g, "ratio": None if g is None else g / wt})
    excluded = frozenset(mtdna)
    result["scores"] = [metrics(rows, c, excluded, view) for view in (False, True) for c in CUTOFFS]
    # Route analysis: WT and complex I KO.
    result["route"] = {}
    for label, ko in (("WT", None), ("complex_I_KO", "R2062")):
        with model:
            if ko:
                model.reactions.get_by_id(ko).bounds = (0.0, 0.0)
            g, s = growth(model)
            result["route"][label] = {"growth": g, "status": s, "fva": fva(model, ROUTE) if g else None,
                                      "pfba": pfba_vector(model) if g else None}
    return result, rows


def verify_cap_fractions(caps):
    rows = {r["ID"]: r for r in csv.DictReader(INPUTS["w29_fluxes"][0].open())}
    for sub in caps["substrates"]:
        checks = {("pyruvate_carboxylase", "fraction"): float(rows["R31"][f"{sub}_UB"]) / 100,
                  ("net_malate_import", "fraction"): float(rows["R47"][f"{sub}_UB"]) / 100,
                  ("glycine_cleavage_R361_as_specified", "fraction"): -float(rows["R51 net"][f"{sub}_LB"]) / 100,
                  ("glycine_cleavage_R361_as_specified", "looser_fraction_used_earlier"): float(rows["R51 net"][f"{sub}_UB"]) / 100,
                  ("glycine_cleavage_module", "fraction"): -float(rows["R51 net"][f"{sub}_LB"]) / 100}
        for (cap, field), value in checks.items():
            assert abs(caps["caps"][cap][field][sub] - value) < 1e-9, (sub, cap, field, value)


def add_caps(model, caps, substrate, gcv):
    uptake = -model.reactions.get_by_id(caps["substrates"][substrate]["exchange"]).flux_expression
    chosen = [("cap_pc", "pyruvate_carboxylase", "fraction"), ("cap_malate", "net_malate_import", "fraction")]
    chosen += {"R361": [("cap_gcv_R361", "glycine_cleavage_R361_as_specified", "fraction")],
               "R361_looser": [("cap_gcv_R361_looser", "glycine_cleavage_R361_as_specified", "looser_fraction_used_earlier")],
               "module": [("cap_gcv_module", "glycine_cleavage_module", "fraction")]}[gcv]
    constraints = []
    for name, cap, field in chosen:
        spec = caps["caps"][cap]
        expression = sum(c * model.reactions.get_by_id(rid).flux_expression for rid, c in spec["reactions"].items())
        constraints.append(model.problem.Constraint(expression - spec[field][substrate] * uptake, ub=0, name=name))
    model.add_cons_vars(constraints)
    return {name: caps["caps"][cap][field][substrate] for name, cap, field in chosen}


def step5(caps):
    variants = (("baseline_as_supplied", "baseline", None), ("candidate_no_caps", "candidate", None),
                ("candidate_caps_R361", "candidate", "R361"), ("candidate_caps_R361_looser", "candidate", "R361_looser"),
                ("candidate_caps_gcv_module", "candidate", "module"))
    table, fluxes = [], {}
    for substrate, info in caps["substrates"].items():
        for name, model_key, gcv in variants:
            sim, model = context(model_key, substrate)
            exchange = info["exchange"]
            assert math.isclose(-model.reactions.get_by_id(exchange).lower_bound, info["measured_uptake"])
            with model:
                applied = add_caps(model, caps, substrate, gcv) if gcv else {}
                g, s = growth(model)
                row = {"substrate": substrate, "variant": name, "model": model_key, "caps": applied, "growth": g,
                       "status": s, "measured_growth": info["measured_growth"],
                       "predicted_over_measured": None if g is None else g / info["measured_growth"]}
                if g:
                    vector = pfba_vector(model)
                    uptake = -vector[exchange] if "error" not in vector else None
                    row.update(uptake_pfba=uptake, fva=fva(model, STEP5_FVA),
                               pfba_per_100_substrate={rid: 100 * vector[rid] / uptake for rid in STEP5_FVA
                                                       if rid in vector} if uptake else None)
                    fluxes[f"{substrate}|{name}"] = vector
            table.append(row)
    # Reconciliation with the recalled values, which came from an amino-acid-free medium (labelled; glucose only).
    for name, model_key, gcv in (("baseline_as_supplied", "baseline", None), ("candidate_caps_R361", "candidate", "R361"),
                                 ("candidate_caps_gcv_module", "candidate", "module")):
        sim, model = context(model_key, "glucose")
        medium = dict(model.medium)
        for rid in AMINO_ACID_ROWS:
            medium.pop(rid, None)
        model.medium = medium
        with model:
            applied = add_caps(model, caps, "glucose", gcv) if gcv else {}
            g, s = growth(model)
        table.append({"substrate": "glucose", "variant": f"RECONCILIATION_aa_free|{name}", "model": model_key,
                      "caps": applied, "growth": g, "status": s, "measured_growth": 0.2663,
                      "predicted_over_measured": None if g is None else g / 0.2663})
    return table, fluxes


def main():
    for key, (path, digest) in {**MODELS, **MEDIA, **INPUTS}.items():
        if digest:
            assert sha(path) == digest, f"Input identity changed: {key}"
    OUT.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    positives, _ = load_experimental_reference(INPUTS["experimental"][0])
    positives = set(positives)
    mtdna = json.loads(INPUTS["mtdna_scope"][0].read_text())["genes"]
    caps = json.loads(INPUTS["caps"][0].read_text())
    verify_cap_fractions(caps)
    record = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "command": "python -B model/reports/r132_mito_adh_20261012/run_task_b.py",
              "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
              "git_status": subprocess.check_output(["git", "status", "--short"], cwd=REPO, text=True),
              "inputs": {k: {"path": str(p), "sha256": sha(p)} for k, (p, _) in {**MODELS, **MEDIA, **INPUTS}.items()},
              "solver": PARAMS, "versions": {"python": platform.python_version(), "cobra": cobra.__version__,
                                             "gurobi": list(gurobipy.gurobi.version())},
              "scoring": "mtDNA genes out of scope (benchmark.md 2026-10-08); negatives are model genes absent from the "
                         "positive list (user-defined class); benchmark rule: no solution = unresolved; second view "
                         "counts no solution as lethal (labelled); development reference, not independent validation"}
    screens = {}
    for key in MODELS:
        result, rows = screen(key, positives, mtdna)
        screens[key] = rows
        record[key] = {k: v for k, v in result.items() if k != "route"}
        record[key]["route"] = {label: {"growth": r["growth"], "status": r["status"], "fva": r["fva"]}
                                for label, r in result["route"].items()}
        (OUT / f"pfba_route_{key}.json").open("x").write(
            json.dumps({label: r["pfba"] for label, r in result["route"].items()}) + "\n")
        with (OUT / f"gene_screen_{key}.tsv").open("x") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        primary = next(s for s in record[key]["scores"] if s["cutoff"] == 0.15 and s["view"] == "benchmark rule")
        print(key, "WT", record[key]["WT"], {k: v["ratio"] for k, v in record[key]["reaction_ko"].items()},
              {k: primary[k] for k in ("TP", "FP", "FN", "TN", "unresolved", "mcc")}, flush=True)
    base = {r["gene_id"]: r for r in screens["baseline"]}
    record["gene_differences"] = [
        {"gene_id": r["gene_id"], "baseline_ratio": base[r["gene_id"]]["ratio"], "candidate_ratio": r["ratio"],
         "baseline_status": base[r["gene_id"]]["status"], "candidate_status": r["status"]}
        for r in screens["candidate"]
        if base[r["gene_id"]]["status"] != r["status"]
        or (r["ratio"] is None) != (base[r["gene_id"]]["ratio"] is None)
        or (r["ratio"] is not None and abs(r["ratio"] - base[r["gene_id"]]["ratio"]) > 1e-6)]
    record["call_changes"] = {
        str(c): [g["gene_id"] for g in record["gene_differences"]
                 if (g["baseline_ratio"] is not None and g["baseline_ratio"] < c) != (g["candidate_ratio"] is not None and g["candidate_ratio"] < c)
                 or (g["baseline_ratio"] is None) != (g["candidate_ratio"] is None)]
        for c in CUTOFFS}
    table, fluxes = step5(caps)
    record["step5"] = table
    (OUT / "pfba_step5.json").open("x").write(json.dumps(fluxes) + "\n")
    for row in table:
        print("step5", row["substrate"], row["variant"], row["status"], row["growth"],
              None if row["predicted_over_measured"] is None else round(row["predicted_over_measured"], 4), flush=True)
    record.update(lps=Counter_.lps, elapsed_seconds=time.monotonic() - started,
                  inputs_unchanged=all(sha(p) == d for p, d in {**MODELS, **MEDIA, **INPUTS}.values() if d),
                  outputs={p.name: sha(p) for p in sorted(OUT.iterdir())})
    with (HERE / "run.json").open("x") as f:
        f.write(json.dumps(record, indent=1, default=str) + "\n")
    print("LPs", Counter_.lps, "elapsed", round(record["elapsed_seconds"]), "s")


if __name__ == "__main__":
    main()
