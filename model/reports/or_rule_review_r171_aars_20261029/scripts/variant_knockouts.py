"""Knockout mechanics of the R171 ATP synthase rule and of the cytosolic aminoacyl-tRNA synthetase rules
(R759 with R760, R723, R673) on candidate 965ab947, plus in-memory rule variants scored at the four
benchmark cutoffs. Read-only: no repository file is changed; rules are edited in memory and restored.

Setup follows reports/essentiality_fn_causes_20261027/run_screen.py: Task A medium + PO1f profile, Gurobi at
FeasibilityTol / OptimalityTol 1e-9, Presolve 0, one thread, basis reset before every solve, and a point is
scorable only if the solve is optimal, |S v| and bound violations are <= 1e-6 and no knocked reaction carries
flux above 1e-9. A failed or infeasible solve is "unresolved", never "dead".

Usage (repository root, IYALI26_RESEARCH_ROOT set):
    <python> -B model/reports/or_rule_review_r171_aars_20261029/scripts/variant_knockouts.py
"""
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TASK = HERE.parent
REPO = TASK.parents[2]
sys.path.insert(0, str(REPO / "platform"))

import cobra  # noqa: E402
import gurobipy  # noqa: E402
from cobra.util.array import create_stoichiometric_matrix  # noqa: E402
from cobra.util.solver import linear_reaction_coefficients  # noqa: E402
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context  # noqa: E402
from scripts.gem_annotate.model_layout import MODEL  # noqa: E402

CUTOFFS = (0.01, 0.05, 0.10, 0.15)
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0, FeasibilityTol=1e-9, OptimalityTol=1e-9, IntFeasTol=1e-9)
KNOCKED_TOL, GENERAL_TOL = 1e-9, 1e-6
LP_BUDGET, WALL_BUDGET = 1500, 1800
MODEL_NAME = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
              "cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_r2171_1to1_20261028.xml")
EXPECTED = dict(
    model="965ab947f2aadc2c69443e7a794358c502acfdf0892e78582278ff07a37429d7",
    medium="6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45",
    profile="35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383",
    experimental="1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b",
    mtdna="16b3e9b4418f24821b8775936fd377ae21f77afa8e8af2b8949a246f62c44996")

AND14 = ["YALI1M00260g", "YALI1M00267g", "YALI1M00318g", "YALI1F04433g", "YALI1F04110g", "YALI1D27826g",
         "YALI1D15734g", "YALI1D14639g", "YALI1B08885g", "YALI1F07384g", "YALI1B05364g", "YALI0D17490g",
         "YALI1F26967g", "YALI1E38131g"]
G, K, E, J = "YALI1B28025g", "YALI1B15733g", "YALI1E38131g", "YALI0D17490g"  # subunits g, k, e, j


def and_rule(genes):
    return "(" + " and ".join(genes) + ")"


ORIGINAL = {
    "R171": f"{and_rule(AND14)} or {G} or {K}",
    "R759": "YALI1A09121g or YALI1F14882g",
    "R760": "YALI1A09121g or YALI1F14882g",
    "R723": "YALI1B07957g or YALI1F03908g",
    "R673": "YALI1D19110g or YALI1E05929g",
}
AARS_SPLIT = {"R759": "YALI1F14882g", "R760": "YALI1A09121g", "R723": "YALI1F03908g", "R673": "YALI1E05929g"}
VARIANTS = {
    "baseline": {},
    "R171_A_drop_g_k": {"R171": and_rule(AND14)},
    "R171_B_and16": {"R171": and_rule(AND14 + [G, K])},
    "R171_C_drop_e_g_k": {"R171": and_rule([x for x in AND14 if x != E])},
    "R171_D_drop_e_j_g_k": {"R171": and_rule([x for x in AND14 if x not in (E, J)])},
    "aaRS_split": dict(AARS_SPLIT),
    "aaRS_split_alt_R673_D19110": {**AARS_SPLIT, "R673": "YALI1D19110g"},
    "combined_A_split": {"R171": and_rule(AND14), **AARS_SPLIT},
    "combined_D_split": {"R171": and_rule([x for x in AND14 if x not in (E, J)]), **AARS_SPLIT},
}
FULL_SWEEP_VARIANT = "combined_A_split"
REACTION_KOS = {"R171": ["R171"], "R759": ["R759"], "R760": ["R760"], "R759+R760": ["R759", "R760"],
                "R723": ["R723"], "R673": ["R673"]}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verdict(row, wt, cutoff):
    if not row["scorable"]:
        return "unresolved"
    return "essential" if row["growth"] / wt < cutoff else "nonessential"


def main():
    started = time.monotonic()
    research = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
    assert all((research / x).is_dir() for x in ("reference", "state", "artifacts"))
    inputs = {
        "model": REPO / "model/candidates" / MODEL_NAME,
        "medium": MODEL.root / "conditions/media/po1f_csm_leu_w29uptake_20261012.csv",
        "profile": research / "state/strain_profiles/po1f_sd_leu.json",
        "experimental": research / "state/essentiality/repository/consensus_essential_genes.csv",
        "mtdna": MODEL.root / "expected/mtdna_genes_out_of_scope.json",
    }
    assert {k: sha(p) for k, p in inputs.items()} == EXPECTED, "Input identity changed"
    # Short-path copy (macOS cobra rejects str paths of 260+ characters); byte identity checked.
    scratch = Path(os.environ.get("OR_RULE_SCRATCH", "/tmp")) / "or_rule_review_model_965ab947.xml"
    scratch.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(inputs["model"], scratch)
    assert sha(scratch) == EXPECTED["model"]
    positives = {r["gene_id"] for r in csv.DictReader(inputs["experimental"].open())}
    mtdna = set(json.loads(inputs["mtdna"].read_text())["genes"])
    assert len(positives) == 1612

    sim = load_effective_simulation_context(model_path=scratch, media_path=inputs["medium"],
                                            strain_profile_path=inputs["profile"])
    model = sim.model
    model.solver = "gurobi"
    for key, value in PARAMS.items():
        model.solver.problem.setParam(key, value)
    assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {"biomass_C": 1.0}
    assert model.reactions.xMAINTENANCE.lower_bound == 7.8625
    for rid, rule in ORIGINAL.items():
        assert model.reactions.get_by_id(rid).gene_reaction_rule == rule, (rid, model.reactions.get_by_id(rid).gene_reaction_rule)
    genes_all = sorted(g.id for g in model.genes if g.id not in sim.excluded_runtime_genes)
    assert len(genes_all) == 1081
    target_genes = sorted({g for rule in ORIGINAL.values() for g in rule.replace("(", " ").replace(")", " ").split() if g.startswith("YALI")})
    assert len(target_genes) == 22
    S = create_stoichiometric_matrix(model, array_type="dense")
    baseline_bounds = {r.id: r.bounds for r in model.reactions}
    calls = 0
    journal = (TASK / "solver_calls.jsonl").open("w", buffering=1)

    def solve(label, knocked=()):
        nonlocal calls
        if calls >= LP_BUDGET or time.monotonic() - started > WALL_BUDGET:
            raise RuntimeError("Declared budget reached")
        calls += 1
        model.solver.problem.reset()  # discard the previous basis: a knockout never starts from the WT optimum
        with warnings.catch_warnings(record=True) as caught:
            solution = model.optimize()
        raw = solution.objective_value
        ok = solution.status == "optimal" and raw is not None and math.isfinite(raw)
        row = dict(label=label, status=solution.status, growth=float(raw) if ok else None,
                   max_residual=None, max_bound_violation=None, max_knocked_flux=None, scorable=False,
                   warnings=[str(w.message) for w in caught])
        if ok:
            flux = solution.fluxes.to_numpy(dtype=float)
            bounds = np.array([r.bounds for r in model.reactions])
            row["max_residual"] = float(np.max(np.abs(S @ flux)))
            row["max_bound_violation"] = float(max(0.0, np.max(bounds[:, 0] - flux), np.max(flux - bounds[:, 1])))
            row["max_knocked_flux"] = float(max((abs(solution.fluxes[r]) for r in knocked), default=0.0))
            row["scorable"] = bool(raw >= 0 and row["max_residual"] <= GENERAL_TOL
                                   and row["max_bound_violation"] <= GENERAL_TOL and row["max_knocked_flux"] <= KNOCKED_TOL)
        journal.write(json.dumps(dict(call=calls, **row)) + "\n")
        return row

    def knockout(gene_id):
        gene = model.genes.get_by_id(gene_id)
        original_functional = gene.functional
        before = {r.id: r.bounds for r in gene.reactions}
        with model:
            gene.knock_out()
            closed = tuple(r.id for r in gene.reactions if before[r.id] != r.bounds)
            row = solve(gene_id, knocked=closed)
        row["closed_reactions"] = list(closed)
        if gene._model is None:  # SBML orphan genes lack a model link, so the context cannot undo them (as in run_screen.py)
            gene.functional = original_functional
        assert gene.functional == original_functional and all(model.reactions.get_by_id(r).bounds == b for r, b in before.items())
        return row

    def apply_rules(rules):
        for rid, rule in rules.items():
            model.reactions.get_by_id(rid).gene_reaction_rule = rule

    def restore_rules():
        apply_rules(ORIGINAL)
        for rid, rule in ORIGINAL.items():
            assert model.reactions.get_by_id(rid).gene_reaction_rule == rule
        assert {r.id: r.bounds for r in model.reactions} == baseline_bounds

    results = dict(status="running", started_utc=datetime.now(timezone.utc).isoformat(),
                   inputs={k: dict(path=str(p), sha256=EXPECTED[k]) for k, p in inputs.items()},
                   script_sha256=sha(__file__), solver=PARAMS,
                   versions=dict(python=platform.python_version(), cobra=cobra.__version__, gurobi=list(gurobipy.gurobi.version())),
                   git=dict(head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
                            branch=subprocess.check_output(["git", "branch", "--show-current"], cwd=REPO, text=True).strip()),
                   simulation_context=sim.provenance(), excluded_runtime_genes=sorted(sim.excluded_runtime_genes),
                   original_rules=ORIGINAL, variants=VARIANTS, target_genes=target_genes)
    try:
        wt = solve("WT")
        assert wt["scorable"] and wt["growth"] > 0
        results["WT"] = wt
        print(f"WT {wt['growth']:.10g}", flush=True)

        # Reaction-level knockouts on the original rules.
        reaction_level = {}
        for label, rids in REACTION_KOS.items():
            with model:
                for rid in rids:
                    model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                reaction_level[label] = solve(f"rxn:{label}", knocked=tuple(rids))
            reaction_level[label]["ratio"] = (reaction_level[label]["growth"] / wt["growth"]) if reaction_level[label]["growth"] is not None else None
        results["reaction_level_knockouts"] = reaction_level
        print("reaction-level:", {k: v["ratio"] for k, v in reaction_level.items()}, flush=True)

        # Gene knockouts of the 22 target genes under every variant.
        per_variant = {}
        for name, rules in VARIANTS.items():
            apply_rules(rules)
            rows = {}
            for gene_id in target_genes:
                row = knockout(gene_id)
                row["ratio"] = row["growth"] / wt["growth"] if row["growth"] is not None else None
                row["classes"] = {str(c): verdict(row, wt["growth"], c) for c in CUTOFFS}
                rows[gene_id] = row
            restore_rules()
            per_variant[name] = rows
            print(f"variant {name}: {sum(1 for r in rows.values() if r['classes']['0.15'] == 'essential')} essential at 15 %", flush=True)
        results["gene_knockouts_by_variant"] = per_variant

        # Class transitions versus the in-run baseline, per variant and cutoff (mtDNA genes out of scope).
        def cls(gene_id, pred):
            if gene_id in mtdna:
                return "mtDNA_out_of_scope"
            if pred == "unresolved":
                return "unresolved"
            pos = gene_id in positives
            return ("TP" if pos else "FP") if pred == "essential" else ("FN" if pos else "TN")
        transitions = {}
        for name in VARIANTS:
            transitions[name] = {}
            for c in CUTOFFS:
                counter = Counter()
                for gene_id in target_genes:
                    b = cls(gene_id, per_variant["baseline"][gene_id]["classes"][str(c)])
                    a = cls(gene_id, per_variant[name][gene_id]["classes"][str(c)])
                    if a != b:
                        counter[f"{b}->{a}"] += 1
                transitions[name][str(c)] = dict(counter)
        results["class_transitions_vs_baseline"] = transitions

        # Confirming full sweep of one variant.
        apply_rules(VARIANTS[FULL_SWEEP_VARIANT])
        sweep_dir = research / "artifacts/tasks" / TASK.name
        sweep_dir.mkdir(parents=True, exist_ok=True)
        sweep_rows = []
        with (sweep_dir / f"gene_screen_{FULL_SWEEP_VARIANT}.tsv").open("w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["gene_id", "experimental_positive", "mtdna_out_of_scope", "status", "growth", "ratio",
                             "scorable", "max_knocked_flux", "closed_reactions", *[f"class_{int(c*100)}pct" for c in CUTOFFS]])
            for index, gene_id in enumerate(genes_all, 1):
                row = knockout(gene_id)
                row["ratio"] = row["growth"] / wt["growth"] if row["growth"] is not None else None
                classes = [cls(gene_id, verdict(row, wt["growth"], c)) for c in CUTOFFS]
                sweep_rows.append((gene_id, classes, row))
                writer.writerow([gene_id, gene_id in positives, gene_id in mtdna, row["status"], row["growth"], row["ratio"],
                                 row["scorable"], row["max_knocked_flux"], ";".join(row["closed_reactions"]), *classes])
                if index % 250 == 0:
                    print(f"sweep {index}/{len(genes_all)}", flush=True)
        restore_rules()
        sweep_metrics = {}
        for i, c in enumerate(CUTOFFS):
            counter = Counter(classes[i] for _, classes, _ in sweep_rows)
            tp, fn, fp, tn = (counter[k] for k in ("TP", "FN", "FP", "TN"))
            d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
            sweep_metrics[str(c)] = dict(TP=tp, FN=fn, FP=fp, TN=tn, unresolved=counter["unresolved"],
                                         mtDNA_out_of_scope=counter["mtDNA_out_of_scope"],
                                         recall=tp / (tp + fn), precision=tp / (tp + fp) if tp + fp else None,
                                         mcc=(tp * tn - fp * fn) / d if d else None)
        results["full_sweep"] = dict(variant=FULL_SWEEP_VARIANT, metrics=sweep_metrics,
                                     unresolved_genes=[g for g, cl, _ in sweep_rows if "unresolved" in cl],
                                     failing_points=[g for g, _, r in sweep_rows if r["status"] == "optimal" and not r["scorable"]],
                                     output=str(sweep_dir / f"gene_screen_{FULL_SWEEP_VARIANT}.tsv"))
        results.update(status="complete")
    except Exception as exc:  # keep the partial record
        results.update(status="stopped", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        journal.close()
        results.update(finished_utc=datetime.now(timezone.utc).isoformat(), lp_calls=calls,
                       elapsed_seconds=time.monotonic() - started)
        (TASK / "variant_knockouts.json").write_text(json.dumps(results, indent=1) + "\n")
    print(json.dumps(results["full_sweep"]["metrics"]["0.15"], indent=1))
    print(json.dumps(results["class_transitions_vs_baseline"], indent=1))


if __name__ == "__main__":
    main()
