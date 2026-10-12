"""One bounded static screen. Run --self-check for a solver-free scoring check."""
import csv
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "platform"))

import cobra
import gurobipy
import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.solver import linear_reaction_coefficients
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL
from tools.dfba_new_fn_essentiality import load_experimental_reference
from tools.diagnose_closed_energy import configure_solver

HERE = Path(__file__).resolve().parent
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0,
              FeasibilityTol=1e-7, OptimalityTol=1e-7, IntFeasTol=1e-7)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def verdict(row, wt, cutoff):
    if not row["scorable"]:
        return "unresolved"
    return "essential" if row["raw_growth"] / wt < cutoff else "nonessential"


def metrics(rows, wt, cutoff):
    counts = Counter(TP=0, FN=0, FP=0, TN=0, unresolved=0)
    for row in rows:
        prediction = verdict(row, wt, cutoff)
        key = "unresolved" if prediction == "unresolved" else (
            ("TP" if row["experimental_positive"] else "FP")
            if prediction == "essential" else
            ("FN" if row["experimental_positive"] else "TN"))
        counts[key] += 1
    def ratio(a, b):
        return a / b if b else None
    return dict(cutoff=cutoff, **counts,
                recall=ratio(counts["TP"], counts["TP"] + counts["FN"]),
                precision=ratio(counts["TP"], counts["TP"] + counts["FP"]),
                specificity=ratio(counts["TN"], counts["TN"] + counts["FP"]))


def self_check():
    rows = [dict(scorable=True, raw_growth=g, experimental_positive=p)
            for g, p in [(0.0, True), (0.3, True), (0.2, False), (2.0, False)]]
    rows += [dict(scorable=False, raw_growth=0.0, experimental_positive=True)]
    result = metrics(rows, 2.0, 0.15)
    assert all(result[k] == 1 for k in ("TP", "FN", "FP", "TN", "unresolved"))
    assert verdict(rows[1], 2.0, 0.15) == "nonessential"  # Exact threshold.
    assert verdict(rows[-1], 2.0, 0.15) == "unresolved"  # Failure is not death.
    print("Scoring self-check passed; no optimization.")


def main():
    research = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
    assert all((research / x).is_dir() for x in ("reference", "state", "artifacts"))
    inputs = {
        "model": MODEL.reports / "coq_r305_candidate_20261005/E5_coq9_alpha_1e-4_R305_qcycle.xml",
        "medium": research / "state/media/sd_leu.csv",
        "profile": research / "state/strain_profiles/po1f_sd_leu.json",
        "experimental": research / "state/essentiality/repository/consensus_essential_genes.csv",
    }
    expected = dict(
        model="33468b94f378f0ab7f29f05e4356fd93b6d3eff23bc112d566cf16d6e4d6041c",
        medium="ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85",
        profile="35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383",
        experimental="1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b")
    assert {k: sha(p) for k, p in inputs.items()} == expected, "Input identity changed"
    positives, reference_metadata = load_experimental_reference(inputs["experimental"])
    reference = {r["gene_id"]: r for r in csv.DictReader(inputs["experimental"].open())}
    sim = load_effective_simulation_context(model_path=inputs["model"],
          media_path=inputs["medium"], strain_profile_path=inputs["profile"])
    model = sim.model
    model.solver = "gurobi"
    configure_solver(model, {"solver": PARAMS})
    assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {"biomass_C": 1.0}
    assert model.objective.direction == "max"
    assert model.reactions.biomass_C.get_coefficient("m468[C_mi]") == -1e-4
    assert model.reactions.R305.get_coefficient("m28[C_mi]") == -2
    assert model.reactions.R305.get_coefficient("m10[C_cy]") == 4
    assert model.reactions.xMAINTENANCE.lower_bound == 7.8625
    genes = sorted(g.id for g in model.genes if g.id not in sim.excluded_runtime_genes)
    assert len(genes) == 1073 and len(positives) == 1612
    code_paths = {Path(__file__).resolve(), MODEL.root / "expected/benchmark.md"}
    for module in list(sys.modules.values()):
        path = getattr(module, "__file__", None)
        if path and Path(path).resolve().is_relative_to(REPO / "platform"):
            code_paths.add(Path(path).resolve())
    code_hashes = {str(p): sha(p) for p in sorted(code_paths)}
    output = research / "artifacts/tasks" / HERE.name
    assert not (HERE / "manifest.json").exists(), "Never overwrite an existing screen"
    output.mkdir(parents=True, exist_ok=False)
    git = lambda *args: subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()
    manifest = dict(status="running", started_utc=datetime.now(timezone.utc).isoformat(),
        command=".venv/bin/python -B model/reports/screen_coq_r305_20261008/run_screen.py",
        branch=git("branch", "--show-current"), head=git("rev-parse", "HEAD"),
        git_status_before=git("status", "--short"),
        inputs={k: dict(path=str(p), sha256=expected[k]) for k, p in inputs.items()},
        code_sha256=code_hashes, reference_model=None, outputs=str(output),
        versions=dict(python=platform.python_version(), cobra=cobra.__version__,
                      gurobi=list(gurobipy.gurobi.version())),
        solver=PARAMS, call_limit=1074, wall_limit_seconds=900,
        simulation_context=sim.provenance(), active_medium=sim.active_medium,
        strain_overlay=sim.strain_overlay_audit,
        excluded_runtime_genes=list(sim.excluded_runtime_genes),
        reference_metadata=reference_metadata,
        scoring_negative_class="user_defined_native_model_genes_absent_from_positive_list",
        validation_status="development_reference_not_independent_validation",
        gene_count=len(genes), experimental_positive_count=len(positives),
        in_model_positive_count=len(set(positives) & set(genes)),
        out_of_scope_positive_count=len(set(positives) - set(genes)),
        negative_growth_policy="unresolved_no_clipping", numerical_tolerance=1e-6)
    write_json(HERE / "manifest.json", manifest)
    baseline = [(r.id, r.bounds, r.gene_reaction_rule) for r in model.reactions]
    write_json(output / "effective_model.json", dict(
        reactions=[dict(id=r.id, bounds=r.bounds, gpr=r.gene_reaction_rule,
                        metabolites={m.id: c for m, c in r.metabolites.items()}) for r in model.reactions],
        objective={"biomass_C": 1.0}, direction="max"))
    with (output / "experimental_coverage.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["gene_id", "source_gene_id", "coverage", "source_putative_function"])
        for gene in positives:
            writer.writerow([gene, reference[gene]["source_gene_id"],
                "in_native_model" if gene in genes else "absent_direct_native_id_no_cross_version_mapping_attempted",
                reference[gene]["function"]])
    matrix = create_stoichiometric_matrix(model, array_type="lil").tocsr()
    started = time.monotonic()
    rows, calls = [], 0
    journal = (output / "solver_calls.jsonl").open("w", buffering=1)

    def solve(label):
        nonlocal calls
        if calls >= 1074 or time.monotonic() - started > 840:
            raise RuntimeError("Declared screening budget reached")
        calls += 1
        journal.write(json.dumps(dict(event="start", call=calls, label=label)) + "\n")
        t = time.monotonic()
        with warnings.catch_warnings(record=True) as warning_list:
            solution = model.optimize()
        raw = float(solution.objective_value)
        flux = solution.fluxes.to_numpy(dtype=float)
        finite = bool(np.isfinite(flux).all())
        residual = float(np.max(np.abs(matrix @ flux))) if finite else None
        bounds = np.array([r.bounds for r in model.reactions])
        violation = float(max(0.0, np.max(bounds[:, 0] - flux), np.max(flux - bounds[:, 1]))) if finite else None
        row = dict(label=label, solver_status=solution.status,
            native_status=int(model.solver.problem.Status),
            raw_growth=raw if math.isfinite(raw) else None, raw_growth_repr=repr(raw),
            seconds=time.monotonic() - t, max_mass_residual=residual,
            max_bound_violation=violation,
            warnings=[str(w.message) for w in warning_list],
            scorable=bool(solution.status == "optimal" and math.isfinite(raw) and raw >= 0
                          and finite and residual <= 1e-6 and violation <= 1e-6))
        journal.write(json.dumps(dict(event="finish", call=calls, **row), allow_nan=False) + "\n")
        if label == "WT":
            solution.fluxes.to_csv(output / "WT_fluxes.tsv", sep="\t", header=["flux"])
        return row

    try:
        wt = solve("WT")
        manifest["WT"] = wt
        assert wt["scorable"] and wt["raw_growth"] > 0, "Invalid WT; stop without KO"
        print(f"WT={wt['raw_growth']:.12g}; screening {len(genes)} native genes", flush=True)
        fields = ["gene_id", "model_name_unverified", "verified_symbol", "source_putative_function",
                  "function_evidence", "model_reaction_ids", "experimental_positive",
                  "runtime_confounded", "solver_status", "native_status", "raw_growth", "raw_growth_repr",
                  "ko_wt_ratio", "scorable", "max_mass_residual", "max_bound_violation", "seconds",
                  "changed_reaction_bounds", *[f"prediction_{int(c*100)}pct" for c in CUTOFFS]]
        with (output / "gene_screen.tsv").open("w", buffering=1) as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            for index, gene_id in enumerate(genes, 1):
                gene = model.genes.get_by_id(gene_id)
                before = {r.id: r.bounds for r in gene.reactions}
                with model:
                    gene.knock_out()
                    changes = {r.id: [before[r.id], r.bounds] for r in gene.reactions if before[r.id] != r.bounds}
                    row = solve(gene_id)
                    row.update(gene_id=gene_id, model_name_unverified=gene.name,
                        verified_symbol="not verified in this screen",
                        source_putative_function=reference.get(gene_id, {}).get("function", "not verified; see model_reaction_ids"),
                        function_evidence="source putative annotation and/or model assignment; not catalytic validation",
                        model_reaction_ids=";".join(sorted(before)), experimental_positive=gene_id in reference,
                        runtime_confounded=any(op.get("gene_id") == gene_id for op in sim.strain_profile["operations"]),
                        ko_wt_ratio=row["raw_growth"] / wt["raw_growth"] if row["raw_growth"] is not None else None,
                        changed_reaction_bounds=json.dumps(changes, sort_keys=True))
                    row.update({f"prediction_{int(c*100)}pct": verdict(row, wt["raw_growth"], c) for c in CUTOFFS})
                    rows.append(row)
                    writer.writerow(row)
                assert gene.functional and all(model.reactions.get_by_id(r).bounds == b for r, b in before.items())
                if index % 100 == 0 or index == len(genes):
                    print(f"{index}/{len(genes)} completed; statuses={dict(Counter(r['solver_status'] for r in rows))}", flush=True)
        assert baseline == [(r.id, r.bounds, r.gene_reaction_rule) for r in model.reactions]
        assert {k: sha(p) for k, p in inputs.items()} == expected
        assert {p: sha(Path(p)) for p in code_hashes} == code_hashes
        manifest.update(status="complete", all_inputs_and_code_unchanged=True,
                        context_restored=True,
                        metrics=[metrics(rows, wt["raw_growth"], c) for c in CUTOFFS])
    except Exception as exc:
        manifest.update(status="stopped", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        journal.close()
        manifest.update(finished_utc=datetime.now(timezone.utc).isoformat(), calls_used=calls,
            genes_completed=len(rows), elapsed_seconds=time.monotonic() - started,
            status_counts=dict(Counter(r["solver_status"] for r in rows)),
            unresolved_genes=[r["gene_id"] for r in rows if not r["scorable"]],
            output_sha256={str(p): sha(p) for p in sorted(output.iterdir()) if p.is_file()})
        write_json(HERE / "manifest.json", manifest)
    print(json.dumps(manifest["metrics"], indent=2), flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["--self-check"]:
        self_check()
    elif not sys.argv[1:]:
        main()
    else:
        raise SystemExit("Only --self-check or no arguments are supported")
