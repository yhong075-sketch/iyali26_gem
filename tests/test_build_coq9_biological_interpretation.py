import hashlib
import json
import math
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from scripts.gem_annotate import build_coq9_biological_interpretation as bundle
from scripts.gem_annotate.coq9_dilution import runtime_coq9_dilution_record
from scripts.gem_annotate.essentiality_simulation_context import sha256_file, sha256_payload


REPO = Path(__file__).parents[1]
SPEC = replace(
    bundle.FORMAL_SPEC,
    panel=("YALI1C26017g", "YALI1F08349g"),
    alphas=(1e-4,), pools=(0.0, 1.0), chunk_count=2,
    hours=0.125, enforce_run_names=False,
)
JOB_IDS = {
    "finite_batch": {"array": "1001", "merge": "1002"},
    "po1f_nonlimiting": {"array": "1003", "merge": "1004"},
}


def _compute_sha(path: str) -> str:
    value = __import__("subprocess").check_output(
        ["git", "show", f"{bundle.COMPUTE_COMMIT}:{path}"], cwd=REPO,
    )
    return hashlib.sha256(value).hexdigest()


def _trace(entity: str, mode: str, alpha: float, pool: float) -> list[dict]:
    biomass = SPEC.initial_biomass
    glucose = 111.0
    uracil = 0.178428
    q9 = alpha * biomass * pool
    rows = []
    for step in range(2):
        growth = 0.1
        glucose_flux = 0.2
        uracil_flux = 0.01
        biomass_end = biomass * (1 + growth * SPEC.dt)
        glucose_end = glucose - glucose_flux * biomass * SPEC.dt
        uracil_end = uracil - uracil_flux * biomass * SPEC.dt
        row = {
            "gene_id": entity, "uracil_mode": mode, "step_index": step,
            "time_h": step * SPEC.dt, "time_end_h": (step + 1) * SPEC.dt,
            "interval_advanced": True,
            "biomass_gDW_L": biomass, "biomass_end_gDW_L": biomass_end,
            "growth_h-1": growth, "biomass_flux_h-1": growth,
            "objective_value": growth,
            "source_free_solver_status": "optimal", "source_free_growth_h-1": growth,
            "q9_pool_mmol_L": q9, "q9_pool_end_mmol_L": q9,
            "q9_source_flux_mmol_gDW_h": 0.0,
            "coq9_dilution_flux_mmol_gDW_h": alpha * growth,
            "uracil_mmol_L": uracil if mode == "finite_batch" else math.nan,
            "uracil_end_mmol_L": uracil_end if mode == "finite_batch" else math.nan,
            "uracil_uptake_flux_mmol_gDW_h": uracil_flux,
            "glucose_mmol_L": glucose, "glucose_end_mmol_L": glucose_end,
            "glucose_uptake_flux_mmol_gDW_h": glucose_flux,
            "oxygen_uptake_flux_mmol_gDW_h": 0.3,
            "atp_maintenance_flux_mmol_gDW_h": 7.8625,
            "status": "optimal",
            **{column: 0.0 for column in bundle.REACTION_FLUX_COLUMNS},
            "alpha_mmol_gDW": alpha, "alpha_replicate_id": math.nan,
            "pool_multiplier": pool,
        }
        row["reaction_flux_R1287_mmol_gDW_h"] = -0.3
        row["reaction_flux_xMAINTENANCE_mmol_gDW_h"] = 7.8625
        rows.append(row)
        biomass, glucose, uracil = biomass_end, glucose_end, uracil_end
    return rows


def _call(gene_id: str, mode: str, alpha: float, pool: float) -> dict:
    trace = _trace(gene_id, mode, alpha, pool)
    final = trace[-1]
    doublings = math.log2(final["biomass_end_gDW_L"] / SPEC.initial_biomass)
    return {
        "gene_id": gene_id, "uracil_mode": mode,
        "final_biomass_gDW_L": final["biomass_end_gDW_L"],
        "dynamic_doublings": doublings, "initial_growth_h-1": 0.1,
        "initial_source_free_solver_status": "optimal",
        "initial_source_free_growth_h-1": 0.1,
        "q9_pool_depleted_h": 0.0 if pool == 0 else math.nan,
        "q9_source_total_mmol_L": 0.0,
        "final_glucose_mmol_L": final["glucose_end_mmol_L"],
        "final_uracil_mmol_L": (
            final["uracil_end_mmol_L"] if mode == "finite_batch" else math.nan
        ),
        "termination_status": "completed", "termination_time_h": SPEC.hours,
        "alpha_mmol_gDW": alpha, "alpha_replicate_id": math.nan,
        "pool_multiplier": pool, "runtime_gpr_scenario": "",
        "runtime_gpr_mapping_sha256": "", "wt_dynamic_doublings": doublings,
        "dynamic_growth_ratio": 1.0, "wt_termination_status": "completed",
        "wt_termination_time_h": SPEC.hours, "experimental_essential": False,
        **{column: False for column in bundle.CUTOFF_COLUMNS},
    }


def _make_run(root: Path, mode: str) -> Path:
    run = root / f"fixture_{mode}"
    run.mkdir()
    fingerprints = []
    calls_paths = []
    for index, gene_id in enumerate(SPEC.panel):
        calls = pd.DataFrame([
            _call(gene_id, mode, alpha, pool)
            for alpha in SPEC.alphas for pool in SPEC.pools
        ])
        trajectories = pd.DataFrame([
            row
            for alpha in SPEC.alphas for pool in SPEC.pools
            for entity in ("WT", gene_id)
            for row in _trace(entity, mode, alpha, pool)
        ], columns=bundle.TRAJECTORY_COLUMNS)
        calls_path = run / f"chunk_{index:03d}_calls.tsv"
        trajectory_path = run / f"chunk_{index:03d}_trajectory.tsv"
        calls.to_csv(calls_path, sep="\t", index=False)
        trajectories.to_csv(trajectory_path, sep="\t", index=False)
        calls_paths.append(calls_path)
        manifest = {
            "workflow": "quinone_dfba_essentiality",
            "schema_version": bundle.SCHEMA_VERSION,
            "run_id": run.name,
            "trajectory_semantics": bundle.TRAJECTORY_SEMANTICS,
            "chunk_index": index, "chunk_count": SPEC.chunk_count,
            "genes": [gene_id], "solver": "gurobi",
            "runtime_versions": {
                "python": "3.11", "cobra": "0.30.0", "gurobipy": "10.0.3",
            },
            "solver_feasibility_tolerance": bundle.SOLVER_TOL,
            "script_sha256": _compute_sha("scripts/gem_annotate/quinone_dfba_essentiality.py"),
            "coq9_dilution_tool_sha256": _compute_sha("scripts/gem_annotate/coq9_dilution.py"),
            "hours": SPEC.hours, "dt_h": SPEC.dt,
            "initial_biomass_gDW_L": SPEC.initial_biomass,
            "alphas_mmol_gDW": list(SPEC.alphas),
            "alpha_sampling": bundle._explicit_alpha_sampling(SPEC),
            "pool_multipliers": list(SPEC.pools), "uracil_mode": mode,
            "initial_pools_mmol_L": {
                **bundle.INITIAL_POOLS_MMOL_L,
                "R1354": None if mode == "po1f_nonlimiting" else 0.178428,
            },
            "base_r1354_bound_mmol_gDW_h": 1000.0,
            "q9_reserve_definition": "alpha * initial_biomass * pool_multiplier mmol/L",
            "q9_reserve_policy": bundle.Q9_RESERVE_POLICY,
            "optimizer": "pfba", "nonoptimal_policy": "terminal_record_and_stop",
            "calibration_status": "sensitivity_only_not_calibrated",
            "coq9_dilution": [
                runtime_coq9_dilution_record(alpha, pool_source_enabled=True)
                for alpha in SPEC.alphas
            ],
            "runtime_topology": {"enabled": False},
            "runtime_gpr_scenario": {"enabled": False},
            "simulation_context": {
                "simulation_context_fingerprint_version": "1",
                "simulation_context_fingerprint": bundle.SIMULATION_CONTEXT_SHA256,
                "strain_overlay_enabled": True,
                "strain_profile_id": "po1f_sd_leu_accrispr_v1",
                "strain_profile_sha256": bundle.PROFILE_SHA256,
                "strain_overlay_effect_fingerprint_version": "1",
                "strain_overlay_effect_sha256": bundle.OVERLAY_EFFECT_SHA256,
            },
            "input_sha256": {
                "model": SPEC.model_sha256, "medium": bundle.MEDIUM_SHA256,
                "profile": bundle.PROFILE_SHA256,
                "experimental": bundle.EXPERIMENTAL_SHA256,
            },
            "output_files": {
                "calls": {"filename": calls_path.name, "sha256": sha256_file(calls_path)},
                "trajectory": {
                    "filename": trajectory_path.name,
                    "sha256": sha256_file(trajectory_path),
                },
            },
            "created_at_utc": f"2026-09-02T00:00:0{index}Z",
        }
        manifest["fingerprint"] = sha256_payload(manifest)
        (run / f"chunk_{index:03d}_manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8",
        )
        fingerprints.append(manifest["fingerprint"])
    pd.concat([pd.read_csv(path, sep="\t") for path in calls_paths], ignore_index=True).to_csv(
        run / "essentiality_dynamic_calls.tsv", sep="\t", index=False,
    )
    (run / "merge_manifest.json").write_text(json.dumps({
        "workflow": "quinone_dfba_essentiality",
        "schema_version": bundle.SCHEMA_VERSION,
        "trajectory_semantics": bundle.TRAJECTORY_SEMANTICS,
        "chunks": SPEC.chunk_count, "chunk_fingerprints": sorted(fingerprints),
        "calls": SPEC.calls_per_mode,
        "calibration_status": "sensitivity_only_not_calibrated",
        "alpha_sampling": bundle._explicit_alpha_sampling(SPEC),
        "coq9_dilution": [
            runtime_coq9_dilution_record(alpha, pool_source_enabled=True)
            for alpha in SPEC.alphas
        ],
        "coq9_dilution_tool_sha256": _compute_sha("scripts/gem_annotate/coq9_dilution.py"),
        "q9_reserve_policy": bundle.Q9_RESERVE_POLICY,
        "runtime_gpr_scenario": {"enabled": False},
    }), encoding="utf-8")
    return run


def _refresh_manifest_fingerprints(run: Path) -> None:
    manifests = [
        json.loads(path.read_text()) for path in sorted(run.glob("chunk_*_manifest.json"))
    ]
    merge_path = run / "merge_manifest.json"
    merge = json.loads(merge_path.read_text())
    merge["chunk_fingerprints"] = sorted(item["fingerprint"] for item in manifests)
    merge_path.write_text(json.dumps(merge), encoding="utf-8")


def _rehash_trajectory(run: Path, index: int) -> None:
    manifest_path = run / f"chunk_{index:03d}_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    trajectory_path = run / f"chunk_{index:03d}_trajectory.tsv"
    manifest["output_files"]["trajectory"]["sha256"] = sha256_file(trajectory_path)
    manifest["fingerprint"] = sha256_payload({
        key: value for key, value in manifest.items() if key != "fingerprint"
    })
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _refresh_manifest_fingerprints(run)


def test_small_dual_mode_fixture_builds_complete_bundle(tmp_path):
    finite = _make_run(tmp_path, "finite_batch")
    nonlimiting = _make_run(tmp_path, "po1f_nonlimiting")
    output, archive = bundle.build_bundle(
        finite, nonlimiting, tmp_path / "bundle",
        compute_commit=bundle.COMPUTE_COMMIT,
        data_sha256=bundle.DATA_SHA256,
        repo_root=REPO, spec=SPEC, job_ids=JOB_IDS,
    )

    expected = {
        "run_manifest.json", "selected_gene_panel.tsv", "reaction_gene_map.tsv",
        "gene_condition_summary.tsv", "stepwise_trajectories.tsv.gz",
        "flux_panel.tsv.gz", "experimental_comparison.tsv", "numerical_checks.json",
        "grid_summary.tsv", "alpha_sensitivity.tsv", "pool_monotonicity.tsv",
        "q9_source_use.tsv", "five_gene_theoretical_controls.tsv",
        "BIOLOGICAL_INTERPRETATION.md", "SHA256SUMS", "plots",
    }
    assert expected <= {path.name for path in output.iterdir()}
    assert archive.is_file()
    assert Path(f"{archive}.sha256").is_file()
    assert len(pd.read_csv(output / "gene_condition_summary.tsv", sep="\t")) == 8
    assert len(pd.read_csv(output / "stepwise_trajectories.tsv.gz", sep="\t")) == 24
    assert json.loads((output / "numerical_checks.json").read_text())["all_checks_pass"]
    run_manifest = json.loads((output / "run_manifest.json").read_text())
    assert run_manifest["runs"]["finite_batch"]["slurm_jobs"] == JOB_IDS["finite_batch"]
    mapping = run_manifest["runtime_mapping"]
    assert mapping["mapping_sha256"] == bundle.sha256_payload({
        key: value for key, value in mapping.items() if key != "mapping_sha256"
    })
    panel = pd.read_csv(output / "selected_gene_panel.tsv", sep="\t")
    assert panel[["name_or_symbol", "function", "evidence_status"]].notna().all().all()
    assert "reaction_name" in pd.read_csv(output / "reaction_gene_map.tsv", sep="\t").columns


def test_merged_calls_require_exact_frozen_serialization(tmp_path):
    run = _make_run(tmp_path, "finite_batch")
    chunk = run / "chunk_000_calls.tsv"
    frame = pd.read_csv(chunk, sep="\t")
    frame.loc[0, "initial_growth_h-1"] = 2.1814287324998794e-11
    frame.to_csv(chunk, sep="\t", index=False)
    chunks = [
        pd.read_csv(run / f"chunk_{index:03d}_calls.tsv", sep="\t")
        for index in range(SPEC.chunk_count)
    ]
    merged_path = run / "essentiality_dynamic_calls.tsv"
    pd.concat(chunks, ignore_index=True).to_csv(merged_path, sep="\t", index=False)
    reparsed = pd.read_csv(merged_path, sep="\t")
    assert chunks[0].loc[0, "initial_growth_h-1"] != reparsed.loc[0, "initial_growth_h-1"]

    manifests = [
        json.loads((run / f"chunk_{index:03d}_manifest.json").read_text())
        for index in range(SPEC.chunk_count)
    ]
    bundle._validate_calls(run, "finite_batch", manifests, SPEC)

    tampered = bytearray(merged_path.read_bytes())
    tampered[-2] = ord("1") if tampered[-2] != ord("1") else ord("2")
    merged_path.write_bytes(tampered)
    with pytest.raises(ValueError, match="merged calls differ"):
        bundle._validate_calls(run, "finite_batch", manifests, SPEC)


def test_valid_hash_cannot_hide_interval_conservation_failure(tmp_path):
    finite = _make_run(tmp_path, "finite_batch")
    nonlimiting = _make_run(tmp_path, "po1f_nonlimiting")
    path = finite / "chunk_000_trajectory.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame.loc[0, "glucose_end_mmol_L"] += 0.01
    frame.to_csv(path, sep="\t", index=False)
    _rehash_trajectory(finite, 0)

    with pytest.raises(ValueError, match="glucose balance"):
        bundle.build_bundle(
            finite, nonlimiting, tmp_path / "bundle",
            compute_commit=bundle.COMPUTE_COMMIT,
            data_sha256=bundle.DATA_SHA256,
            repo_root=REPO, spec=SPEC, job_ids=JOB_IDS,
        )


def test_incomplete_initial_inventory_is_rejected(tmp_path):
    finite = _make_run(tmp_path, "finite_batch")
    nonlimiting = _make_run(tmp_path, "po1f_nonlimiting")
    path = finite / "chunk_000_manifest.json"
    manifest = json.loads(path.read_text())
    manifest["initial_pools_mmol_L"].pop("R1003")
    manifest["fingerprint"] = sha256_payload({
        key: value for key, value in manifest.items() if key != "fingerprint"
    })
    path.write_text(json.dumps(manifest), encoding="utf-8")
    _refresh_manifest_fingerprints(finite)

    with pytest.raises(ValueError, match="initial inventory declaration mismatch"):
        bundle.build_bundle(
            finite, nonlimiting, tmp_path / "bundle",
            compute_commit=bundle.COMPUTE_COMMIT,
            data_sha256=bundle.DATA_SHA256,
            repo_root=REPO, spec=SPEC, job_ids=JOB_IDS,
        )


def test_matching_but_wrong_cross_mode_context_is_rejected(tmp_path):
    runs = [
        _make_run(tmp_path, "finite_batch"),
        _make_run(tmp_path, "po1f_nonlimiting"),
    ]
    for run in runs:
        for path in sorted(run.glob("chunk_*_manifest.json")):
            manifest = json.loads(path.read_text())
            manifest["simulation_context"]["simulation_context_fingerprint"] = "0" * 64
            manifest["fingerprint"] = sha256_payload({
                key: value for key, value in manifest.items() if key != "fingerprint"
            })
            path.write_text(json.dumps(manifest), encoding="utf-8")
        _refresh_manifest_fingerprints(run)

    with pytest.raises(ValueError, match="frozen PO1f simulation context mismatch"):
        bundle.build_bundle(
            runs[0], runs[1], tmp_path / "bundle",
            compute_commit=bundle.COMPUTE_COMMIT,
            data_sha256=bundle.DATA_SHA256,
            repo_root=REPO, spec=SPEC, job_ids=JOB_IDS,
        )
