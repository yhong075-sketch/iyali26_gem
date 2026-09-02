"""Validate two frozen CoQ9 dFBA runs and build a read-only interpretation bundle."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import hashlib
import html
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any
import zipfile
import xml.etree.ElementTree as ET

import pandas as pd

from .coq9_dilution import runtime_coq9_dilution_record
from .essentiality_simulation_context import sha256_file, sha256_payload
from .quinone_dfba_essentiality import (
    INITIAL_POOLS_MMOL_L,
    INTERPRETATION_FLUX_REACTIONS,
    Q9_RESERVE_POLICY,
    SCHEMA_VERSION,
    TRAJECTORY_SEMANTICS,
)
from .summarize_quinone_dfba import FIVE_COQ_CONTROLS, NINE_COQ_CANDIDATES
from .validate_essential_genes import DEFAULT_CUTOFFS


COMPUTE_COMMIT = "36bb6f0735e4c6458bd53c0ceb01952b116b8be7"
MODEL_SHA256 = "bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee"
DATA_SHA256 = "5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d"
MEDIUM_SHA256 = "ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85"
PROFILE_SHA256 = "35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383"
OVERLAY_EFFECT_SHA256 = "d15acbde9438f5d2391c4da23705a34a3585833062d616517d5af052088606c2"
SIMULATION_CONTEXT_SHA256 = "c243b23e7344e3f1e2b4962be25f0f2a38980990c6fe88ee32d3aa4f7af90e30"
EXPERIMENTAL_SHA256 = "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b"
RUNTIME_VERSIONS = {"python": "3.11.15", "cobra": "0.32.1", "gurobipy": "10.0.3"}
ALPHAS = (1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
POOLS = (0.0, 0.5, 1.0, 2.0, 10.0)
GENE_PANEL = (
    "YALI1C26017g", "YALI1F08349g", "YALI1B20835g", "YALI1F34625g",
    "YALI1C25352g", "YALI1A08781g", "YALI1E18269g", "YALI1B20527g",
    "YALI1F34675g", "YALI1A03322g", "YALI1A14736g", "YALI1A17570g",
    "YALI1B02487g", "YALI1C17100g", "YALI1D11769g", "YALI1E40361g",
    "YALI1F02758g", "YALI1F12013g", "YALI1F32010g", "YALI1M00390g",
    "YALI1A14721g", "YALI1D14260g", "YALI1D30256g", "YALI1E35082g",
    "YALI1A21711g", "YALI1B00908g", "YALI1B19507g", "YALI1B26679g",
    "YALI1C04281g", "YALI1D00766g", "YALI1D06302g", "YALI1D07089g",
    "YALI1D09203g", "YALI1D18037g", "YALI1D24109g", "YALI1D32550g",
    "YALI1E06573g", "YALI1E27218g", "YALI1E37603g", "YALI1F01456g",
    "YALI1F03371g", "YALI1F09003g", "YALI1F22993g", "YALI1F24343g",
    "YALI1F32476g", "YALI1M00056g", "YALI1M00064g", "YALI1M00296g",
    "YALI1M00335g", "YALI1M00338r", "YALI1M00458g", "YALIfMp29",
    "YALI1A00309g", "YALI1B02345g", "YALI1B21070g", "YALI1A14287g",
    "YALI1A21372g", "YALI1B26153g", "YALI1F31153g", "YALI1B01084g",
)
MODE_PREFIXES = {
    "finite_batch": "coq9_bio_s18_finite_a9p5_pfba_tol1e9_dt00625_",
    "po1f_nonlimiting": "coq9_bio_s18_nonlimiting_a9p5_pfba_tol1e9_dt00625_",
}
STATE_TOL = 1e-12
SOLVER_TOL = 1e-9
SOURCE_FREE_GROWTH_TOL = 1e-9

CORE_TRAJECTORY_COLUMNS = (
    "gene_id", "uracil_mode", "step_index", "time_h", "time_end_h",
    "interval_advanced", "biomass_gDW_L", "biomass_end_gDW_L", "growth_h-1",
    "biomass_flux_h-1", "objective_value", "source_free_solver_status",
    "source_free_growth_h-1", "q9_pool_mmol_L", "q9_pool_end_mmol_L",
    "q9_source_flux_mmol_gDW_h", "coq9_dilution_flux_mmol_gDW_h",
    "uracil_mmol_L", "uracil_end_mmol_L", "uracil_uptake_flux_mmol_gDW_h",
    "glucose_mmol_L", "glucose_end_mmol_L", "glucose_uptake_flux_mmol_gDW_h",
    "oxygen_uptake_flux_mmol_gDW_h", "atp_maintenance_flux_mmol_gDW_h",
    "status",
)
REACTION_FLUX_COLUMNS = tuple(
    f"reaction_flux_{reaction_id}_mmol_gDW_h"
    for reaction_id in INTERPRETATION_FLUX_REACTIONS
)
GRID_COLUMNS = ("alpha_mmol_gDW", "alpha_replicate_id", "pool_multiplier")
TRAJECTORY_COLUMNS = CORE_TRAJECTORY_COLUMNS + REACTION_FLUX_COLUMNS + GRID_COLUMNS
TRACE_KEY = ["gene_id", "alpha_mmol_gDW", "pool_multiplier"]
CUTOFF_COLUMNS = tuple(f"essential_at_{cutoff * 100:g}pct" for cutoff in DEFAULT_CUTOFFS)
CALL_COLUMNS = (
    "gene_id", "uracil_mode", "final_biomass_gDW_L", "dynamic_doublings",
    "initial_growth_h-1", "initial_source_free_solver_status",
    "initial_source_free_growth_h-1", "q9_pool_depleted_h",
    "q9_source_total_mmol_L", "final_glucose_mmol_L", "final_uracil_mmol_L",
    "termination_status", "termination_time_h", "alpha_mmol_gDW",
    "alpha_replicate_id", "pool_multiplier", "runtime_gpr_scenario",
    "runtime_gpr_mapping_sha256", "wt_dynamic_doublings", "dynamic_growth_ratio",
    "wt_termination_status", "wt_termination_time_h", "experimental_essential",
    *CUTOFF_COLUMNS,
)

COQ_EVIDENCE_CATEGORY = {
    "YALI1C26017g": "uncharacterized; heterologous catalytic-core support",
    "YALI1F08349g": "curated annotation",
    "YALI1B20835g": "model/GPR assignment only",
    "YALI1F34625g": "model/GPR assignment only",
    "YALI1C25352g": "model/GPR assignment only",
    "YALI1A08781g": "uncharacterized; comparative/structure evidence only",
    "YALI1E18269g": "model/GPR assignment only",
    "YALI1B20527g": "uncharacterized; comparative/domain evidence only",
    "YALI1F34675g": "uncharacterized; comparative/domain evidence only",
}
SAFE_PANEL_ANNOTATIONS = {
    "YALI1F32476g": {
        "name_or_symbol": "NDH2",
        "function": "external alternative NADH:ubiquinone oxidoreductase",
        "evidence_status": "experimentally verified",
        "evidence_note": "direct Yarrowia functional evidence; applies to the external NDH2 role",
    },
    "YALI1A21372g": {
        "name_or_symbol": "LIP2",
        "function": "secreted extracellular triacylglycerol lipase",
        "evidence_status": "experimentally verified",
        "evidence_note": "direct Yarrowia functional evidence",
    },
    "YALIfMp29": {
        "name_or_symbol": "unresolved legacy identifier",
        "function": "protein function unverified",
        "evidence_status": "identifier conflict; function not inferred",
        "evidence_note": "legacy YALIfMp29 conflicts with YALI1M00472r; resolve before curation",
    },
    "YALI1D18037g": {
        "name_or_symbol": "no established gene name",
        "function": "protein function unverified",
        "evidence_status": "identifier/function conflict; function not inferred",
        "evidence_note": "ACPM identity conflict remains unresolved; canonical GPR roles are reported without protein-name inference",
    },
}


@dataclass(frozen=True)
class BundleSpec:
    panel: tuple[str, ...] = GENE_PANEL
    alphas: tuple[float, ...] = ALPHAS
    pools: tuple[float, ...] = POOLS
    chunk_count: int = 15
    hours: float = 24.0
    dt: float = 0.0625
    initial_biomass: float = 0.01
    model_sha256: str = MODEL_SHA256
    data_sha256: str = DATA_SHA256
    enforce_run_names: bool = True

    @property
    def combinations(self) -> int:
        return len(self.alphas) * len(self.pools)

    @property
    def calls_per_mode(self) -> int:
        return len(self.panel) * self.combinations


FORMAL_SPEC = BundleSpec()


@dataclass
class ValidatedRun:
    mode: str
    path: Path
    manifests: list[dict[str, Any]]
    calls: pd.DataFrame
    wt_traces: dict[tuple[float, float], pd.DataFrame]
    checks: dict[str, Any]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read JSON artifact {path}: {error}") from error


def _sha_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _git_blob(repo: Path, commit: str, path: str) -> bytes:
    try:
        return subprocess.check_output(
            ["git", "show", f"{commit}:{path}"], cwd=repo, stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as error:
        raise ValueError(f"Missing compute input {commit}:{path}") from error


def _git_head(repo: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=repo, text=True,
    ).strip()


def _strict_bool(series: pd.Series, label: str) -> pd.Series:
    def convert(value: Any) -> bool:
        if isinstance(value, bool):
            return value
        if isinstance(value, str) and value.strip().lower() in {"true", "false"}:
            return value.strip().lower() == "true"
        raise ValueError(f"{label} contains a non-boolean value: {value!r}")

    return series.map(convert)


def _numeric(series: pd.Series, label: str, *, finite: bool = True) -> pd.Series:
    values = pd.to_numeric(series, errors="raise")
    if finite:
        _require(values.map(math.isfinite).all(), f"{label} contains non-finite values")
    return values


def _max_error(left: pd.Series, right: pd.Series) -> float:
    return float((left.astype(float) - right.astype(float)).abs().max()) if len(left) else 0.0


def _assert_close(left: pd.Series, right: pd.Series, tolerance: float, label: str) -> float:
    error = _max_error(left, right)
    _require(error <= tolerance, f"{label} error {error:g} exceeds {tolerance:g}")
    return error


def _explicit_alpha_sampling(spec: BundleSpec) -> list[dict[str, Any]]:
    return [
        {
            "sampler_id": "explicit_alpha_v1", "base_seed": None,
            "replicate_id": None, "distribution": "explicit",
            "low_mmol_gDW": None, "high_mmol_gDW": None,
            "alpha_mmol_gDW": alpha,
        }
        for alpha in spec.alphas
    ]


def _validate_repository(
    repo: Path, compute_commit: str, data_sha256: str, spec: BundleSpec,
) -> dict[str, str]:
    _require(compute_commit == COMPUTE_COMMIT, "compute commit is not the frozen S18 commit")
    _require(data_sha256 == spec.data_sha256, "declared data/iyali26.xml SHA is not frozen")
    model_blob = _git_blob(repo, compute_commit, "model.xml")
    data_blob = _git_blob(repo, compute_commit, "data/iyali26.xml")
    runner_blob = _git_blob(
        repo, compute_commit, "scripts/gem_annotate/quinone_dfba_essentiality.py",
    )
    dilution_blob = _git_blob(repo, compute_commit, "scripts/gem_annotate/coq9_dilution.py")
    _require(_sha_bytes(model_blob) == spec.model_sha256, "compute-commit model SHA mismatch")
    _require(_sha_bytes(data_blob) == data_sha256, "compute-commit data SHA mismatch")
    _require(sha256_file(repo / "model.xml") == spec.model_sha256, "working model.xml changed")
    _require(sha256_file(repo / "data/iyali26.xml") == data_sha256, "working data/iyali26.xml changed")
    postprocessor = Path(__file__).resolve()
    status = subprocess.check_output(
        ["git", "status", "--porcelain", "--", str(postprocessor.relative_to(repo))],
        cwd=repo, text=True,
    ).strip()
    return {
        "compute_commit": compute_commit,
        "postprocessor_commit": _git_head(repo),
        "model_sha256": spec.model_sha256,
        "data_iyali26_sha256": data_sha256,
        "compute_runner_sha256": _sha_bytes(runner_blob),
        "compute_dilution_tool_sha256": _sha_bytes(dilution_blob),
        "postprocessor_sha256": sha256_file(postprocessor),
        "postprocessor_git_state": "committed_clean" if not status else status,
    }


def _validate_manifests(
    run_dir: Path, mode: str, provenance: dict[str, str], spec: BundleSpec,
) -> list[dict[str, Any]]:
    paths = sorted(run_dir.glob("chunk_*_manifest.json"))
    _require(len(paths) == spec.chunk_count, f"{mode}: expected {spec.chunk_count} manifests")
    manifests = [_read_json(path) for path in paths]
    manifests.sort(key=lambda item: item.get("chunk_index", -1))
    _require(
        [item.get("chunk_index") for item in manifests] == list(range(spec.chunk_count)),
        f"{mode}: chunk indices are incomplete",
    )
    expected_files = {"essentiality_dynamic_calls.tsv", "merge_manifest.json"}
    common_payloads = []
    for index, manifest in enumerate(manifests):
        expected_fingerprint = sha256_payload({
            key: value for key, value in manifest.items() if key != "fingerprint"
        })
        _require(manifest.get("fingerprint") == expected_fingerprint, f"{mode} chunk {index}: bad fingerprint")
        expected = {
            "workflow": "quinone_dfba_essentiality",
            "schema_version": SCHEMA_VERSION,
            "trajectory_semantics": TRAJECTORY_SEMANTICS,
            "chunk_index": index,
            "chunk_count": spec.chunk_count,
            "genes": list(spec.panel[index::spec.chunk_count]),
            "solver": "gurobi",
            "solver_feasibility_tolerance": SOLVER_TOL,
            "hours": spec.hours,
            "dt_h": spec.dt,
            "initial_biomass_gDW_L": spec.initial_biomass,
            "alphas_mmol_gDW": list(spec.alphas),
            "alpha_sampling": _explicit_alpha_sampling(spec),
            "pool_multipliers": list(spec.pools),
            "uracil_mode": mode,
            "optimizer": "pfba",
            "nonoptimal_policy": "terminal_record_and_stop",
            "calibration_status": "sensitivity_only_not_calibrated",
            "q9_reserve_definition": "alpha * initial_biomass * pool_multiplier mmol/L",
            "q9_reserve_policy": Q9_RESERVE_POLICY,
            "runtime_topology": {"enabled": False},
            "runtime_gpr_scenario": {"enabled": False},
        }
        for key, value in expected.items():
            _require(manifest.get(key) == value, f"{mode} chunk {index}: manifest {key} mismatch")
        expected_versions = (
            RUNTIME_VERSIONS if spec.enforce_run_names
            else {"gurobipy": RUNTIME_VERSIONS["gurobipy"]}
        )
        _require(
            all(manifest.get("runtime_versions", {}).get(key) == value for key, value in expected_versions.items()),
            f"{mode} chunk {index}: runtime version mismatch",
        )
        context = manifest.get("simulation_context", {})
        expected_context = {
            "simulation_context_fingerprint_version": "1",
            "simulation_context_fingerprint": SIMULATION_CONTEXT_SHA256,
            "strain_overlay_enabled": True,
            "strain_profile_id": "po1f_sd_leu_accrispr_v1",
            "strain_profile_sha256": PROFILE_SHA256,
            "strain_overlay_effect_fingerprint_version": "1",
            "strain_overlay_effect_sha256": OVERLAY_EFFECT_SHA256,
        }
        _require(
            context == expected_context,
            f"{mode} chunk {index}: frozen PO1f simulation context mismatch",
        )
        input_sha = manifest.get("input_sha256", {})
        _require(
            input_sha == {
                "model": spec.model_sha256,
                "medium": MEDIUM_SHA256,
                "profile": PROFILE_SHA256,
                "experimental": EXPERIMENTAL_SHA256,
            },
            f"{mode} chunk {index}: frozen input SHA ledger mismatch",
        )
        base_uracil = manifest.get("base_r1354_bound_mmol_gDW_h")
        _require(
            isinstance(base_uracil, (int, float)) and math.isfinite(base_uracil)
            and base_uracil > 0,
            f"{mode} chunk {index}: invalid base R1354 bound",
        )
        _require(
            manifest.get("script_sha256") == provenance["compute_runner_sha256"],
            f"{mode} chunk {index}: runner SHA mismatch",
        )
        _require(
            manifest.get("coq9_dilution_tool_sha256") == provenance["compute_dilution_tool_sha256"],
            f"{mode} chunk {index}: dilution-tool SHA mismatch",
        )
        _require(
            manifest.get("coq9_dilution") == [
                runtime_coq9_dilution_record(alpha, pool_source_enabled=True)
                for alpha in spec.alphas
            ],
            f"{mode} chunk {index}: dilution records mismatch",
        )
        expected_initial_pools = dict(INITIAL_POOLS_MMOL_L)
        if mode == "po1f_nonlimiting":
            expected_initial_pools["R1354"] = None
        _require(
            manifest.get("initial_pools_mmol_L") == expected_initial_pools,
            f"{mode} chunk {index}: initial inventory declaration mismatch",
        )
        outputs = manifest.get("output_files")
        _require(isinstance(outputs, dict) and set(outputs) == {"calls", "trajectory"}, f"{mode} chunk {index}: output descriptors mismatch")
        for kind in ("calls", "trajectory"):
            filename = f"chunk_{index:03d}_{kind}.tsv"
            descriptor = outputs[kind]
            _require(descriptor.get("filename") == filename, f"{mode} chunk {index}: {kind} filename mismatch")
            path = run_dir / filename
            _require(path.is_file(), f"{mode} chunk {index}: missing {filename}")
            _require(sha256_file(path) == descriptor.get("sha256"), f"{mode} chunk {index}: {kind} SHA mismatch")
            expected_files.add(filename)
        expected_files.add(f"chunk_{index:03d}_manifest.json")
        common_payloads.append({
            key: value for key, value in manifest.items()
            if key not in {"chunk_index", "genes", "output_files", "created_at_utc", "fingerprint"}
        })
    _require(
        len({sha256_payload(item) for item in common_payloads}) == 1,
        f"{mode}: chunk contexts differ",
    )
    observed_files = {path.name for path in run_dir.iterdir() if path.is_file()}
    _require(observed_files == expected_files, f"{mode}: unexpected or missing raw run files")

    merge = _read_json(run_dir / "merge_manifest.json")
    _require(merge.get("workflow") == "quinone_dfba_essentiality", f"{mode}: merge workflow mismatch")
    _require(merge.get("schema_version") == SCHEMA_VERSION, f"{mode}: merge schema mismatch")
    _require(merge.get("trajectory_semantics") == TRAJECTORY_SEMANTICS, f"{mode}: merge trajectory semantics mismatch")
    _require(merge.get("chunks") == spec.chunk_count, f"{mode}: merge chunk count mismatch")
    _require(merge.get("calls") == spec.calls_per_mode, f"{mode}: merge calls count mismatch")
    _require(
        merge.get("chunk_fingerprints") == sorted(item["fingerprint"] for item in manifests),
        f"{mode}: merge fingerprints mismatch",
    )
    merge_expected = {
        "calibration_status": "sensitivity_only_not_calibrated",
        "alpha_sampling": manifests[0]["alpha_sampling"],
        "coq9_dilution": manifests[0]["coq9_dilution"],
        "coq9_dilution_tool_sha256": manifests[0]["coq9_dilution_tool_sha256"],
        "q9_reserve_policy": manifests[0]["q9_reserve_policy"],
        "runtime_gpr_scenario": {"enabled": False},
    }
    for key, value in merge_expected.items():
        _require(merge.get(key) == value, f"{mode}: merge {key} mismatch")
    return manifests


def _validate_calls(
    run_dir: Path, mode: str, manifests: list[dict[str, Any]], spec: BundleSpec,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    frames = []
    for index, manifest in enumerate(manifests):
        frame = pd.read_csv(run_dir / f"chunk_{index:03d}_calls.tsv", sep="\t")
        _require(tuple(frame.columns) == CALL_COLUMNS, f"{mode} chunk {index}: calls fields/order mismatch")
        _require(len(frame) == len(manifest["genes"]) * spec.combinations, f"{mode} chunk {index}: calls shape mismatch")
        _require(set(frame["gene_id"]) == set(manifest["genes"]), f"{mode} chunk {index}: calls genes mismatch")
        _require(not frame.duplicated(["gene_id", "alpha_mmol_gDW", "pool_multiplier"]).any(), f"{mode} chunk {index}: duplicate calls")
        frames.append(frame)
    calls = pd.concat(frames, ignore_index=True)
    _require(len(calls) == spec.calls_per_mode, f"{mode}: calls row count mismatch")
    _require(set(calls["gene_id"]) == set(spec.panel), f"{mode}: calls panel mismatch")
    _require(set(calls["alpha_mmol_gDW"]) == set(spec.alphas), f"{mode}: calls alpha grid mismatch")
    _require(set(calls["pool_multiplier"]) == set(spec.pools), f"{mode}: calls pool grid mismatch")
    _require(not calls.duplicated(["gene_id", "alpha_mmol_gDW", "pool_multiplier"]).any(), f"{mode}: duplicate merged calls")
    _require(calls["alpha_replicate_id"].isna().all(), f"{mode}: explicit alpha replicate IDs must be blank")
    _require(calls["runtime_gpr_scenario"].isna().all(), f"{mode}: runtime GPR scenario is not blank")
    _require(calls["runtime_gpr_mapping_sha256"].isna().all(), f"{mode}: runtime GPR mapping SHA is not blank")
    _require(set(calls["uracil_mode"]) == {mode}, f"{mode}: calls uracil mode mismatch")

    for column in (
        "final_biomass_gDW_L", "dynamic_doublings", "initial_growth_h-1",
        "initial_source_free_growth_h-1", "q9_source_total_mmol_L",
        "final_glucose_mmol_L", "termination_time_h", "wt_dynamic_doublings",
        "dynamic_growth_ratio", "wt_termination_time_h",
    ):
        calls[column] = _numeric(calls[column], f"{mode} {column}")
    _require((calls["final_biomass_gDW_L"] > 0).all(), f"{mode}: nonpositive final biomass")
    _require((calls[["dynamic_doublings", "q9_source_total_mmol_L", "final_glucose_mmol_L", "dynamic_growth_ratio"]] >= 0).all().all(), f"{mode}: negative endpoint")
    _require((calls["wt_dynamic_doublings"] > 0).all(), f"{mode}: nonpositive WT doublings")
    ratio_error = _assert_close(
        calls["dynamic_growth_ratio"],
        calls["dynamic_doublings"] / calls["wt_dynamic_doublings"],
        STATE_TOL, f"{mode} calls ratio",
    )
    for cutoff, column in zip(DEFAULT_CUTOFFS, CUTOFF_COLUMNS):
        calls[column] = _strict_bool(calls[column], f"{mode} {column}")
        _require((calls[column] == (calls["dynamic_growth_ratio"] < cutoff)).all(), f"{mode}: stored {column} mismatch")
    calls["experimental_essential"] = _strict_bool(calls["experimental_essential"], f"{mode} experimental_essential")
    _require(
        not calls.groupby("gene_id")["experimental_essential"].nunique(dropna=False).gt(1).any(),
        f"{mode}: experimental membership changes by condition",
    )
    wt_columns = ["wt_dynamic_doublings", "wt_termination_status", "wt_termination_time_h"]
    _require(
        not calls.groupby(["alpha_mmol_gDW", "pool_multiplier"])[wt_columns].nunique(dropna=False).gt(1).any().any(),
        f"{mode}: WT call endpoints differ",
    )
    reserve = calls["alpha_mmol_gDW"] * spec.initial_biomass * calls["pool_multiplier"]
    source_excess = float((calls["q9_source_total_mmol_L"] - reserve).max())
    _require(source_excess <= STATE_TOL, f"{mode}: calls Q9 source exceeds reserve")
    _require(
        not (calls.loc[calls["pool_multiplier"] == 0, "q9_source_total_mmol_L"] > STATE_TOL).any(),
        f"{mode}: pool-zero calls use Q9 source",
    )
    if mode == "po1f_nonlimiting":
        _require(calls["final_uracil_mmol_L"].isna().all(), f"{mode}: final uracil must be NaN")
    else:
        calls["final_uracil_mmol_L"] = _numeric(calls["final_uracil_mmol_L"], f"{mode} final uracil")
        _require((calls["final_uracil_mmol_L"] >= 0).all(), f"{mode}: negative final uracil")

    pool_violations = 0
    for _, frame in calls.groupby(["gene_id", "alpha_mmol_gDW"], sort=False):
        frame = frame.sort_values("pool_multiplier")
        if (frame["dynamic_doublings"].diff().dropna() < -1e-8).any() or (frame["dynamic_growth_ratio"].diff().dropna() < -1e-8).any():
            pool_violations += 1
        for column in CUTOFF_COLUMNS:
            values = frame[column].tolist()
            if any((not lower) and upper for lower, upper in zip(values, values[1:])):
                pool_violations += 1
    _require(pool_violations == 0, f"{mode}: pool monotonicity failed")
    theory_error = None
    controls = {item["gene_id"] for item in FIVE_COQ_CONTROLS}
    if controls <= set(spec.panel):
        five = calls[calls["gene_id"].isin(controls)].copy()
        expected = five["pool_multiplier"].map(lambda value: math.log2(1.0 + value))
        theory_error = _max_error(five["dynamic_doublings"], expected)
        _require(theory_error <= 1e-6, f"{mode}: five-gene theoretical control failed")
    merged = pd.read_csv(run_dir / "essentiality_dynamic_calls.tsv", sep="\t")
    try:
        pd.testing.assert_frame_equal(
            merged, pd.concat(frames, ignore_index=True), check_dtype=False, check_exact=True,
        )
    except AssertionError as error:
        raise ValueError(f"{mode}: merged calls differ from manifest-bound chunks") from error
    return calls, {
        "calls_rows": len(calls), "calls_ratio_max_abs_error": ratio_error,
        "q9_source_max_excess_mmol_L": source_excess,
        "pool_monotonicity_violations": pool_violations,
        "five_control_max_abs_theory_error_doublings": theory_error,
    }


def _validate_trace_group(
    frame: pd.DataFrame, *, mode: str, manifest: dict[str, Any], spec: BundleSpec,
) -> dict[str, Any]:
    frame = frame.sort_values("step_index").reset_index(drop=True)
    gene_id = str(frame.loc[0, "gene_id"])
    alpha = float(frame.loc[0, "alpha_mmol_gDW"])
    pool = float(frame.loc[0, "pool_multiplier"])
    label = f"{mode}/{gene_id}/{alpha:g}/{pool:g}"
    _require(frame["step_index"].tolist() == list(range(len(frame))), f"{label}: step indices are not contiguous")
    _require(len(frame) <= round(spec.hours / spec.dt), f"{label}: too many intervals")
    _require(float(frame.loc[0, "time_h"]) == 0.0, f"{label}: first time is not zero")
    _require(abs(float(frame.loc[0, "biomass_gDW_L"]) - spec.initial_biomass) <= STATE_TOL, f"{label}: initial biomass mismatch")
    initial_q9 = alpha * spec.initial_biomass * pool
    _require(abs(float(frame.loc[0, "q9_pool_mmol_L"]) - initial_q9) <= STATE_TOL, f"{label}: initial Q9 pool mismatch")
    initial_pools = manifest["initial_pools_mmol_L"]
    _require(abs(float(frame.loc[0, "glucose_mmol_L"]) - float(initial_pools["R1070"])) <= STATE_TOL, f"{label}: initial glucose mismatch")
    advanced = _strict_bool(frame["interval_advanced"], f"{label} interval_advanced")
    terminal_count = int((~advanced).sum())
    _require(terminal_count <= 1, f"{label}: multiple terminal rows")
    if terminal_count:
        _require(not bool(advanced.iloc[-1]) and advanced.iloc[:-1].all(), f"{label}: terminal row is not uniquely last")
        _require(str(frame.iloc[-1]["status"]) != "optimal", f"{label}: terminal row is optimal")
    else:
        _require(advanced.all() and len(frame) == round(spec.hours / spec.dt), f"{label}: completed trajectory has wrong length")
    a = frame.loc[advanced].copy()
    _require((a["status"] == "optimal").all(), f"{label}: advanced row is non-optimal")
    finite_fluxes = [
        "growth_h-1", "biomass_flux_h-1", "objective_value",
        "source_free_growth_h-1", "q9_source_flux_mmol_gDW_h",
        "coq9_dilution_flux_mmol_gDW_h", "uracil_uptake_flux_mmol_gDW_h",
        "glucose_uptake_flux_mmol_gDW_h", "oxygen_uptake_flux_mmol_gDW_h",
        "atp_maintenance_flux_mmol_gDW_h", *REACTION_FLUX_COLUMNS,
    ]
    for column in finite_fluxes:
        a[column] = _numeric(a[column], f"{label} {column}")
    _require((a["source_free_solver_status"] == "optimal").all(), f"{label}: advanced source-free solve is non-optimal")
    nonnegative = [
        "biomass_gDW_L", "biomass_end_gDW_L", "growth_h-1", "biomass_flux_h-1",
        "objective_value", "source_free_growth_h-1", "q9_pool_mmol_L",
        "q9_pool_end_mmol_L", "q9_source_flux_mmol_gDW_h",
        "coq9_dilution_flux_mmol_gDW_h", "uracil_uptake_flux_mmol_gDW_h",
        "glucose_mmol_L", "glucose_end_mmol_L", "glucose_uptake_flux_mmol_gDW_h",
        "oxygen_uptake_flux_mmol_gDW_h", "atp_maintenance_flux_mmol_gDW_h",
    ]
    _require((a[nonnegative] >= 0).all().all(), f"{label}: negative state or directed flux")
    errors = {
        "time": _assert_close(a["time_end_h"], a["time_h"] + spec.dt, STATE_TOL, f"{label} time Euler"),
        "biomass": _assert_close(
            a["biomass_end_gDW_L"],
            a["biomass_gDW_L"] * (1.0 + a["growth_h-1"] * spec.dt),
            STATE_TOL, f"{label} biomass Euler",
        ),
        "biomass_flux": _assert_close(a["biomass_flux_h-1"], a["growth_h-1"], STATE_TOL, f"{label} biomass flux"),
        "objective": _assert_close(a["objective_value"], a["growth_h-1"], STATE_TOL, f"{label} objective"),
        "dilution": _assert_close(
            a["coq9_dilution_flux_mmol_gDW_h"],
            alpha * a["biomass_flux_h-1"], SOLVER_TOL, f"{label} dilution coupling",
        ),
        "glucose": _assert_close(
            a["glucose_end_mmol_L"],
            (a["glucose_mmol_L"] - a["glucose_uptake_flux_mmol_gDW_h"] * a["biomass_gDW_L"] * spec.dt).clip(lower=0),
            STATE_TOL, f"{label} glucose balance",
        ),
    }
    raw_q9_end = (a["q9_pool_mmol_L"] - a["q9_source_flux_mmol_gDW_h"] * a["biomass_gDW_L"] * spec.dt).clip(lower=0)
    q9_floor = max(1e-12, initial_q9 * 1e-9)
    expected_q9_end = raw_q9_end.where(raw_q9_end > q9_floor, 0.0)
    errors["q9"] = _assert_close(a["q9_pool_end_mmol_L"], expected_q9_end, STATE_TOL, f"{label} Q9 balance")
    source_excess = float((a["q9_source_flux_mmol_gDW_h"] * a["biomass_gDW_L"] * spec.dt - a["q9_pool_mmol_L"]).max())
    _require(source_excess <= STATE_TOL, f"{label}: interval Q9 source exceeds pool")
    _require(
        not (a.loc[a["source_free_growth_h-1"] > SOURCE_FREE_GROWTH_TOL, "q9_source_flux_mmol_gDW_h"] > STATE_TOL).any(),
        f"{label}: Q9 source used despite source-free growth",
    )
    r1287 = "reaction_flux_R1287_mmol_gDW_h"
    maintenance = "reaction_flux_xMAINTENANCE_mmol_gDW_h"
    errors["oxygen"] = _assert_close(a["oxygen_uptake_flux_mmol_gDW_h"], (-a[r1287]).clip(lower=0), STATE_TOL, f"{label} oxygen flux")
    errors["maintenance"] = _assert_close(a["atp_maintenance_flux_mmol_gDW_h"], a[maintenance], STATE_TOL, f"{label} maintenance flux")
    if mode == "finite_batch":
        _require(frame["uracil_mmol_L"].notna().all() and frame["uracil_end_mmol_L"].notna().all(), f"{label}: finite uracil state missing")
        _require(abs(float(frame.loc[0, "uracil_mmol_L"]) - float(initial_pools["R1354"])) <= STATE_TOL, f"{label}: initial uracil mismatch")
        _require((frame[["uracil_mmol_L", "uracil_end_mmol_L"]] >= 0).all().all(), f"{label}: negative uracil")
        errors["uracil"] = _assert_close(
            a["uracil_end_mmol_L"],
            (a["uracil_mmol_L"] - a["uracil_uptake_flux_mmol_gDW_h"] * a["biomass_gDW_L"] * spec.dt).clip(lower=0),
            STATE_TOL, f"{label} uracil balance",
        )
    else:
        _require(frame["uracil_mmol_L"].isna().all() and frame["uracil_end_mmol_L"].isna().all(), f"{label}: nonlimiting uracil states are not NaN")
        errors["uracil"] = 0.0
    for previous, current in zip(frame.iloc[:-1].itertuples(index=False), frame.iloc[1:].itertuples(index=False)):
        _require(abs(previous.time_end_h - current.time_h) <= STATE_TOL, f"{label}: time discontinuity")
        for end_column, start_column in (
            ("biomass_end_gDW_L", "biomass_gDW_L"),
            ("q9_pool_end_mmol_L", "q9_pool_mmol_L"),
            ("glucose_end_mmol_L", "glucose_mmol_L"),
        ):
            _require(abs(getattr(previous, end_column) - getattr(current, start_column)) <= STATE_TOL, f"{label}: {start_column} discontinuity")
        if mode == "finite_batch":
            _require(abs(previous.uracil_end_mmol_L - current.uracil_mmol_L) <= STATE_TOL, f"{label}: uracil discontinuity")
    if terminal_count:
        terminal = frame.iloc[-1]
        _require(float(terminal["time_end_h"]) == float(terminal["time_h"]), f"{label}: terminal interval advanced")
        for start, end in (
            ("biomass_gDW_L", "biomass_end_gDW_L"),
            ("q9_pool_mmol_L", "q9_pool_end_mmol_L"),
            ("glucose_mmol_L", "glucose_end_mmol_L"),
        ):
            _require(float(terminal[start]) == float(terminal[end]), f"{label}: terminal {start} changed")
        terminal_nan = [
            "growth_h-1", "biomass_flux_h-1", "objective_value",
            "q9_source_flux_mmol_gDW_h", "coq9_dilution_flux_mmol_gDW_h",
            "uracil_uptake_flux_mmol_gDW_h", "glucose_uptake_flux_mmol_gDW_h",
            "oxygen_uptake_flux_mmol_gDW_h", "atp_maintenance_flux_mmol_gDW_h",
            *REACTION_FLUX_COLUMNS,
        ]
        _require(terminal[terminal_nan].isna().all(), f"{label}: terminal row contains interval flux")
        depleted = float(terminal["glucose_mmol_L"]) <= SOLVER_TOL
        if mode == "finite_batch":
            _require(float(terminal["uracil_mmol_L"]) == float(terminal["uracil_end_mmol_L"]), f"{label}: terminal uracil changed")
            depleted = depleted or float(terminal["uracil_mmol_L"]) <= SOLVER_TOL
        _require(depleted, f"{label}: terminal row has no observable depleted batch inventory")
    return {
        "frame": frame, "errors": errors, "source_excess": source_excess,
        "termination_status": str(frame.iloc[-1]["status"]) if terminal_count else "completed",
        "termination_time_h": float(frame.iloc[-1]["time_h"]) if terminal_count else spec.hours,
        "final_biomass_gDW_L": float(frame.iloc[-1]["biomass_end_gDW_L"]),
        "final_glucose_mmol_L": float(frame.iloc[-1]["glucose_end_mmol_L"]),
        "final_uracil_mmol_L": float(frame.iloc[-1]["uracil_end_mmol_L"]) if mode == "finite_batch" else math.nan,
        "final_q9_pool_mmol_L": float(frame.iloc[-1]["q9_pool_end_mmol_L"]),
    }


def _validate_trajectories(
    run_dir: Path, mode: str, manifests: list[dict[str, Any]], calls: pd.DataFrame,
    spec: BundleSpec,
) -> tuple[dict[tuple[float, float], pd.DataFrame], dict[str, Any]]:
    call_index = calls.set_index(["gene_id", "alpha_mmol_gDW", "pool_multiplier"])
    wt_calls = calls.groupby(["alpha_mmol_gDW", "pool_multiplier"], sort=False).first()
    wt_traces: dict[tuple[float, float], pd.DataFrame] = {}
    max_errors = {key: 0.0 for key in ("time", "biomass", "biomass_flux", "objective", "dilution", "glucose", "q9", "oxygen", "maintenance", "uracil")}
    max_source_excess = -math.inf
    groups_seen = 0
    terminal_rows = 0
    for index, manifest in enumerate(manifests):
        path = run_dir / f"chunk_{index:03d}_trajectory.tsv"
        frame = pd.read_csv(path, sep="\t")
        _require(tuple(frame.columns) == TRAJECTORY_COLUMNS, f"{mode} chunk {index}: trajectory fields/order mismatch")
        _require(frame["alpha_replicate_id"].isna().all(), f"{mode} chunk {index}: replicate IDs must be blank")
        _require(set(frame["gene_id"]) == {"WT", *manifest["genes"]}, f"{mode} chunk {index}: trajectory genes mismatch")
        _require(set(frame["uracil_mode"]) == {mode}, f"{mode} chunk {index}: trajectory mode mismatch")
        observed_keys = set(zip(frame["gene_id"], frame["alpha_mmol_gDW"], frame["pool_multiplier"]))
        expected_keys = {
            (gene_id, alpha, pool)
            for gene_id in ("WT", *manifest["genes"])
            for alpha in spec.alphas for pool in spec.pools
        }
        _require(observed_keys == expected_keys, f"{mode} chunk {index}: trajectory groups mismatch")
        for key, group in frame.groupby(TRACE_KEY, sort=False, dropna=False):
            result = _validate_trace_group(group, mode=mode, manifest=manifest, spec=spec)
            groups_seen += 1
            terminal_rows += int(result["termination_status"] != "completed")
            for name, value in result["errors"].items():
                max_errors[name] = max(max_errors[name], value)
            max_source_excess = max(max_source_excess, result["source_excess"])
            gene_id, alpha, pool = key
            if gene_id == "WT":
                wt_key = (float(alpha), float(pool))
                canonical = wt_traces.get(wt_key)
                trace = result["frame"].loc[:, TRAJECTORY_COLUMNS]
                advanced_trace = trace[_strict_bool(trace["interval_advanced"], f"{mode}/WT interval_advanced")]
                _require(
                    not (advanced_trace["q9_source_flux_mmol_gDW_h"] > STATE_TOL).any(),
                    f"{mode}: WT uses Q9 reserve at {wt_key}",
                )
                if canonical is None:
                    wt_traces[wt_key] = trace
                else:
                    try:
                        pd.testing.assert_frame_equal(canonical, trace, check_dtype=False, check_exact=True)
                    except AssertionError as error:
                        raise ValueError(f"{mode}: WT trajectory differs across chunks at {wt_key}") from error
                wt_call = wt_calls.loc[wt_key]
                _require(abs(math.log2(result["final_biomass_gDW_L"] / spec.initial_biomass) - float(wt_call["wt_dynamic_doublings"])) <= STATE_TOL, f"{mode}: WT doublings endpoint mismatch")
                _require(result["termination_status"] == wt_call["wt_termination_status"], f"{mode}: WT termination status mismatch")
                _require(abs(result["termination_time_h"] - float(wt_call["wt_termination_time_h"])) <= STATE_TOL, f"{mode}: WT termination time mismatch")
                continue
            call = call_index.loc[(gene_id, float(alpha), float(pool))]
            final_doublings = math.log2(result["final_biomass_gDW_L"] / spec.initial_biomass)
            _require(abs(final_doublings - float(call["dynamic_doublings"])) <= STATE_TOL, f"{mode}/{key}: dynamic doublings endpoint mismatch")
            _require(abs(result["final_biomass_gDW_L"] - float(call["final_biomass_gDW_L"])) <= STATE_TOL, f"{mode}/{key}: final biomass mismatch")
            _require(abs(result["final_glucose_mmol_L"] - float(call["final_glucose_mmol_L"])) <= STATE_TOL, f"{mode}/{key}: final glucose mismatch")
            if mode == "finite_batch":
                _require(abs(result["final_uracil_mmol_L"] - float(call["final_uracil_mmol_L"])) <= STATE_TOL, f"{mode}/{key}: final uracil mismatch")
            _require(result["termination_status"] == call["termination_status"], f"{mode}/{key}: termination status mismatch")
            _require(abs(result["termination_time_h"] - float(call["termination_time_h"])) <= STATE_TOL, f"{mode}/{key}: termination time mismatch")
            initial_q9 = float(alpha) * spec.initial_biomass * float(pool)
            used_q9 = initial_q9 - result["final_q9_pool_mmol_L"]
            _require(abs(used_q9 - float(call["q9_source_total_mmol_L"])) <= STATE_TOL, f"{mode}/{key}: Q9 endpoint mismatch")
            first = result["frame"].iloc[0]
            _require(abs(float(first["growth_h-1"]) - float(call["initial_growth_h-1"])) <= STATE_TOL, f"{mode}/{key}: initial growth mismatch")
            _require(str(first["source_free_solver_status"]) == str(call["initial_source_free_solver_status"]), f"{mode}/{key}: initial source-free status mismatch")
            _require(abs(float(first["source_free_growth_h-1"]) - float(call["initial_source_free_growth_h-1"])) <= STATE_TOL, f"{mode}/{key}: initial source-free growth mismatch")
            depleted = 0.0 if initial_q9 == 0 else math.nan
            if initial_q9 > 0:
                depleted_rows = result["frame"][(result["frame"]["q9_pool_mmol_L"] > 0) & (result["frame"]["q9_pool_end_mmol_L"] == 0)]
                if not depleted_rows.empty:
                    depleted = float(depleted_rows.iloc[0]["time_end_h"])
            observed_depleted = call["q9_pool_depleted_h"]
            _require(
                (pd.isna(depleted) and pd.isna(observed_depleted))
                or (not pd.isna(depleted) and abs(float(observed_depleted) - depleted) <= STATE_TOL),
                f"{mode}/{key}: Q9 depletion time mismatch",
            )
    expected_groups = spec.chunk_count * (1 + len(spec.panel) // spec.chunk_count) * spec.combinations
    _require(groups_seen == expected_groups, f"{mode}: trajectory group count mismatch")
    _require(len(wt_traces) == spec.combinations, f"{mode}: WT trajectory grid incomplete")
    return wt_traces, {
        "raw_trajectory_groups": groups_seen,
        "unique_trajectory_groups_after_wt_dedup": (len(spec.panel) + 1) * spec.combinations,
        "wt_duplicate_copies_verified_per_condition": spec.chunk_count,
        "terminal_rows": terminal_rows,
        "wt_q9_source_zero": True,
        "max_abs_errors": max_errors,
        "max_interval_q9_source_excess_mmol_L": max_source_excess,
        "trajectory_pass": True,
    }


def _validate_run(
    path: Path, mode: str, provenance: dict[str, str], spec: BundleSpec,
) -> ValidatedRun:
    _require(path.is_dir(), f"{mode}: run directory does not exist")
    manifests = _validate_manifests(path, mode, provenance, spec)
    calls, call_checks = _validate_calls(path, mode, manifests, spec)
    wt_traces, trajectory_checks = _validate_trajectories(path, mode, manifests, calls, spec)
    return ValidatedRun(
        mode=mode, path=path, manifests=manifests, calls=calls,
        wt_traces=wt_traces, checks={**call_checks, **trajectory_checks},
    )


def _validate_cross_mode(runs: list[ValidatedRun], spec: BundleSpec) -> str:
    _require({run.mode for run in runs} == set(MODE_PREFIXES), "exactly the two uracil modes are required")
    contexts = []
    stamps = []
    context_fields = (
        "schema_version", "trajectory_semantics", "solver", "runtime_versions",
        "solver_feasibility_tolerance", "script_sha256", "coq9_dilution_tool_sha256",
        "hours", "dt_h", "initial_biomass_gDW_L", "alphas_mmol_gDW",
        "alpha_sampling", "pool_multipliers", "q9_reserve_definition",
        "q9_reserve_policy", "optimizer", "nonoptimal_policy", "calibration_status",
        "coq9_dilution", "runtime_topology", "runtime_gpr_scenario",
        "simulation_context", "input_sha256", "base_r1354_bound_mmol_gDW_h",
    )
    for run in runs:
        manifest = run.manifests[0]
        contexts.append({field: manifest.get(field) for field in context_fields})
        prefix = MODE_PREFIXES[run.mode]
        if spec.enforce_run_names:
            _require(str(manifest["run_id"]).startswith(prefix), f"{run.mode}: run ID prefix mismatch")
            stamps.append(str(manifest["run_id"])[len(prefix):])
        _require(run.path.name == manifest["run_id"], f"{run.mode}: run directory/name mismatch")
    _require(len({sha256_payload(value) for value in contexts}) == 1, "cross-mode frozen contexts differ")
    finite = next(run for run in runs if run.mode == "finite_batch")
    nonlimiting = next(run for run in runs if run.mode == "po1f_nonlimiting")
    finite_pools = dict(finite.manifests[0]["initial_pools_mmol_L"])
    nonlimiting_pools = dict(nonlimiting.manifests[0]["initial_pools_mmol_L"])
    finite_uracil = finite_pools.pop("R1354", None)
    nonlimiting_uracil = nonlimiting_pools.pop("R1354", "missing")
    _require(finite_uracil == 0.178428 and nonlimiting_uracil is None, "cross-mode uracil inventory contract mismatch")
    _require(finite_pools == nonlimiting_pools, "cross-mode non-uracil initial pools differ")
    finite_membership = finite.calls.groupby("gene_id")["experimental_essential"].first()
    nonlimiting_membership = nonlimiting.calls.groupby("gene_id")["experimental_essential"].first()
    _require(finite_membership.equals(nonlimiting_membership), "experimental membership differs across modes")
    if spec.enforce_run_names:
        _require(len(set(stamps)) == 1 and bool(stamps[0]), "cross-mode run IDs do not share one UTC stamp")
    return sha256_payload(contexts[0])


def _load_coq_evidence(repo: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = repo / "docs" / "research" / "quinone_gpr_synthome_2026-08-17"
    manifest = _read_json(root / "research_manifest.json")
    files = manifest.get("files", manifest)
    for filename in ("gene_evidence_matrix.tsv", "fitness_threshold_matrix.tsv"):
        _require(
            files.get(filename) == sha256_file(root / filename),
            f"committed evidence SHA mismatch: {filename}",
        )
    evidence = pd.read_csv(root / "gene_evidence_matrix.tsv", sep="\t")
    fitness = pd.read_csv(root / "fitness_threshold_matrix.tsv", sep="\t")
    _require(set(evidence["gene_id"]) == set(NINE_COQ_CANDIDATES), "CoQ evidence matrix does not cover nine candidates")
    _require(set(fitness["gene_id"]) == set(NINE_COQ_CANDIDATES), "fitness matrix does not cover nine candidates")
    _require(not evidence["gene_id"].duplicated().any() and not fitness["gene_id"].duplicated().any(), "CoQ evidence matrices contain duplicate genes")
    return evidence, fitness


def _selected_panel(spec: BundleSpec, evidence: pd.DataFrame) -> pd.DataFrame:
    metadata = evidence.set_index("gene_id")
    rows = []
    for order, gene_id in enumerate(spec.panel, 1):
        row = {"panel_order": order, "gene_id": gene_id, "is_coq_candidate": gene_id in NINE_COQ_CANDIDATES}
        if gene_id in metadata.index:
            record = metadata.loc[gene_id]
            row.update({
                "name_or_symbol": record["name_or_symbol"],
                "function": record["function"],
                "evidence_status": record["evidence_status"],
                "evidence_class": COQ_EVIDENCE_CATEGORY[gene_id],
                "evidence_note": record["evidence_status"],
            })
        elif gene_id in SAFE_PANEL_ANNOTATIONS:
            row.update({
                **SAFE_PANEL_ANNOTATIONS[gene_id],
                "evidence_class": SAFE_PANEL_ANNOTATIONS[gene_id]["evidence_status"],
            })
        else:
            row.update({
                "name_or_symbol": "no established gene name",
                "function": "protein function unverified",
                "evidence_status": "canonical model/GPR role only; protein identity not inferred",
                "evidence_class": "model/GPR assignment only",
                "evidence_note": "Only canonical reaction associations are reported.",
            })
        rows.append(row)
    return pd.DataFrame(rows)


def _render_gpr(node: ET.Element) -> str:
    tag = node.tag.rsplit("}", 1)[-1]
    if tag == "geneProductRef":
        value = next(value for key, value in node.attrib.items() if key.endswith("geneProduct"))
        return value.removeprefix("G_")
    children = [_render_gpr(child) for child in node if child.tag.rsplit("}", 1)[-1] in {"and", "or", "geneProductRef"}]
    if not children:
        return ""
    if tag in {"and", "or"}:
        return "(" + f" {tag} ".join(children) + ")"
    return children[0]


def _reaction_gene_map(repo: Path, spec: BundleSpec, panel: pd.DataFrame) -> pd.DataFrame:
    root = ET.parse(repo / "model.xml").getroot()
    metadata = panel.set_index("gene_id")
    associations: dict[str, list[tuple[str, str, str]]] = {gene_id: [] for gene_id in spec.panel}
    for reaction in root.iter():
        if reaction.tag.rsplit("}", 1)[-1] != "reaction":
            continue
        reaction_id = reaction.attrib.get("id", "").removeprefix("R_")
        reaction_name = reaction.attrib.get("name", "")
        association = next(
            (child for child in reaction if child.tag.rsplit("}", 1)[-1] == "geneProductAssociation"),
            None,
        )
        if association is None:
            continue
        rule = _render_gpr(association)
        gene_refs = {
            value.removeprefix("G_") for element in association.iter()
            for key, value in element.attrib.items() if key.endswith("geneProduct")
        }
        for gene_id in spec.panel:
            if gene_id in gene_refs:
                associations[gene_id].append((reaction_id, reaction_name, rule))
    rows = []
    for gene_id in spec.panel:
        reactions = sorted(associations[gene_id])
        if not reactions:
            reactions = [("", "", "")]
        for reaction_id, reaction_name, rule in reactions:
            mapping_evidence = (
                "canonical model/GPR; not re-curated by this runtime-only bundle"
                if reaction_id else "no canonical reaction association"
            )
            if gene_id == "YALI1F32476g" and reaction_id != "R570":
                mapping_evidence = "canonical model/GPR assignment only; NDH2 experimental function does not verify this reaction membership"
            if gene_id == "YALI1F32476g" and reaction_id == "R570":
                mapping_evidence = "experimentally verified Yarrowia external NDH2 reaction association"
            if gene_id == "YALI1A21372g" and reaction_id != "R2274":
                mapping_evidence = "canonical model/GPR assignment only; LIP2 evidence does not independently verify this reaction membership"
            if gene_id == "YALI1A21372g" and reaction_id == "R2274":
                mapping_evidence = "experimentally verified Yarrowia LIP2 reaction association"
            if gene_id in {"YALIfMp29", "YALI1D18037g"}:
                mapping_evidence += "; identity conflict unresolved"
            rows.append({
                "gene_id": gene_id,
                "name_or_symbol": metadata.loc[gene_id, "name_or_symbol"],
                "function": metadata.loc[gene_id, "function"],
                "evidence_status": metadata.loc[gene_id, "evidence_status"],
                "reaction_id": reaction_id,
                "reaction_name": reaction_name,
                "gene_reaction_rule": rule,
                "mapping_evidence": mapping_evidence,
            })
    return pd.DataFrame(rows)


def _calls_tables(conditions: pd.DataFrame, spec: BundleSpec) -> dict[str, pd.DataFrame]:
    grid_rows = []
    for (mode, alpha, pool), frame in conditions.groupby(
        ["uracil_mode", "alpha_mmol_gDW", "pool_multiplier"], sort=True,
    ):
        source = frame["q9_source_total_mmol_L"] > STATE_TOL
        row = {
            "uracil_mode": mode, "alpha_mmol_gDW": alpha,
            "pool_multiplier": pool, "genes_n": frame["gene_id"].nunique(),
            "ratio_min": frame["dynamic_growth_ratio"].min(),
            "ratio_median": frame["dynamic_growth_ratio"].median(),
            "ratio_max": frame["dynamic_growth_ratio"].max(),
            "q9_source_user_genes_n": frame.loc[source, "gene_id"].nunique(),
        }
        row.update({f"{column}_n": int(frame[column].sum()) for column in CUTOFF_COLUMNS})
        grid_rows.append(row)
    grid = pd.DataFrame(grid_rows)

    alpha_rows = []
    for (mode, pool), frame in conditions.groupby(["uracil_mode", "pool_multiplier"], sort=True):
        per_gene = frame.groupby("gene_id")["dynamic_growth_ratio"].agg(lambda values: values.max() - values.min())
        row = {
            "uracil_mode": mode, "pool_multiplier": pool,
            "genes_n": len(per_gene), "ratio_range_median": per_gene.median(),
            "ratio_range_max": per_gene.max(),
        }
        for column in CUTOFF_COLUMNS:
            row[f"{column}_changed_genes_n"] = int(
                frame.groupby("gene_id")[column].nunique().gt(1).sum()
            )
        alpha_rows.append(row)

    pool_rows = []
    for (mode, gene_id, alpha), frame in conditions.groupby(
        ["uracil_mode", "gene_id", "alpha_mmol_gDW"], sort=True,
    ):
        frame = frame.sort_values("pool_multiplier")
        pool_rows.append({
            "uracil_mode": mode, "gene_id": gene_id, "alpha_mmol_gDW": alpha,
            "minimum_doublings_delta": frame["dynamic_doublings"].diff().dropna().min(),
            "minimum_ratio_delta": frame["dynamic_growth_ratio"].diff().dropna().min(),
            "pass": True,
        })

    source = conditions.copy()
    reserve = source["alpha_mmol_gDW"] * spec.initial_biomass * source["pool_multiplier"]
    source["source_fraction_of_initial_reserve"] = source["q9_source_total_mmol_L"].where(
        reserve > 0, 0.0,
    ) / reserve.where(reserve > 0, 1.0)
    source = source[
        ["uracil_mode", "gene_id", "alpha_mmol_gDW", "pool_multiplier",
         "q9_source_total_mmol_L", "source_fraction_of_initial_reserve"]
    ]

    controls = {item["gene_id"] for item in FIVE_COQ_CONTROLS}
    five = conditions[conditions["gene_id"].isin(controls)].copy()
    five["theoretical_doublings"] = five["pool_multiplier"].map(lambda value: math.log2(1.0 + value))
    five["theory_error_doublings"] = five["dynamic_doublings"] - five["theoretical_doublings"]
    five = five[
        ["uracil_mode", "gene_id", "alpha_mmol_gDW", "pool_multiplier",
         "dynamic_doublings", "theoretical_doublings", "theory_error_doublings",
         "dynamic_growth_ratio", *CUTOFF_COLUMNS]
    ]
    return {
        "grid_summary": grid,
        "alpha_sensitivity": pd.DataFrame(alpha_rows),
        "pool_monotonicity": pd.DataFrame(pool_rows),
        "q9_source_use": source,
        "five_gene_theoretical_controls": five,
    }


def _validated_job_ids(
    job_ids: dict[str, dict[str, str]] | None, spec: BundleSpec,
) -> dict[str, dict[str, str]]:
    if job_ids is None:
        _require(not spec.enforce_run_names, "formal bundle requires four SLURM array/merge job IDs")
        return {}
    _require(set(job_ids) == set(MODE_PREFIXES), "SLURM job ledger must cover both uracil modes")
    validated: dict[str, dict[str, str]] = {}
    for mode in MODE_PREFIXES:
        jobs = job_ids[mode]
        _require(set(jobs) == {"array", "merge"}, f"{mode}: SLURM job ledger must contain array and merge")
        validated[mode] = {key: str(value) for key, value in jobs.items()}
        _require(
            all(value.isdigit() for value in validated[mode].values()),
            f"{mode}: SLURM job IDs must be numeric",
        )
    return validated


def _write_trajectory_outputs(runs: list[ValidatedRun], output: Path) -> tuple[int, int]:
    core_path = output / "stepwise_trajectories.tsv.gz"
    flux_path = output / "flux_panel.tsv.gz"
    core_header = True
    flux_header = True
    rows = 0
    with gzip.open(core_path, "wt", encoding="utf-8", newline="") as core_handle, gzip.open(flux_path, "wt", encoding="utf-8", newline="") as flux_handle:
        for run in runs:
            for index in range(len(run.manifests)):
                path = run.path / f"chunk_{index:03d}_trajectory.tsv"
                _require(
                    sha256_file(path) == run.manifests[index]["output_files"]["trajectory"]["sha256"],
                    f"{run.mode} chunk {index}: trajectory changed after validation",
                )
                frame = pd.read_csv(path, sep="\t")
                if index:
                    frame = frame[frame["gene_id"] != "WT"].copy()
                frame.insert(0, "source_chunk", index)
                frame.insert(0, "mode", run.mode)
                keys = ["mode", "source_chunk", "gene_id", "alpha_mmol_gDW", "pool_multiplier", "step_index", "time_h", "time_end_h"]
                core = frame[[
                    "mode", "source_chunk", *CORE_TRAJECTORY_COLUMNS,
                    *REACTION_FLUX_COLUMNS, *GRID_COLUMNS,
                ]]
                flux = frame[[*keys, *REACTION_FLUX_COLUMNS]]
                core.to_csv(core_handle, sep="\t", index=False, header=core_header, na_rep="NaN")
                flux.to_csv(flux_handle, sep="\t", index=False, header=flux_header, na_rep="NaN")
                core_header = flux_header = False
                rows += len(frame)
    return rows, rows


def _write_plot(conditions: pd.DataFrame, plots: Path) -> None:
    plots.mkdir()
    counts = conditions.groupby("uracil_mode")[list(CUTOFF_COLUMNS)].sum().T
    modes = list(counts.columns)
    maximum = max(1, int(counts.max().max()))
    colors = ("#4C78A8", "#F58518")
    bars = []
    for cutoff_index, cutoff in enumerate(CUTOFF_COLUMNS):
        for mode_index, mode in enumerate(modes):
            value = int(counts.loc[cutoff, mode])
            x = 85 + cutoff_index * 145 + mode_index * 50
            height = 250 * value / maximum
            bars.append(
                f'<rect x="{x}" y="{330-height:.2f}" width="38" height="{height:.2f}" fill="{colors[mode_index]}"/>'
                f'<text x="{x+19}" y="{320-height:.2f}" text-anchor="middle" font-size="11">{value}</text>'
            )
    labels = "".join(
        f'<text x="{110 + index*145}" y="355" text-anchor="middle" font-size="12">{html.escape(column.removeprefix("essential_at_").removesuffix("pct"))}%</text>'
        for index, column in enumerate(CUTOFF_COLUMNS)
    )
    legend = "".join(
        f'<rect x="{430}" y="{35 + i*20}" width="12" height="12" fill="{colors[i]}"/>'
        f'<text x="448" y="{46 + i*20}" font-size="11">{html.escape(mode)}</text>'
        for i, mode in enumerate(modes)
    )
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="700" height="390" viewBox="0 0 700 390">'
        '<rect width="100%" height="100%" fill="white"/>'
        '<text x="350" y="22" text-anchor="middle" font-size="16">CoQ9 H-Q9-1 cutoff calls — full frozen grid</text>'
        '<line x1="65" y1="330" x2="660" y2="330" stroke="black"/>'
        + "".join(bars) + labels + legend
        + '<text x="350" y="382" text-anchor="middle" font-size="11">No parameter was selected by recall.</text></svg>'
    )
    (plots / "cutoff_call_counts.svg").write_text(svg, encoding="utf-8")


def _write_report(
    path: Path, runs: list[ValidatedRun], conditions: pd.DataFrame,
    evidence: pd.DataFrame, context_sha: str, spec: BundleSpec,
) -> None:
    meta = evidence.set_index("gene_id")
    coq = conditions[conditions["gene_id"].isin(NINE_COQ_CANDIDATES)]
    lines = [
        "# CoQ9 biological interpretation bundle", "",
        "## Scope", "",
        "This is a **runtime_only** H-Q9-1 sensitivity analysis with status "
        "**sensitivity_only_not_calibrated**. It does not change `model.xml`, "
        "`data/iyali26.xml`, GPRs, bounds, stoichiometry, curated data, or the FN dossier. "
        "No alpha or pool multiplier was selected by experimental recall.", "",
        "## Validated computation", "",
        f"- Frozen compute commit: `{COMPUTE_COMMIT}`",
        f"- Cross-mode context SHA256: `{context_sha}`",
        f"- Grid: {len(spec.alphas)} explicit alpha values × {len(spec.pools)} pool multipliers × 2 uracil modes",
        f"- Gene-condition rows: {len(conditions)}",
        "- Numerical result: every manifest, input/output hash, interval balance, trajectory continuity, "
        "dilution coupling, source-free policy, terminal record, endpoint, and repeated-WT comparison passed for all recorded states.", "",
        "## Full-grid interpretation checks", "",
        "The bundle reports every frozen condition. It includes cutoff counts, alpha sensitivity, pool monotonicity, "
        "Q9-source use, and the five-gene theoretical reserve controls as separate TSV tables. These are sensitivity summaries, not calibration targets.", "",
        "## Nine CoQ candidates", "",
        "A binary model association describes the imposed model dependency; it does not by itself verify protein function.", "",
        "| Gene — name/symbol — function (evidence status) | Mode | KO/WT ratio range | Calls at 1/5/10/15% |",
        "|---|---:|---:|---:|",
    ]
    for gene_id in NINE_COQ_CANDIDATES:
        record = meta.loc[gene_id]
        identity = f"{gene_id} — {record['name_or_symbol']} — {record['function']} ({record['evidence_status']})"
        for mode in MODE_PREFIXES:
            frame = coq[(coq["gene_id"] == gene_id) & (coq["uracil_mode"] == mode)]
            calls = "/".join(str(int(frame[column].sum())) for column in CUTOFF_COLUMNS)
            lines.append(
                f"| {identity} | {mode} | {frame['dynamic_growth_ratio'].min():.6g}–{frame['dynamic_growth_ratio'].max():.6g} | {calls} of {spec.combinations} |"
            )
    lines.extend([
        "", "## Interpretation limits", "",
        "The grid asks what follows under the imposed growth-coupled net CoQ9-demand and finite-reserve assumptions. "
        "It does not estimate a native CoQ9 pool, turnover rate, degradation half-life, synthome stoichiometry, or respiratory capacity.", "",
        "Cas9/Cas12a agreement or disagreement is shown for comparison only. It was not used to tune alpha, pool size, "
        "reaction mappings, or cutoffs. Results require separate human review before any curation or publication claim.", "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_sha256s(root: Path) -> None:
    paths = sorted(path for path in root.rglob("*") if path.is_file() and path.name != "SHA256SUMS")
    (root / "SHA256SUMS").write_text(
        "".join(f"{sha256_file(path)}  {path.relative_to(root)}\n" for path in paths),
        encoding="utf-8",
    )


def build_bundle(
    finite_dir: Path, nonlimiting_dir: Path, output_dir: Path, *,
    compute_commit: str, data_sha256: str, repo_root: Path | None = None,
    spec: BundleSpec = FORMAL_SPEC,
    job_ids: dict[str, dict[str, str]] | None = None,
) -> tuple[Path, Path]:
    """Fail fast on any mismatch, then atomically create the read-only bundle."""
    repo = (repo_root or Path(__file__).resolve().parents[2]).resolve()
    finite_dir = finite_dir.resolve()
    nonlimiting_dir = nonlimiting_dir.resolve()
    output_dir = output_dir.resolve()
    zip_path = Path(f"{output_dir}.zip")
    zip_sha_path = Path(f"{zip_path}.sha256")
    _require(len(spec.panel) == len(set(spec.panel)), "selected gene panel contains duplicates")
    _require(spec.chunk_count > 0 and len(spec.panel) % spec.chunk_count == 0, "panel must divide evenly across chunks")
    _require(set(NINE_COQ_CANDIDATES) <= set(spec.panel) or spec is not FORMAL_SPEC, "formal panel omits a CoQ candidate")
    if output_dir.exists() or zip_path.exists() or zip_sha_path.exists():
        raise FileExistsError("output directory/archive already exists")
    validated_job_ids = _validated_job_ids(job_ids, spec)
    provenance = _validate_repository(repo, compute_commit, data_sha256, spec)
    runs = [
        _validate_run(finite_dir, "finite_batch", provenance, spec),
        _validate_run(nonlimiting_dir, "po1f_nonlimiting", provenance, spec),
    ]
    context_sha = _validate_cross_mode(runs, spec)
    evidence, fitness = _load_coq_evidence(repo)
    conditions = pd.concat(
        [run.calls.assign(source_run_id=run.manifests[0]["run_id"]) for run in runs],
        ignore_index=True,
    )
    _require(len(conditions) == 2 * spec.calls_per_mode, "combined condition table shape mismatch")
    panel = _selected_panel(spec, evidence)
    reaction_map = _reaction_gene_map(repo, spec, panel)
    roles = reaction_map.groupby("gene_id")["reaction_id"].agg(
        lambda values: ";".join(value for value in values if value) or "none",
    )
    named_roles = {
        gene_id: "; ".join(
            f"{row.reaction_id} — {row.reaction_name or 'unnamed reaction'}"
            for row in frame.itertuples() if row.reaction_id
        ) or "no canonical reaction association"
        for gene_id, frame in reaction_map.groupby("gene_id")
    }
    panel["canonical_model_reactions"] = panel["gene_id"].map(roles)
    panel["canonical_model_role"] = panel["canonical_model_reactions"].map(
        lambda value: "no canonical reaction association" if value == "none" else f"canonical GPR association: {value}",
    )
    unverified = panel["function"].eq("protein function unverified")
    panel.loc[unverified, "function"] = panel.loc[unverified, "gene_id"].map(
        lambda gene_id: f"protein function unverified; canonical model/GPR role only: {named_roles[gene_id]}",
    )
    panel_metadata = panel.set_index("gene_id")
    for column in ("name_or_symbol", "function", "evidence_status"):
        reaction_map[column] = reaction_map["gene_id"].map(panel_metadata[column])
    tables = _calls_tables(conditions, spec)
    mapping_payload = {
        "mapping_id": "canonical_gpr_h_q9_1_no_optional_topology_v1",
        "scope": "runtime_only",
        "canonical_model_sha256": spec.model_sha256,
        "q9_reserve_definition": runs[0].manifests[0]["q9_reserve_definition"],
        "q9_reserve_policy": runs[0].manifests[0]["q9_reserve_policy"],
        "coq9_dilution_tool_sha256": provenance["compute_dilution_tool_sha256"],
        "runtime_topology": runs[0].manifests[0]["runtime_topology"],
        "runtime_gpr_scenario": runs[0].manifests[0]["runtime_gpr_scenario"],
    }
    runtime_mapping = {
        **mapping_payload,
        "mapping_sha256": sha256_payload(mapping_payload),
    }
    experimental = conditions.rename(columns={"experimental_essential": "positive_only_consensus_member"})
    experimental = experimental.merge(fitness, on="gene_id", how="left", validate="many_to_one")
    experimental["experimental_evidence_scope"] = experimental["cas9_paper_call"].notna().map({
        True: "Cas9/Cas12a evidence matrix plus positive-only consensus membership",
        False: "positive-only consensus membership only",
    })

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    zip_fd, zip_name = tempfile.mkstemp(prefix=f".{zip_path.name}.", dir=output_dir.parent)
    os.close(zip_fd)
    temporary_zip = Path(zip_name)
    temporary_zip.unlink()
    temporary_zip_sha = Path(f"{temporary_zip}.sha256")
    try:
        panel.to_csv(temporary / "selected_gene_panel.tsv", sep="\t", index=False)
        reaction_map.to_csv(temporary / "reaction_gene_map.tsv", sep="\t", index=False)
        conditions.to_csv(temporary / "gene_condition_summary.tsv", sep="\t", index=False)
        experimental.to_csv(temporary / "experimental_comparison.tsv", sep="\t", index=False)
        for name, table in tables.items():
            table.to_csv(temporary / f"{name}.tsv", sep="\t", index=False)
        trajectory_rows, flux_rows = _write_trajectory_outputs(runs, temporary)
        _write_plot(conditions, temporary / "plots")
        _write_report(
            temporary / "BIOLOGICAL_INTERPRETATION.md", runs, conditions,
            evidence, context_sha, spec,
        )
        numerical = {
            "scope": "runtime_only",
            "calibration_status": "sensitivity_only_not_calibrated",
            "all_checks_pass": True,
            "cross_mode_context_sha256": context_sha,
            "runtime_mapping": runtime_mapping,
            "simulation_context_fingerprint": runs[0].manifests[0]["simulation_context"]["simulation_context_fingerprint"],
            "input_sha256": runs[0].manifests[0]["input_sha256"],
            "trajectory_rows_after_wt_dedup": trajectory_rows,
            "flux_panel_rows": flux_rows,
            "runs": {run.mode: run.checks for run in runs},
            "tolerances": {
                "state_and_inventory_abs": STATE_TOL,
                "solver_coupling_abs": SOLVER_TOL,
                "source_free_growth": SOURCE_FREE_GROWTH_TOL,
            },
        }
        (temporary / "numerical_checks.json").write_text(
            json.dumps(numerical, indent=2, sort_keys=True) + "\n", encoding="utf-8",
        )
        run_manifest = {
            **provenance,
            "scope": "runtime_only",
            "calibration_status": "sensitivity_only_not_calibrated",
            "cross_mode_context_sha256": context_sha,
            "runtime_mapping": runtime_mapping,
            "runs": {
                run.mode: {
                    "run_id": run.manifests[0]["run_id"],
                    "run_directory": str(run.path),
                    "slurm_jobs": validated_job_ids.get(run.mode),
                    "chunk_manifest_sha256": [
                        sha256_file(run.path / f"chunk_{index:03d}_manifest.json")
                        for index in range(spec.chunk_count)
                    ],
                    "merged_calls_sha256": sha256_file(run.path / "essentiality_dynamic_calls.tsv"),
                }
                for run in runs
            },
            "simulation_context": runs[0].manifests[0]["simulation_context"],
            "input_sha256": runs[0].manifests[0]["input_sha256"],
            "grid": {
                "genes": list(spec.panel), "alphas_mmol_gDW": list(spec.alphas),
                "pool_multipliers": list(spec.pools), "uracil_modes": list(MODE_PREFIXES),
                "hours": spec.hours, "dt_h": spec.dt,
                "initial_biomass_gDW_L": spec.initial_biomass,
            },
            "evidence_inputs": {
                "gene_evidence_matrix_sha256": sha256_file(repo / "docs/research/quinone_gpr_synthome_2026-08-17/gene_evidence_matrix.tsv"),
                "fitness_threshold_matrix_sha256": sha256_file(repo / "docs/research/quinone_gpr_synthome_2026-08-17/fitness_threshold_matrix.tsv"),
            },
            "parameter_selection": "none; full frozen grid reported without recall tuning",
        }
        (temporary / "run_manifest.json").write_text(
            json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8",
        )
        _write_sha256s(temporary)
        with zipfile.ZipFile(temporary_zip, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(item for item in temporary.rglob("*") if item.is_file()):
                archive.write(path, Path(output_dir.name) / path.relative_to(temporary))
        temporary_zip_sha.write_text(
            f"{sha256_file(temporary_zip)}  {zip_path.name}\n", encoding="utf-8",
        )
        temporary.replace(output_dir)
        temporary_zip.replace(zip_path)
        temporary_zip_sha.replace(zip_sha_path)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        temporary_zip.unlink(missing_ok=True)
        temporary_zip_sha.unlink(missing_ok=True)
        shutil.rmtree(output_dir, ignore_errors=True)
        zip_path.unlink(missing_ok=True)
        zip_sha_path.unlink(missing_ok=True)
        raise
    return output_dir, zip_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--finite-run-dir", type=Path, required=True)
    parser.add_argument("--nonlimiting-run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--compute-commit", required=True)
    parser.add_argument("--data-sha256", required=True)
    parser.add_argument("--finite-array-job-id", required=True)
    parser.add_argument("--finite-merge-job-id", required=True)
    parser.add_argument("--nonlimiting-array-job-id", required=True)
    parser.add_argument("--nonlimiting-merge-job-id", required=True)
    args = parser.parse_args(argv)
    output, archive = build_bundle(
        args.finite_run_dir, args.nonlimiting_run_dir, args.output_dir,
        compute_commit=args.compute_commit, data_sha256=args.data_sha256,
        job_ids={
            "finite_batch": {
                "array": args.finite_array_job_id, "merge": args.finite_merge_job_id,
            },
            "po1f_nonlimiting": {
                "array": args.nonlimiting_array_job_id,
                "merge": args.nonlimiting_merge_job_id,
            },
        },
    )
    print(output)
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
