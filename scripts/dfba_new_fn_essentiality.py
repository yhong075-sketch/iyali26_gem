#!/usr/bin/env python3
"""Compare baseline and candidate dFBA essentiality for experimental essentials.

The dynamic-medium CSV must contain:

    reaction_id,compound,initial_concentration_mmol_l,pool_mode,
    initial_concentration_status,max_uptake_mmol_gdw_h,
    uptake_evidence_status,uptake_basis,source_locator,source_accessed_on,notes

A blank initial concentration means the exchange is available but is not
depleted as a finite extracellular pool.  It never means zero concentration.

Only experimentally essential genes are simulated.  A "new FN" is a gene
that is essential in the supplied experiment, essential in the baseline dFBA
prediction, and non-essential in the candidate dFBA prediction.

The Ramesh et al. consensus file is a positive-only reference.  Genes absent
from that file are unknown, not experimentally non-essential.  Listed genes
outside either model are reported as untested rather than silently relabeled.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import cobra
from cobra.flux_analysis import pfba
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients


REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from scripts.gem_annotate.validate_essential_genes import load_experimental  # noqa: E402
from scripts.lp_sn12_candidate import source_fingerprint  # noqa: E402


MEDIUM_COLUMNS = (
    "reaction_id",
    "compound",
    "initial_concentration_mmol_l",
    "pool_mode",
    "initial_concentration_status",
    "max_uptake_mmol_gdw_h",
    "uptake_evidence_status",
    "uptake_basis",
    "source_locator",
    "source_accessed_on",
    "notes",
)
POOL_MODE_INITIAL_STATUSES = {
    "finite": {"nominal_formulation", "formulation_proxy"},
    "closed": {"formulation_absent"},
    "nondepleting": {"unresolved"},
    "boundary": {"environmental_boundary"},
}
POOL_MODE_UPTAKE_STATUSES = {
    "finite": {"inferred_upper_bound", "model_default_upper_bound"},
    "closed": {"closed"},
    "nondepleting": {"permissive_upper_bound"},
    "boundary": {"permissive_upper_bound"},
}
UPTAKE_BASIS_BY_STATUS = {
    "inferred_upper_bound": "initial_concentration_mmol_l_times_10_divided_by_111_rounded",
    "model_default_upper_bound": "model_default_static_bound_not_measured_kinetics",
    "permissive_upper_bound": "permissive_sensitivity_assumption_not_measured_kinetics",
    "closed": "zero_by_formulation",
}
POSITIVE_ONLY_REFERENCE_COLUMNS = (
    "gene_id",
    "source_gene_id",
    "function",
    "source",
    "confidence",
)
POSITIVE_ONLY_REFERENCE_SOURCE = "https://doi.org/10.1038/s42003-023-04996-8"
POSITIVE_ONLY_REFERENCE_CONFIDENCE = "consensus_essential_in_at_least_2_of_3_screens"
SCHEMA_VERSION = 4
ALGORITHM = "first-order Euler dFBA with parsimonious FBA flux selection"
RESULT_COLUMNS = (
    "gene_id",
    "baseline_ko_biomass_gain_gdw_l",
    "candidate_ko_biomass_gain_gdw_l",
    "baseline_ko_to_wt_gain_ratio",
    "candidate_ko_to_wt_gain_ratio",
    "baseline_predicted_essential",
    "candidate_predicted_essential",
    "new_false_negative",
    "baseline_reaction_ids",
    "candidate_reaction_ids",
    "gpr_evidence_status",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def git_head() -> str:
    return subprocess.run(
        ["git", "-C", str(REPOSITORY), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def input_records(args: argparse.Namespace) -> dict:
    paths = {
        "runner_script": Path(__file__),
        "essentiality_loader": REPOSITORY / "scripts/gem_annotate/validate_essential_genes.py",
        "fingerprint_helper": REPOSITORY / "scripts/lp_sn12_candidate.py",
        "baseline_model": args.baseline,
        "candidate_model": args.candidate,
        "experimental": args.experimental,
        "dynamic_medium": args.dynamic_medium,
    }
    return {
        name: {"path": str(path.resolve()), "sha256": sha256(path)}
        for name, path in paths.items()
    }


def _id_digest(values: list[str]) -> str:
    payload = "".join(f"{value}\n" for value in sorted(values)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def shard_gene_ids(
    gene_ids: list[str], *, shard_index: int, shard_count: int
) -> list[str]:
    """Return one deterministic, non-empty slice of a unique gene universe."""
    ordered = sorted(gene_ids)
    if not ordered or len(ordered) != len(set(ordered)):
        raise ValueError("gene shard universe must be non-empty and unique")
    if shard_count < 1 or shard_count > len(ordered):
        raise ValueError("shard_count must be in [1, number of genes]")
    if not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must be in [0, shard_count)")
    return ordered[shard_index::shard_count]


def _strict_bool(value: object, *, field: str) -> bool:
    if type(value) is bool:
        return value
    if type(value) is str and value in {"True", "False"}:
        return value == "True"
    raise ValueError(f"{field} must be exactly True or False")


def load_experimental_reference(path: Path) -> tuple[list[str], dict]:
    """Load labeled data or the Ramesh consensus positive-only reference."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = tuple(field.strip() for field in (reader.fieldnames or []))
        if "essential" in {field.lower() for field in fields}:
            experimental = load_experimental(path)
            positive_ids = sorted(
                experimental.loc[experimental["essential"], "gene_id"].tolist()
            )
            if not positive_ids:
                raise ValueError("experimental CSV contains no essential genes")
            return positive_ids, {
                "mode": "labeled",
                "positive_only": False,
                "row_count": len(experimental),
                "positive_gene_count": len(positive_ids),
                "negative_label_count": int((~experimental["essential"]).sum()),
                "positive_gene_id_sha256": _id_digest(positive_ids),
                "unlisted_gene_semantics": "unknown_not_nonessential",
            }
        if fields != POSITIVE_ONLY_REFERENCE_COLUMNS:
            raise ValueError(
                "experimental reference must have an essential column or the exact "
                f"positive-only columns: {list(POSITIVE_ONLY_REFERENCE_COLUMNS)}"
            )
        reader.fieldnames = list(fields)
        rows = []
        seen = set()
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise ValueError(f"extra positive-only reference value at line {line_number}")
            row = {
                str(key).strip(): "" if value is None else str(value).strip()
                for key, value in raw.items()
            }
            if any(not row[column] for column in POSITIVE_ONLY_REFERENCE_COLUMNS):
                raise ValueError(f"blank positive-only reference value at line {line_number}")
            gene_id = row["gene_id"]
            if gene_id in seen:
                raise ValueError(f"duplicate positive-only gene_id at line {line_number}: {gene_id}")
            if row["source_gene_id"].replace("YALI1_", "YALI1", 1) != gene_id:
                raise ValueError(f"gene ID normalization mismatch at line {line_number}: {gene_id}")
            if row["source"] != POSITIVE_ONLY_REFERENCE_SOURCE:
                raise ValueError(f"unexpected positive-only source at line {line_number}")
            if row["confidence"] != POSITIVE_ONLY_REFERENCE_CONFIDENCE:
                raise ValueError(f"unexpected positive-only confidence at line {line_number}")
            seen.add(gene_id)
            rows.append(row)
    if not rows:
        raise ValueError("positive-only experimental reference is empty")
    positive_ids = sorted(seen)
    return positive_ids, {
        "mode": "positive_only_consensus",
        "positive_only": True,
        "row_count": len(rows),
        "positive_gene_count": len(positive_ids),
        "negative_label_count": 0,
        "positive_gene_id_sha256": _id_digest(positive_ids),
        "source": POSITIVE_ONLY_REFERENCE_SOURCE,
        "confidence_definition": POSITIVE_ONLY_REFERENCE_CONFIDENCE,
        "listed_gene_semantics": "experimentally_essential",
        "unlisted_gene_semantics": "unknown_not_nonessential",
        "source_gene_id_normalization": "remove_the_underscore_after_YALI1",
    }


def experimental_reference_coverage(
    positive_gene_ids: list[str], baseline, candidate
) -> tuple[list[str], dict]:
    """Partition positive genes without turning unlisted genes into negatives."""
    reference = set(positive_gene_ids)
    baseline_ids = {gene.id for gene in baseline.genes}
    candidate_ids = {gene.id for gene in candidate.genes}
    in_both = sorted(reference & baseline_ids & candidate_ids)
    baseline_only = sorted((reference & baseline_ids) - candidate_ids)
    candidate_only = sorted((reference & candidate_ids) - baseline_ids)
    in_neither = sorted(reference - baseline_ids - candidate_ids)
    model_common_unlisted = sorted((baseline_ids & candidate_ids) - reference)
    if not in_both:
        raise ValueError("no positive-only reference genes are jointly testable")
    partitions = {
        "reference_gene_ids_in_both_models": in_both,
        "reference_gene_ids_in_baseline_only": baseline_only,
        "reference_gene_ids_in_candidate_only": candidate_only,
        "reference_gene_ids_in_neither_model": in_neither,
        "model_common_gene_ids_absent_from_reference": model_common_unlisted,
    }
    return in_both, {
        "reference_positive_gene_count": len(reference),
        "jointly_testable_positive_gene_count": len(in_both),
        "untested_positive_gene_count": len(reference) - len(in_both),
        "model_common_unlisted_gene_count": len(model_common_unlisted),
        "absence_semantics": "unknown_not_nonessential",
        "classification_scope": "jointly_testable_positive_genes_only",
        "partitions": {
            name: {
                "count": len(values),
                "gene_id_sha256": _id_digest(values),
                "gene_ids": values,
            }
            for name, values in partitions.items()
        },
    }


def load_dynamic_medium(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = tuple(field.strip() for field in (reader.fieldnames or []))
        if fields != MEDIUM_COLUMNS:
            raise ValueError(
                "dynamic-medium CSV must have the exact columns in order: "
                f"{list(MEDIUM_COLUMNS)}"
            )
        reader.fieldnames = list(fields)

        rows = []
        seen = set()
        for line_number, raw in enumerate(reader, start=2):
            if None in raw:
                raise ValueError(f"extra dynamic-medium value at line {line_number}")
            row = {
                str(key).strip(): "" if value is None else str(value).strip()
                for key, value in raw.items()
            }
            reaction_id = row["reaction_id"]
            if not reaction_id or reaction_id in seen:
                raise ValueError(f"invalid or duplicate reaction_id at line {line_number}: {reaction_id!r}")
            try:
                concentration = (
                    float(row["initial_concentration_mmol_l"])
                    if row["initial_concentration_mmol_l"]
                    else None
                )
                uptake = float(row["max_uptake_mmol_gdw_h"])
            except ValueError as error:
                raise ValueError(f"non-numeric dynamic-medium value at line {line_number}") from error
            if concentration is not None and (not math.isfinite(concentration) or concentration < 0):
                raise ValueError(
                    f"dynamic-medium concentration must be finite and non-negative at line {line_number}"
                )
            if not math.isfinite(uptake) or uptake < 0:
                raise ValueError(f"dynamic-medium values must be finite and non-negative at line {line_number}")
            pool_mode = row["pool_mode"]
            initial_status = row["initial_concentration_status"]
            uptake_status = row["uptake_evidence_status"]
            if pool_mode not in POOL_MODE_INITIAL_STATUSES:
                raise ValueError(f"invalid pool_mode at line {line_number}: {pool_mode!r}")
            if initial_status not in POOL_MODE_INITIAL_STATUSES[pool_mode]:
                raise ValueError(
                    f"initial_concentration_status conflicts with pool_mode at line {line_number}"
                )
            if uptake_status not in POOL_MODE_UPTAKE_STATUSES[pool_mode]:
                raise ValueError(
                    f"uptake_evidence_status conflicts with pool_mode at line {line_number}"
                )
            if not row["compound"] or not row["source_locator"]:
                raise ValueError(f"compound and source_locator are required at line {line_number}")
            if row["uptake_basis"] != UPTAKE_BASIS_BY_STATUS[uptake_status]:
                raise ValueError(f"uptake_basis conflicts with evidence status at line {line_number}")
            try:
                datetime.strptime(row["source_accessed_on"], "%Y-%m-%d")
            except ValueError as error:
                raise ValueError(f"invalid source_accessed_on at line {line_number}") from error
            if pool_mode in {"finite", "closed"} and concentration is None:
                raise ValueError(f"{pool_mode} pool requires a concentration at line {line_number}")
            if pool_mode in {"nondepleting", "boundary"} and concentration is not None:
                raise ValueError(f"{pool_mode} pool requires a blank concentration at line {line_number}")
            if pool_mode == "closed" and not (
                concentration == 0 and uptake == 0 and uptake_status == "closed"
            ):
                raise ValueError(f"closed pool must have zero concentration and uptake at line {line_number}")
            seen.add(reaction_id)
            rows.append({
                "reaction_id": reaction_id,
                "compound": row["compound"],
                "initial_concentration_mmol_l": concentration,
                "pool_mode": pool_mode,
                "initial_concentration_status": initial_status,
                "max_uptake_mmol_gdw_h": uptake,
                "uptake_evidence_status": uptake_status,
                "uptake_basis": row["uptake_basis"],
                "source_locator": row["source_locator"],
                "source_accessed_on": row["source_accessed_on"],
                "notes": row["notes"],
            })
    if not rows or not any(row["max_uptake_mmol_gdw_h"] > 0 for row in rows):
        raise ValueError("dynamic-medium CSV must contain at least one positive uptake limit")
    return rows


def _medium_contract(medium: list[dict]) -> dict:
    return {
        "row_count": len(medium),
        "pool_mode_reaction_ids": {
            mode: sorted(row["reaction_id"] for row in medium if row["pool_mode"] == mode)
            for mode in sorted(POOL_MODE_INITIAL_STATUSES)
        },
        "uptake_evidence_status_by_reaction": {
            row["reaction_id"]: row["uptake_evidence_status"]
            for row in sorted(medium, key=lambda item: item["reaction_id"])
        },
    }


def _exchange_details(model, medium: list[dict]) -> list[dict]:
    exchange_ids = {reaction.id for reaction in model.exchanges}
    missing = [row["reaction_id"] for row in medium if row["reaction_id"] not in exchange_ids]
    if missing:
        raise ValueError(f"dynamic-medium reactions are not model exchanges: {missing}")

    details = []
    metabolite_ids = set()
    for row in medium:
        reaction = model.reactions.get_by_id(row["reaction_id"])
        if len(reaction.metabolites) != 1:
            raise ValueError(f"exchange {reaction.id} must contain exactly one metabolite")
        coefficient = next(iter(reaction.metabolites.values()))
        metabolite_id = next(iter(reaction.metabolites)).id
        if coefficient not in {-1.0, 1.0}:
            raise ValueError(f"exchange {reaction.id} must have coefficient -1 or +1")
        if metabolite_id in metabolite_ids:
            raise ValueError(f"multiple dynamic-medium reactions share metabolite {metabolite_id}")
        metabolite_ids.add(metabolite_id)
        details.append({
            **row,
            "reaction": reaction,
            "metabolite_id": metabolite_id,
            "coefficient": coefficient,
        })
    return details


def _medium_signature(model, medium: list[dict]) -> list[tuple[str, str, float]]:
    return [
        (item["reaction_id"], item["metabolite_id"], item["coefficient"])
        for item in _exchange_details(model, medium)
    ]


def _growth_objective(model) -> dict:
    coefficients = linear_reaction_coefficients(model)
    observed = {reaction.id: coefficient for reaction, coefficient in coefficients.items()}
    if model.objective.direction != "max" or observed != {"biomass_C": 1.0}:
        raise ValueError(
            "dFBA requires the singleton maximization objective biomass_C with coefficient +1"
        )
    return coefficients


def _model_contract(model, medium: list[dict]) -> dict:
    _growth_objective(model)
    return {
        "semantic_fingerprint": source_fingerprint(model),
        "metabolite_count": len(model.metabolites),
        "reaction_count": len(model.reactions),
        "gene_count": len(model.genes),
        "objective": "maximize biomass_C",
        "dynamic_medium_signature": [
            {"reaction_id": reaction_id, "metabolite_id": metabolite_id, "coefficient": coefficient}
            for reaction_id, metabolite_id, coefficient in _medium_signature(model, medium)
        ],
    }


def simulate_dfba(
    model,
    medium: list[dict],
    *,
    hours: float,
    step_hours: float,
    initial_biomass_gdw_l: float,
) -> dict:
    """Run first-order Euler dFBA without mutating the caller's model."""
    if not all(
        math.isfinite(value) and value > 0
        for value in (hours, step_hours, initial_biomass_gdw_l)
    ):
        raise ValueError("hours, step-hours, and initial biomass must be finite and positive")
    concentrations = {
        row["reaction_id"]: row["initial_concentration_mmol_l"]
        for row in medium
        if row["initial_concentration_mmol_l"] is not None
    }
    biomass = initial_biomass_gdw_l
    time_hours = 0.0
    steps = 0

    with model:
        objective = _growth_objective(model)
        model.medium = {
            row["reaction_id"]: row["max_uptake_mmol_gdw_h"] for row in medium
        }
        exchanges = _exchange_details(model, medium)

        while time_hours < hours - 1e-12:
            step = min(step_hours, hours - time_hours)
            for item in exchanges:
                reaction = item["reaction"]
                coefficient = item["coefficient"]
                uptake = item["max_uptake_mmol_gdw_h"]
                if reaction.id in concentrations:
                    available = concentrations[reaction.id] / (biomass * step * abs(coefficient))
                    uptake = min(uptake, available)
                if coefficient < 0:
                    reaction.lower_bound = -uptake
                else:
                    reaction.upper_bound = uptake

            solution = pfba(model)
            if solution.status != "optimal":
                raise RuntimeError(f"dFBA solve failed at t={time_hours:g} h: {solution.status}")
            growth = float(sum(
                coefficient * solution.fluxes[reaction.id]
                for reaction, coefficient in objective.items()
            ))
            if not math.isfinite(growth) or growth < -1e-9:
                raise RuntimeError(f"invalid dFBA growth rate at t={time_hours:g} h: {growth}")
            growth = max(0.0, growth)

            for item in exchanges:
                reaction = item["reaction"]
                if reaction.id not in concentrations:
                    continue
                flux = float(solution.fluxes[reaction.id])
                next_concentration = (
                    concentrations[reaction.id]
                    - item["coefficient"] * flux * biomass * step
                )
                if not math.isfinite(next_concentration) or next_concentration < -1e-8:
                    raise RuntimeError(
                        f"negative/non-finite concentration for {reaction.id} at t={time_hours:g} h"
                    )
                concentrations[reaction.id] = max(0.0, next_concentration)

            biomass += growth * biomass * step
            if not math.isfinite(biomass) or biomass <= 0:
                raise RuntimeError(f"invalid biomass at t={time_hours + step:g} h: {biomass}")
            time_hours += step
            steps += 1

    return {
        "status": "optimal",
        "steps": steps,
        "final_biomass_gdw_l": biomass,
        "biomass_gain_gdw_l": biomass - initial_biomass_gdw_l,
        "final_concentrations_mmol_l": concentrations,
    }


def _knockout_simulation(model, gene_id: str, medium: list[dict], settings: dict) -> dict:
    with model:
        model.genes.get_by_id(gene_id).knock_out()
        return simulate_dfba(model, medium, **settings)


def compare_models(
    baseline,
    candidate,
    gene_ids: list[str],
    medium: list[dict],
    settings: dict,
    growth_cutoff: float,
    output_tsv: Path,
) -> tuple[list[dict], dict, dict]:
    if not math.isfinite(growth_cutoff) or not 0 < growth_cutoff <= 1:
        raise ValueError("growth_cutoff must be finite and in (0, 1]")
    if output_tsv.exists():
        raise FileExistsError(output_tsv)
    missing_baseline = sorted(set(gene_ids) - {gene.id for gene in baseline.genes})
    missing_candidate = sorted(set(gene_ids) - {gene.id for gene in candidate.genes})
    if missing_baseline or missing_candidate:
        raise ValueError(
            "experimental essential genes are missing from a model; "
            f"baseline={missing_baseline}, candidate={missing_candidate}"
        )
    baseline_medium = _medium_signature(baseline, medium)
    candidate_medium = _medium_signature(candidate, medium)
    if baseline_medium != candidate_medium:
        raise ValueError(
            f"baseline/candidate dynamic-medium exchange tuples differ: "
            f"baseline={baseline_medium}, candidate={candidate_medium}"
        )
    _growth_objective(baseline)
    _growth_objective(candidate)

    baseline_wt = simulate_dfba(baseline, medium, **settings)
    candidate_wt = simulate_dfba(candidate, medium, **settings)
    if baseline_wt["biomass_gain_gdw_l"] <= 1e-12 or candidate_wt["biomass_gain_gdw_l"] <= 1e-12:
        raise RuntimeError("wild-type dFBA biomass gain must be positive in both models")

    rows = []
    with output_tsv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for index, gene_id in enumerate(sorted(gene_ids), start=1):
            baseline_ko = _knockout_simulation(baseline, gene_id, medium, settings)
            candidate_ko = _knockout_simulation(candidate, gene_id, medium, settings)
            baseline_ratio = baseline_ko["biomass_gain_gdw_l"] / baseline_wt["biomass_gain_gdw_l"]
            candidate_ratio = candidate_ko["biomass_gain_gdw_l"] / candidate_wt["biomass_gain_gdw_l"]
            baseline_essential = baseline_ratio < growth_cutoff
            candidate_essential = candidate_ratio < growth_cutoff
            row = {
                "gene_id": gene_id,
                "baseline_ko_biomass_gain_gdw_l": baseline_ko["biomass_gain_gdw_l"],
                "candidate_ko_biomass_gain_gdw_l": candidate_ko["biomass_gain_gdw_l"],
                "baseline_ko_to_wt_gain_ratio": baseline_ratio,
                "candidate_ko_to_wt_gain_ratio": candidate_ratio,
                "baseline_predicted_essential": baseline_essential,
                "candidate_predicted_essential": candidate_essential,
                "new_false_negative": baseline_essential and not candidate_essential,
                "baseline_reaction_ids": ";".join(sorted(
                    reaction.id for reaction in baseline.genes.get_by_id(gene_id).reactions
                )),
                "candidate_reaction_ids": ";".join(sorted(
                    reaction.id for reaction in candidate.genes.get_by_id(gene_id).reactions
                )),
                "gpr_evidence_status": "model/GPR assignment only",
            }
            writer.writerow(row)
            stream.flush()
            rows.append(row)
            print(f"[{index}/{len(gene_ids)}] {gene_id}: new_FN={row['new_false_negative']}", flush=True)
    return rows, baseline_wt, candidate_wt


def _validate_arguments(args: argparse.Namespace) -> None:
    for path in (args.baseline, args.candidate, args.experimental, args.dynamic_medium):
        if not path.is_file():
            raise FileNotFoundError(path)
    if not all(
        math.isfinite(value) and value > 0
        for value in (args.hours, args.step_hours, args.initial_biomass)
    ):
        raise ValueError("hours, step-hours, and initial-biomass must be positive")
    if not math.isfinite(args.growth_cutoff) or not 0 < args.growth_cutoff <= 1:
        raise ValueError("growth-cutoff must be finite and in (0, 1]")


def _slurm_record() -> dict:
    return {
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "node_list": os.environ.get("SLURM_JOB_NODELIST"),
        "partition": os.environ.get("SLURM_JOB_PARTITION"),
        "account": os.environ.get("SLURM_JOB_ACCOUNT"),
        "cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "memory_per_node": os.environ.get("SLURM_MEM_PER_NODE"),
        "time_limit": os.environ.get("SLURM_TIMELIMIT"),
    }


def _prepare_context(args: argparse.Namespace) -> dict:
    _validate_arguments(args)
    inputs = input_records(args)
    medium = load_dynamic_medium(args.dynamic_medium)
    positive_gene_ids, reference_contract = load_experimental_reference(args.experimental)
    baseline = read_sbml_model(str(args.baseline))
    candidate = read_sbml_model(str(args.candidate))
    baseline.solver = args.solver
    candidate.solver = args.solver
    gene_ids, coverage_contract = experimental_reference_coverage(
        positive_gene_ids, baseline, candidate
    )
    settings = {
        "hours": args.hours,
        "step_hours": args.step_hours,
        "initial_biomass_gdw_l": args.initial_biomass,
        "growth_cutoff": args.growth_cutoff,
        "solver": args.solver,
    }
    return {
        "git_commit": git_head(),
        "inputs": inputs,
        "medium": medium,
        "positive_gene_ids": positive_gene_ids,
        "gene_ids": gene_ids,
        "baseline": baseline,
        "candidate": candidate,
        "settings": settings,
        "model_contracts": {
            "baseline": _model_contract(baseline, medium),
            "candidate": _model_contract(candidate, medium),
        },
        "dynamic_medium_contract": _medium_contract(medium),
        "experimental_reference_contract": reference_contract,
        "experimental_reference_coverage": coverage_contract,
        "software": {
            "python": platform.python_version(),
            "cobra": cobra.__version__,
            "solver_interface": baseline.solver.interface.__name__,
        },
    }


def _shared_contract(context: dict) -> dict:
    return {
        "git_commit": context["git_commit"],
        "inputs": context["inputs"],
        "settings": context["settings"],
        "software": context["software"],
        "model_contracts": context["model_contracts"],
        "dynamic_medium_contract": context["dynamic_medium_contract"],
        "experimental_reference_contract": context["experimental_reference_contract"],
        "experimental_reference_coverage": context["experimental_reference_coverage"],
    }


def _assert_context_unchanged(context: dict, args: argparse.Namespace) -> None:
    if input_records(args) != context["inputs"] or git_head() != context["git_commit"]:
        raise RuntimeError("repository or input files changed during dFBA evaluation")


def _execution_contract(context: dict, args: argparse.Namespace) -> dict:
    selected = shard_gene_ids(
        context["gene_ids"],
        shard_index=args.shard_index,
        shard_count=args.shard_count,
    )
    deferred = args.shard_count > 1
    return {
        "mode": "shard" if deferred else "single",
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "evaluated_gene_ids": selected,
        "evaluated_gene_id_sha256": _id_digest(selected),
        "universe_gene_count": len(context["gene_ids"]),
        "universe_gene_id_sha256": _id_digest(context["gene_ids"]),
        "final_scientific_gate_deferred_to_aggregate": deferred,
    }


def _governance() -> dict:
    return {
        "model_mutation_performed": False,
        "production_gate_passed": False,
        "human_review_required": True,
        "positive_only_reference_absence_is_not_nonessential": True,
        "gene_identity_function_followup_required_for_new_false_negatives": True,
    }


def _write_results(path: Path, rows: list[dict]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def run(args: argparse.Namespace) -> dict:
    context = _prepare_context(args)
    execution = _execution_contract(context, args)
    running = json.loads((args.output_dir / "summary.json").read_text(encoding="utf-8"))
    running.update({
        **_shared_contract(context),
        "execution": execution,
        "production_gate_passed": False,
        "human_review_required": True,
    })
    write_json(args.output_dir / "summary.json", running)

    output_tsv = args.output_dir / "gene_results.tsv"
    partial_tsv = args.output_dir / "gene_results.partial.tsv"
    rows, baseline_wt, candidate_wt = compare_models(
        context["baseline"],
        context["candidate"],
        execution["evaluated_gene_ids"],
        context["medium"],
        {
            "hours": context["settings"]["hours"],
            "step_hours": context["settings"]["step_hours"],
            "initial_biomass_gdw_l": context["settings"]["initial_biomass_gdw_l"],
        },
        context["settings"]["growth_cutoff"],
        partial_tsv,
    )
    _assert_context_unchanged(context, args)
    partial_tsv.replace(output_tsv)
    new_false_negatives = sorted(row["gene_id"] for row in rows if row["new_false_negative"])
    local_gate = not new_false_negatives
    final_gate = local_gate if not execution["final_scientific_gate_deferred_to_aggregate"] else None
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "definition": (
            "new FN = experimentally essential AND baseline dFBA predicts essential "
            "AND candidate dFBA predicts non-essential, among positive-reference "
            "genes present in both models"
        ),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv,
        "algorithm": ALGORITHM,
        **_shared_contract(context),
        "slurm": _slurm_record(),
        "wild_type": {"baseline": baseline_wt, "candidate": candidate_wt},
        "execution": execution,
        "production_gate_passed": False,
        "human_review_required": True,
        "experimental_essential_gene_count": len(context["positive_gene_ids"]),
        "evaluated_gene_count": len(rows),
        "new_false_negative_count": len(new_false_negatives),
        "new_false_negative_gene_ids": new_false_negatives,
        "no_new_false_negative_gate_passed": local_gate,
        "final_scientific_gate_passed": final_gate,
        "gene_results": {"path": str(output_tsv.resolve()), "sha256": sha256(output_tsv)},
        "governance": _governance(),
    }


def _read_shard_results(path: Path, expected_gene_ids: list[str]) -> list[dict]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or []) != RESULT_COLUMNS:
            raise ValueError(f"unexpected result columns in {path}")
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or set(row) != set(RESULT_COLUMNS):
                raise ValueError(f"malformed result row at {path}:{line_number}")
            for field in (
                "baseline_predicted_essential",
                "candidate_predicted_essential",
                "new_false_negative",
            ):
                _strict_bool(row[field], field=f"{path}:{line_number}:{field}")
            expected_new = (
                _strict_bool(row["baseline_predicted_essential"], field="baseline_predicted_essential")
                and not _strict_bool(
                    row["candidate_predicted_essential"], field="candidate_predicted_essential"
                )
            )
            if _strict_bool(row["new_false_negative"], field="new_false_negative") != expected_new:
                raise ValueError(f"new-FN classification is inconsistent at {path}:{line_number}")
            rows.append(row)
    actual_gene_ids = [row["gene_id"] for row in rows]
    if actual_gene_ids != expected_gene_ids:
        raise ValueError(
            f"result gene IDs do not match the expected shard: {path}"
        )
    return rows


def _worker_summary(
    path: Path, expected_execution: dict, context: dict
) -> tuple[dict, list[dict], dict]:
    summary_path = path / "summary.json"
    results_path = path / "gene_results.tsv"
    partial_path = path / "gene_results.partial.tsv"
    wrapper_path = path / "wrapper_status.json"
    if partial_path.exists():
        raise ValueError(f"partial result remains in shard: {partial_path}")
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid shard summary: {summary_path}") from error
    if summary.get("status") != "complete":
        raise ValueError(f"shard is not complete: {path}")
    if summary.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"shard schema version differs: {path}")
    if summary.get("execution") != expected_execution:
        raise ValueError(f"shard execution contract differs: {path}")
    for key, value in _shared_contract(context).items():
        if summary.get(key) != value:
            raise ValueError(f"shard {key} differs: {path}")
    descriptor = summary.get("gene_results")
    expected_descriptor = {"path": str(results_path.resolve()), "sha256": sha256(results_path)}
    if descriptor != expected_descriptor:
        raise ValueError(f"shard result descriptor differs: {path}")
    if not wrapper_path.is_file():
        raise FileNotFoundError(wrapper_path)
    try:
        wrapper = json.loads(wrapper_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid shard wrapper status: {wrapper_path}") from error
    if (
        wrapper.get("mode") != "shard"
        or wrapper.get("array_task_id") != str(expected_execution["shard_index"])
        or wrapper.get("shard_count") != str(expected_execution["shard_count"])
        or wrapper.get("preflight_exit_code") != 0
        or wrapper.get("runner_exit_code") != 0
        or wrapper.get("postflight_exit_code") != 0
        or wrapper.get("technical_execution_completed") is not True
        or wrapper.get("final_scientific_gate_deferred_to_aggregate") is not True
        or wrapper.get("wrapper_gate_passed") is not True
        or wrapper.get("summary_sha256") != sha256(summary_path)
        or wrapper.get("gene_results_sha256") != sha256(results_path)
    ):
        raise ValueError(f"shard wrapper status is inconsistent: {path}")
    rows = _read_shard_results(results_path, expected_execution["evaluated_gene_ids"])
    new_ids = [
        row["gene_id"]
        for row in rows
        if _strict_bool(row["new_false_negative"], field="new_false_negative")
    ]
    if (
        summary.get("evaluated_gene_count") != len(rows)
        or summary.get("new_false_negative_count") != len(new_ids)
        or summary.get("new_false_negative_gene_ids") != new_ids
        or _strict_bool(
            summary.get("no_new_false_negative_gate_passed"),
            field="no_new_false_negative_gate_passed",
        ) != (not new_ids)
        or summary.get("final_scientific_gate_passed") is not None
    ):
        raise ValueError(f"shard result summary is inconsistent: {path}")
    return summary, rows, wrapper


def aggregate_shards(args: argparse.Namespace) -> dict:
    context = _prepare_context(args)
    shard_root = args.aggregate_shards.resolve()
    expected_shard_names = {f"{index:03d}" for index in range(args.shard_count)}
    observed_shard_names = {path.name for path in shard_root.iterdir()}
    missing_shard_names = expected_shard_names - observed_shard_names
    if missing_shard_names:
        raise FileNotFoundError(
            f"missing expected shard directories: {sorted(missing_shard_names)}"
        )
    if observed_shard_names != expected_shard_names:
        raise ValueError("shards directory must contain exactly the expected shard directories")
    all_rows = []
    manifests = []
    common_wild_type = None
    for shard_index in range(args.shard_count):
        expected_execution = {
            "mode": "shard",
            "shard_index": shard_index,
            "shard_count": args.shard_count,
            "evaluated_gene_ids": shard_gene_ids(
                context["gene_ids"], shard_index=shard_index, shard_count=args.shard_count
            ),
            "evaluated_gene_id_sha256": _id_digest(shard_gene_ids(
                context["gene_ids"], shard_index=shard_index, shard_count=args.shard_count
            )),
            "universe_gene_count": len(context["gene_ids"]),
            "universe_gene_id_sha256": _id_digest(context["gene_ids"]),
            "final_scientific_gate_deferred_to_aggregate": True,
        }
        shard_path = shard_root / f"{shard_index:03d}"
        summary, rows, wrapper = _worker_summary(shard_path, expected_execution, context)
        if common_wild_type is None:
            common_wild_type = summary["wild_type"]
        elif summary.get("wild_type") != common_wild_type:
            raise ValueError(f"wild-type result differs between shards: {shard_path}")
        all_rows.extend(rows)
        manifests.append({
            "shard_index": shard_index,
            "summary_path": str((shard_path / "summary.json").resolve()),
            "summary_sha256": sha256(shard_path / "summary.json"),
            "wrapper_status_path": str((shard_path / "wrapper_status.json").resolve()),
            "wrapper_status_sha256": sha256(shard_path / "wrapper_status.json"),
            "gene_results_path": str((shard_path / "gene_results.tsv").resolve()),
            "gene_results_sha256": sha256(shard_path / "gene_results.tsv"),
            "evaluated_gene_id_sha256": expected_execution["evaluated_gene_id_sha256"],
        })
    by_gene_id = {row["gene_id"]: row for row in all_rows}
    if len(by_gene_id) != len(all_rows) or sorted(by_gene_id) != context["gene_ids"]:
        raise ValueError("shards do not provide one result for every expected gene")
    ordered_rows = [by_gene_id[gene_id] for gene_id in context["gene_ids"]]
    output_tsv = args.output_dir / "gene_results.tsv"
    partial_tsv = args.output_dir / "gene_results.partial.tsv"
    _write_results(partial_tsv, ordered_rows)
    _assert_context_unchanged(context, args)
    partial_tsv.replace(output_tsv)
    new_false_negatives = [
        row["gene_id"]
        for row in ordered_rows
        if _strict_bool(row["new_false_negative"], field="new_false_negative")
    ]
    final_gate = not new_false_negatives
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "complete",
        "definition": (
            "new FN = experimentally essential AND baseline dFBA predicts essential "
            "AND candidate dFBA predicts non-essential, among positive-reference "
            "genes present in both models"
        ),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv,
        "algorithm": ALGORITHM,
        **_shared_contract(context),
        "slurm": _slurm_record(),
        "wild_type": common_wild_type,
        "execution": {
            "mode": "aggregate",
            "shard_count": args.shard_count,
            "shards_root": str(shard_root),
            "universe_gene_count": len(context["gene_ids"]),
            "universe_gene_id_sha256": _id_digest(context["gene_ids"]),
        },
        "shard_manifest": manifests,
        "production_gate_passed": False,
        "human_review_required": True,
        "experimental_essential_gene_count": len(context["positive_gene_ids"]),
        "evaluated_gene_count": len(ordered_rows),
        "new_false_negative_count": len(new_false_negatives),
        "new_false_negative_gene_ids": new_false_negatives,
        "no_new_false_negative_gate_passed": final_gate,
        "final_scientific_gate_passed": final_gate,
        "gene_results": {"path": str(output_tsv.resolve()), "sha256": sha256(output_tsv)},
        "governance": _governance(),
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--baseline", required=True, type=Path)
    result.add_argument("--candidate", required=True, type=Path)
    result.add_argument("--experimental", required=True, type=Path)
    result.add_argument("--dynamic-medium", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--hours", type=float, default=24.0)
    result.add_argument("--step-hours", type=float, default=0.1)
    result.add_argument("--initial-biomass", type=float, default=0.05)
    result.add_argument("--growth-cutoff", type=float, default=0.01)
    result.add_argument("--solver", default="glpk")
    result.add_argument("--shard-index", type=int, default=0)
    result.add_argument("--shard-count", type=int, default=1)
    result.add_argument("--aggregate-shards", type=Path)
    return result


def main() -> None:
    args = parser().parse_args()
    if args.aggregate_shards:
        if args.shard_index != 0:
            raise SystemExit("aggregate mode requires the default shard index 0")
        if args.aggregate_shards.resolve() != (args.output_dir / "shards").resolve():
            raise SystemExit("aggregate shards must be exactly output-dir/shards")
        if not args.aggregate_shards.is_dir():
            raise SystemExit(f"missing shards directory: {args.aggregate_shards}")
        if set(path.name for path in args.output_dir.iterdir()) != {"shards"}:
            raise SystemExit("aggregate output directory may contain only shards before aggregation")
    else:
        if args.output_dir.exists() and any(args.output_dir.iterdir()):
            raise SystemExit(f"output directory must be empty: {args.output_dir}")
        args.output_dir.mkdir(parents=True, exist_ok=True)
    write_json(args.output_dir / "summary.json", {
        "schema_version": SCHEMA_VERSION,
        "status": "running",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "production_gate_passed": False,
        "human_review_required": True,
    })
    try:
        summary = aggregate_shards(args) if args.aggregate_shards else run(args)
    except Exception as error:
        summary_path = args.output_dir / "summary.json"
        failure = json.loads(summary_path.read_text())
        failure.update({
            "status": "failed",
            "error_type": type(error).__name__,
            "error": str(error),
            "production_gate_passed": False,
            "human_review_required": True,
        })
        write_json(summary_path, failure)
        if args.aggregate_shards:
            raise SystemExit(2) from error
        raise
    write_json(args.output_dir / "summary.json", summary)
    print(json.dumps({
        "new_false_negative_count": summary["new_false_negative_count"],
        "no_new_false_negative_gate_passed": summary["no_new_false_negative_gate_passed"],
        "summary": str((args.output_dir / "summary.json").resolve()),
    }, sort_keys=True))
    if summary["final_scientific_gate_passed"] is False:
        raise SystemExit(3)


if __name__ == "__main__":
    main()
