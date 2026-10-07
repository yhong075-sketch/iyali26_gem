#!/usr/bin/env python3
"""User-authorized E5 CoQ9 dilution sensitivity: 21 alpha values + one control."""

import argparse
import copy
import csv
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cobra.util.solver import linear_reaction_coefficients
from scripts.build_coq_biomass_candidate import (
    BIOMASS, Q9, SPEC_PATH, apply_coq_biomass, pool_balance,
)
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.execution import execution_limits
from scripts.validate_energy_candidates import SolverBudget, load_model, sha, signature, table, write

ALPHAS = [0.0] + [10 ** (-4 + i / 10) for i in range(21)]
ALPHA_SOURCE = "User-authorized 2026-10-05 sensitivity range, 1e-4 to 1e-2 mmol/gDW; not a measured abundance"


def plot(rows, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter

    baseline = rows[0]["growth_h_inverse"]
    x = [r["alpha_mmol_per_gDW"] for r in rows[1:]]
    y = [r["growth_h_inverse"] if r["status"] == "optimal" else float("nan") for r in rows[1:]]
    loss = [100 * (1 - mu / baseline) for mu in y]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlabel(r"CoQ9 biomass requirement $\alpha$ (mmol/gDW)")
        ax.grid(True, which="major", alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].plot(x, y, "o-", color="#177c87", markersize=4, label="CoQ9-coupled E5")
    axes[0].axhline(baseline, color="#6d737b", linestyle="--", label=r"Original E5 ($\alpha=0$)")
    axes[0].set_ylim(0, max([baseline] + [v for v in y if v == v]) * 1.08)
    axes[0].set_ylabel(r"Maximum growth rate $\mu$ (h$^{-1}$)")
    axes[0].set_title("Growth response")
    axes[0].legend(frameon=False, loc="lower left", fontsize=9)
    axes[1].plot(x, loss, "o-", color="#b4653e", markersize=4)
    axes[1].axhline(0, color="#6d737b", linestyle="--", linewidth=1)
    axes[1].set_ylabel("Growth reduction vs. original E5 (%)")
    axes[1].set_title("Change from the control")
    axes[1].yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    fig.suptitle("E5 | CoQ9 growth-dilution sensitivity", fontsize=15)
    fig.supxlabel("PO1f / static SD-Leu; 21 tested values + control. Alpha values are uncalibrated assumptions.", fontsize=9)
    fig.savefig(output / "growth_curve.png", dpi=200)
    fig.savefig(output / "growth_curve.pdf")
    plt.close(fig)
    return {"python": sys.version, "matplotlib": matplotlib.__version__, "code_sha256": sha(__file__)}


def run(output, reuse_control=None):
    output = (ROOT / output).resolve()
    if not output.is_relative_to(ROOT):
        raise ValueError("Output must remain in the project workspace")
    output.mkdir(parents=True, exist_ok=False)
    spec = json.loads(SPEC_PATH.read_text())
    prior_config = ROOT / "artifacts/atp_candidate_repair_20260924/config.json"
    config = json.loads(prior_config.read_text())
    config.update(model=spec["source_path"], model_sha256=spec["source_sha256"])
    calls = len(ALPHAS) - int(reuse_control is not None)
    config["limits"] = {"solves": calls, "solve_wall_seconds": calls * 60, "per_solve_seconds": 60}
    assert sha(ROOT / config["model"]) == config["model_sha256"]
    paths = [ROOT / config[k] for k in ("model", "media", "strain_profile")]
    paths += [SPEC_PATH, prior_config, Path(__file__).resolve()]
    if reuse_control is not None:
        reuse_control = (ROOT / reuse_control).resolve()
        if not reuse_control.is_relative_to(ROOT):
            raise ValueError("Control must remain in the project workspace")
        paths.append(reuse_control)
    paths += [Path(m.__file__).resolve() for m in list(sys.modules.values())
              if getattr(m, "__file__", None) and Path(m.__file__).resolve().is_relative_to(ROOT / "scripts")]
    identities = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    manifest = {
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Growth sensitivity only; no KO, FVA, energy retest, dFBA or default model change",
        "config": config, "alpha_source": ALPHA_SOURCE, "alphas_mmol_per_gDW": ALPHAS,
        "input_and_code_sha256": identities,
        "software": {"python": sys.version, **{p: version(p) for p in ("cobra", "optlang", "gurobipy")}},
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "initial_git_dirty": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True).splitlines(),
        "historical_environment_reconstructed": False,
        "reused_control": str(reuse_control) if reuse_control is not None else None,
    }
    write(output / "manifest.json", manifest)
    rows = []
    try:
        sim = load_model(config["model"], config)
        model = sim.model
        assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {BIOMASS: 1.0}
        assert pool_balance(model) == {"R385": 1.0}
        assert model.reactions.get_by_id("R385").bounds == (0.0, 1000.0)
        before = model_definition(model)
        bio = model.reactions.get_by_id(BIOMASS)
        notes = copy.deepcopy(bio.notes)
        manifest.update(runtime_context=sim.provenance(), active_medium=sim.active_medium,
                        strain_overlay_audit=sim.strain_overlay_audit,
                        NGAM_bounds=list(model.reactions.get_by_id(config["maintenance"]).bounds))
        budget = SolverBudget(output / "solver_budget.json", config)
        for i, alpha in enumerate(ALPHAS):
            with model:
                try:
                    expected = copy.deepcopy(before)
                    if alpha:
                        apply_coq_biomass(model, alpha, ALPHA_SOURCE)
                        expected["reactions"][BIOMASS]["stoichiometry"][Q9] = -alpha
                    assert model_definition(model) == expected
                    if i == 0 and reuse_control is not None:
                        saved = json.loads(reuse_control.read_text())
                        assert saved["problem"] == signature(model), "Saved control optimization problem differs"
                        result, flux = saved["result"], saved["fluxes"]
                        assert result["status"] == "optimal" and result["postprocessing"] == "complete"
                        assert max(result[k] for k in ("max_mass_residual", "max_bound_violation", "max_constraint_violation")) <= config["tolerance"]
                        write(output / "point_00.json", saved)
                    else:
                        result, flux = budget.solve(model, f"point_{i:02d}", output / f"point_{i:02d}.json",
                            {"kind": "maximum_growth", "alpha_mmol_per_gDW": alpha})
                    mu = flux[BIOMASS] if flux is not None else None
                    if i == 0 and (result["status"] != "optimal" or mu <= config["tolerance"]):
                        raise RuntimeError("Positive optimal control required; do not normalize this run")
                    residual = flux["R385"] - alpha * mu if flux is not None else None
                    if flux is not None:
                        assert abs(mu - result["objective"]) <= config["tolerance"]
                        assert abs(residual) <= 2 * config["tolerance"]
                    baseline = rows[0]["growth_h_inverse"] if rows else mu
                    rows.append({"alpha_mmol_per_gDW": alpha, "status": result["status"],
                        "growth_h_inverse": mu, "growth_percent_control": 100 * mu / baseline if mu is not None else None,
                        "R385_flux": flux["R385"] if flux is not None else None,
                        "alpha_times_growth": alpha * mu if mu is not None else None,
                        "pool_balance_residual": residual})
                    table(output / "growth_curve.tsv", rows)
                finally:
                    bio.notes = copy.deepcopy(notes)
            assert model_definition(model) == before
        assert all(sha(ROOT / p) == h for p, h in identities.items())
        assert len(budget.record["calls"]) == calls and len(rows) == len(ALPHAS) == 22
        manifest.update(status="complete", source_and_code_unchanged=True,
                        runtime_restored=True, actual_optimization_calls=calls,
                        max_pool_residual=max(abs(r["pool_balance_residual"]) for r in rows if r["pool_balance_residual"] is not None))
    except BaseException as exc:
        manifest.update(status="failed", error=repr(exc))
        raise
    finally:
        manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
        write(output / "manifest.json", manifest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New project output directory")
    parser.add_argument("--plot-only", action="store_true", help="Plot saved results; never run optimization")
    parser.add_argument("--reuse-control", type=Path, help="Reuse a saved optimal alpha=0 result only after exact problem comparison")
    args = parser.parse_args()
    with execution_limits(allow_network=False):
        if args.plot_only:
            out = (ROOT / args.output).resolve()
            if not out.is_relative_to(ROOT):
                raise ValueError("Output must remain in the project workspace")
            manifest = json.loads((out / "manifest.json").read_text())
            assert manifest["status"] == "complete"
            assert manifest["input_and_code_sha256"][str(Path(__file__).resolve().relative_to(ROOT))] == sha(__file__)
            with (out / "growth_curve.tsv").open() as f:
                rows = [{k: v if k == "status" else float(v) if v else None for k, v in row.items()}
                        for row in csv.DictReader(f, delimiter="\t")]
            with execution_limits(no_solve=True, allow_network=False):
                manifest["plot_runtime"] = plot(rows, out)
            write(out / "manifest.json", manifest)
        else:
            run(args.output, args.reuse_control)
