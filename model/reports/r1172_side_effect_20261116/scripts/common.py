"""Shared helpers for the R1172 side-effect investigation (in memory only; the SBML file is never written).

Conditions
- "open": raw model, every exchange at (-1000, 1000) (the gap census definition).
- "taskA": Task A medium + PO1f profile through the platform's load_effective_simulation_context.

Solvers
- "gurobi": feasibility and optimality tolerance 1e-9, solution and basis discarded before every LP.
- "glpk": basis reset with glp_adv_basis before every LP (stale-basis guard).

Every LP goes through lp(), which counts it and returns the raw status next to the value. A non-optimal status
returns value None; callers must never read None as zero.
"""
import hashlib
import json
import math
import os
import sys
import warnings
from pathlib import Path

import cobra

REPO = Path(__file__).resolve().parents[4]
MODEL = REPO / "model/candidates/E5_353cefdc_xref_20261115.xml"
MODEL_SHA = "7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146"
MEDIUM = REPO / "model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
MEDIUM_SHA = "6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45"
RR = Path(os.environ["IYALI26_RESEARCH_ROOT"])
PROFILE = RR / "state/strain_profiles/po1f_sd_leu.json"
PROFILE_SHA = "35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383"
ART = RR / "artifacts/tasks/r1172_side_effect_20261116"

NINE = ["R25", "R552", "R95", "R2074", "R2075", "R2076", "R2077", "R2112", "R834"]
# Reactions with a CoA-containing species in C_mi and another outside C_mi (static scan of 7ccbb568).
BOUNDARY = ["R834", "R1172", "R2025", "R2077", "R2083", "R2084", "R2113", "R2118"]
R1172_OLD_BOUNDS = (-1000.0, 1000.0)
TOL = 1e-9

LP_COUNT = {"n": 0}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check_inputs():
    for path, want in [(MODEL, MODEL_SHA), (MEDIUM, MEDIUM_SHA), (PROFILE, PROFILE_SHA)]:
        got = sha256(path)
        if got != want:
            sys.exit(f"hash mismatch for {path}: {got} != {want}")


def load(condition, solver):
    check_inputs()
    if condition == "open":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = cobra.io.read_sbml_model(MODEL)  # Path object: no 260-character check
        for r in model.exchanges:
            r.bounds = (-1000.0, 1000.0)
    elif condition == "taskA":
        sys.path.insert(0, str(REPO / "platform"))
        from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
        ctx = load_effective_simulation_context(model_path=MODEL, media_path=MEDIUM, strain_profile_path=PROFILE)
        if ctx.canonical_model_sha256 != MODEL_SHA:
            sys.exit("context model hash mismatch")
        model = ctx.model
    else:
        raise ValueError(condition)
    configure(model, solver)
    return model


def configure(model, solver):
    model.solver = solver
    if solver == "gurobi":
        model.solver.configuration.tolerances.feasibility = 1e-9
        model.solver.configuration.tolerances.optimality = 1e-9
    model.solver.configuration.timeout = 120


def reset(model):
    iface = model.solver.interface.__name__
    if "gurobi" in iface:
        model.solver.problem.reset()
    elif "glpk" in iface:
        import swiglpk
        swiglpk.glp_adv_basis(model.solver.problem, 0)


def lp(model):
    """One counted LP on the model's current objective. Returns (value or None, raw status)."""
    reset(model)
    LP_COUNT["n"] += 1
    value = model.slim_optimize(error_value=float("nan"))
    status = model.solver.status
    if status != "optimal" or value is None or math.isnan(value):
        return None, status
    return value, status


def extreme(model, objective, sense):
    """Max or min of a reaction id or of a {reaction_id: coefficient} expression, no growth floor."""
    with model:
        if isinstance(objective, str):
            model.objective = model.reactions.get_by_id(objective)
        else:
            model.objective = {model.reactions.get_by_id(k): v for k, v in objective.items()}
        model.objective_direction = sense
        return lp(model)


def fva(model, rxn_ids):
    out = {}
    for rid in rxn_ids:
        lo, slo = extreme(model, rid, "min")
        hi, shi = extreme(model, rid, "max")
        blocked = None
        if lo is not None and hi is not None:
            blocked = abs(lo) <= TOL and abs(hi) <= TOL
        out[rid] = {"min": lo, "min_status": slo, "max": hi, "max_status": shi, "blocked": blocked}
    return out


def flux_vector(model):
    primal = model.solver.primal_values
    return {r.id: primal[r.forward_variable.name] - primal[r.reverse_variable.name] for r in model.reactions}


def dump(obj, name):
    ART.joinpath("out").mkdir(parents=True, exist_ok=True)
    path = ART / "out" / name
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True)
    return path
