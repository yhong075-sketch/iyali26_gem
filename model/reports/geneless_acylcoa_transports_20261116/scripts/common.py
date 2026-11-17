"""Shared helpers for the gene-less CoA-ester transport review (in memory only; no SBML file is ever written).

Adapted from model/reports/r1172_side_effect_20261116/scripts/common.py (worktree bold-bartik-794e7e).

Conditions
- "open": raw model, every exchange at (-1000, 1000) (the gap census definition).
- "taskA": Task A medium + PO1f profile through the platform's load_effective_simulation_context.

Solvers
- "gurobi": feasibility and optimality tolerance 1e-9, solution and basis discarded before every LP.
- "glpk": basis reset with glp_adv_basis before every LP (stale-basis guard).

Every LP goes through lp(), which counts it and returns the raw status next to the value. A non-optimal status
returns value None; callers must never read None as zero. The LP count of each script is written to its output.
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
LIST = RR / "state/essentiality/repository/consensus_essential_genes.csv"
LIST_SHA = "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b"
ART = RR / "artifacts/tasks/geneless_acylcoa_transports_20261116"

# The eight reactions of the task.
EIGHT = ["R2025", "R2113", "R2077", "R2083", "R2084", "R2118", "R1291", "R2117"]
# Neighbours: the reactions sharing a crossing species, the branches they serve, and the matrix CoA reactions.
NEIGH = ["R2024", "R2023", "R2022", "R2281", "R2282", "R2283", "R2284", "R2285", "R968", "R2112", "R2074",
         "R2075", "R2076", "R1969", "R2114", "R2115", "R2116", "R2082", "R847", "R2140", "R2141", "R2142",
         "R1292", "R265", "R834", "R1971", "R1972", "R25", "R95", "R552", "R1968", "R1970"]
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
TOL = 1e-9

# Genes on the reactions that share a crossing species with one of the eight, or that make or use those species.
KO_GENES = sorted(set([
    "YALI1F08845g",  # R2024 propanoate:CoA ligase (also R90 acetyl-CoA synthetase)
    "YALI1C02550g",  # R2112 (also R2119, R2121)
    "YALI0E18568g",  # R2285 thiolase (YALI0 tag of YALI1E22238g per W29 audit)
    "YALI0B10406g",  # R2283 hydratase (YALI0 tag of YALI1B14044g per W29 audit)
    "YALI1D19252g",  # R2117, R2282
    "YALI1D11159g", "YALI1F07630g", "YALI1D26360g", "YALI1D30866g",  # BCKDH rule (R1968-R1970, R2281)
    "YALI1B14044g",  # R2116, R2141, R1978 and peroxisomal isomerases
    "YALI1F03883g",  # R2114
    "YALI1D07981g",  # R2140 (list annotation: similar to S. cerevisiae EHD3)
    "YALI1E15555g",  # R2082
    "YALI0D12573g",  # R1971 (no YALI1 match in the W29 audit)
    "YALI1B19217g",  # R1972
    "YALI1B02194g", "YALI1C04093g", "YALI1E00588g", "YALI1F31215g",  # R2074
    "YALI1E01019g", "YALI1E03300g",  # R25 (also citrate synthase)
    "YALI1F03803g",  # R95
    "YALI1C24124g", "YALI1F39620g",  # R552
    "YALI1F20788g",  # R834 CoA carrier
    "YALI1F13202g",  # R265 and R1092 dephospho-CoA kinase
    "YALI1D24649g",  # R615 cytosolic PPAT
    "YALI1D33780g",  # pantothenate kinase R616, R1908, R1909
    "YALI1A16988g",  # R847 FADH2 transport
    "YALI1F31153g",  # R1977 ETF:ubiquinone oxidoreductase
]))

LP_COUNT = {"n": 0}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check_inputs():
    for path, want in [(MODEL, MODEL_SHA), (MEDIUM, MEDIUM_SHA), (PROFILE, PROFILE_SHA), (LIST, LIST_SHA)]:
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
    """Max or min of a reaction id or of a {reaction_id: coefficient} expression."""
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
        if rid not in model.reactions:
            out[rid] = {"absent": True}
            continue
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


def min_total_flux(model):
    """One counted LP: minimise the sum of absolute fluxes under the model's current bounds.

    Same construction as cobra's add_pfba (sloppy Zero objective, then linear coefficients), which the R1172 round
    found to be the form optlang accepts. Returns (total or None, status, flux or None).
    """
    from optlang.symbolics import Zero
    with model:
        coefficients = {}
        for r in model.reactions:
            coefficients[r.forward_variable] = 1.0
            coefficients[r.reverse_variable] = 1.0
        model.objective = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective.set_linear_coefficients(coefficients)
        total, status = lp(model)
        flux = flux_vector(model) if total is not None else None
    return total, status, flux


def pfba(model, growth_rxn="biomass_C"):
    """Two counted LPs: maximise growth, then minimise total flux with growth fixed at that optimum (1e-9 slack)."""
    growth, s = lp(model)
    if growth is None:
        return None, s, None
    with model:
        model.reactions.get_by_id(growth_rxn).lower_bound = growth * (1 - 1e-9)
        total, s2, flux = min_total_flux(model)
    return growth, s2, flux


def classify(ko, wt):
    if ko is None or wt is None or wt == 0:
        return {c: None for c in CUTOFFS}
    ratio = ko / wt
    return {c: ("essential" if ratio < c else "non-essential") for c in CUTOFFS}


def energy_check(model):
    """Every boundary lower bound and every positive lower bound (xMAINTENANCE) set to 0; max R_NTP1 and biomass_C."""
    out = {}
    with model:
        for r in model.boundary:
            r.lower_bound = max(r.lower_bound, 0.0)
        for r in model.reactions:
            if r.lower_bound > 0:
                r.lower_bound = 0.0
        for obj in ("R_NTP1", "biomass_C"):
            with model:
                model.objective = obj
                model.objective_direction = "max"
                v, s = lp(model)
            out[obj] = {"value": v, "status": s}
    return out


def dump(obj, name):
    ART.joinpath("out").mkdir(parents=True, exist_ok=True)
    path = ART / "out" / name
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=True, default=str)
    return path
