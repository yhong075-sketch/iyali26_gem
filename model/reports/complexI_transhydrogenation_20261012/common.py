"""Shared setup for Task D (complex I transhydrogenation module). Diagnostic only: nothing here writes a model.

Inputs are pinned by SHA-256 and read from the tracked files at HEAD 4d97b82 (the Task B candidate, the Task A
medium and the ASSUMPTION files are on main). Medium logic (A = Task A medium + PO1f profile; V2 = W29-type
defined medium) and the caps helper follow Task C's common.py; the overlay is Task C's file, byte-identical.
"""
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "platform"))
sys.path.insert(0, str(HERE))

from cobra.flux_analysis import flux_variability_analysis, pfba  # noqa: E402
from cobra.util.solver import linear_reaction_coefficients  # noqa: E402
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context  # noqa: E402
from tools.diagnose_closed_energy import configure_solver  # noqa: E402
from holistic_overlay import apply_overlay  # noqa: E402

RESEARCH = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
ART = RESEARCH / "artifacts/tasks" / HERE.name
PINS = {
    "model": (REPO / "model/candidates/E5_coq9_resp_pkg_r132_20261012.xml",
              "a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515"),
    "medium": (REPO / "model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv",
               "6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45"),
    "caps": (REPO / "model/conditions/assumptions/w29_13c_upper_caps_20261012.json",
             "df1e33a8aa72c2d1b83129cb6a0e0d5e6eb283e3aec8feeecd844bcc2dc4257b"),
    "r825": (REPO / "model/conditions/assumptions/r825_reversible_20261013.json",
             "e45cab424c30ee9aec00823cb67392742ac6c36ba7091664249765a2533efcaa"),
    "overlay": (HERE / "holistic_overlay.py", "bf6c7ce54d8338c365d98c82ca345ff2bd8c6055b38d402e8144a16b4abf30c4"),
    "profile": (RESEARCH / "state/strain_profiles/po1f_sd_leu.json",
                "35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383"),
    "w29_fluxes": (RESEARCH / "reference/w29_worland2024/glucose_central_fluxes.csv",
                   "844615d0edfd8ffb033c1bce15ed746e9bb8ea02cbbe75016771b5c5fd08389f"),
}
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0, FeasibilityTol=1e-7, OptimalityTol=1e-7, IntFeasTol=1e-7)
UPTAKE = 2.414
AMINO_ACID_ROWS = ("R1003", "R1202", "R1204", "R1215", "R1217", "R1220", "R1222", "R1223", "R1231", "R1232",
                   "R1233", "R1234")
KO = "R2062"  # complex I
CUTOFF = 0.15
TOL = 1e-6
NODES = {"NADH_cy": "m123[C_cy]", "NADPH_cy": "m40[C_cy]", "NADH_mi": "m30[C_mi]", "NADPH_mi": "m176[C_mi]"}
CONFIGS = ("A", "V2")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_pins():
    for key, (path, digest) in PINS.items():
        assert sha(path) == digest, f"Input identity changed: {key}"
    return {k: {"path": str(p), "sha256": d} for k, (p, d) in PINS.items()}


class Counter_:
    lps = 0


def load(config):
    """config 'A': Task A medium + PO1f profile, glucose uptake <= 2.414. 'V2': W29-type defined, uracil closed,
    glucose fixed at 2.414. Both get the holistic overlay (diagnostic) and the 13C caps (module form)."""
    sim = load_effective_simulation_context(model_path=PINS["model"][0], media_path=PINS["medium"][0],
                                            strain_profile_path=PINS["profile"][0])
    model = sim.model
    model.solver = "gurobi"
    configure_solver(model, {"solver": PARAMS})
    assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {"biomass_C": 1.0}
    assert math.isclose(-model.reactions.R1070.lower_bound, UPTAKE)
    if config == "V2":
        current = dict(model.medium)
        for rid in AMINO_ACID_ROWS:
            assert rid in current, rid
            current.pop(rid)
        current.pop("R1354")
        model.medium = current
        r612 = model.reactions.R612
        assert tuple(r612.bounds) == (0.0, 0.0), "PO1f URA3 overlay not found"
        r612.bounds = (0.0, 1000.0)
        model.reactions.R1070.bounds = (-UPTAKE, -UPTAKE)
    else:
        assert config == "A", config
    overlay = apply_overlay(model)
    caps = add_caps(model, json.loads(PINS["caps"][0].read_text()))
    return model, {"overlay": overlay, "caps": caps}


def add_caps(model, caps, substrate="glucose"):
    """13C WT-state caps (ASSUMPTION, diagnostic): PC, net malate import, glycine-cleavage module form."""
    uptake = -model.reactions.get_by_id(caps["substrates"][substrate]["exchange"]).flux_expression
    chosen = [("cap_pc", "pyruvate_carboxylase"), ("cap_malate", "net_malate_import"),
              ("cap_gcv_module", "glycine_cleavage_module")]
    constraints = []
    for name, cap in chosen:
        spec = caps["caps"][cap]
        expression = sum(c * model.reactions.get_by_id(rid).flux_expression for rid, c in spec["reactions"].items())
        constraints.append(model.problem.Constraint(expression - spec["fraction"][substrate] * uptake, ub=0, name=name))
    model.add_cons_vars(constraints)
    return {name: caps["caps"][cap]["fraction"][substrate] for name, cap in chosen}


def growth(model):
    Counter_.lps += 1
    solution = model.optimize()
    value = solution.objective_value
    ok = solution.status == "optimal" and value is not None and math.isfinite(value) and value >= 0
    return (float(value) if ok else None), solution.status


def fva(model, reactions, fraction=1.0):
    present = [r for r in reactions if r in model.reactions]
    Counter_.lps += 1 + 2 * len(present)
    frame = flux_variability_analysis(model, reaction_list=present, fraction_of_optimum=fraction, processes=1)
    return {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in present}


def pfba_vector(model):
    Counter_.lps += 2
    return pfba(model).fluxes.to_dict()


def per100(value, uptake=UPTAKE):
    return 100.0 * value / uptake


def close(model, closures):
    """closures: {reaction: 'fwd' | 'rev' | 'whole'}. Intersects the current bounds (never opens anything)."""
    for rid, how in closures.items():
        r = model.reactions.get_by_id(rid)
        lb, ub = r.bounds
        if how in ("fwd", "whole"):
            ub = min(ub, 0.0)
        if how in ("rev", "whole"):
            lb = max(lb, 0.0)
        if lb > ub:
            raise ValueError(f"closing {rid} {how} conflicts with a forced flux {r.bounds}")
        r.bounds = (lb, ub)


def ko(model):
    model.reactions.get_by_id(KO).bounds = (0.0, 0.0)


def mfa_rows():
    return {r["ID"]: r for r in csv.DictReader(PINS["w29_fluxes"][0].open())}


def r825_assumption(model):
    """Apply the R825 ASSUMPTION file (step 5 only), checking the reaction's current state against the file."""
    spec = json.loads(PINS["r825"][0].read_text())
    r = model.reactions.get_by_id(spec["reaction"])
    assert list(r.bounds) == spec["original_bounds"] and r.gene_reaction_rule == spec["gpr"]
    setattr(r, spec["field"], spec["value"])
    return spec["condition_id"]


def write_new(path, text):
    """Never overwrite an output."""
    with open(path, "x") as handle:
        handle.write(text)
