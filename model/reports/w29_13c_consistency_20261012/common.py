"""Shared setup for Task C (W29 13C consistency). Diagnostic only: nothing here writes a model.

Inputs are pinned by SHA-256. The base model, medium, caps and overlay come from commit e985dcd
(Task B, branch claude/task-b-e1ff88), extracted with `git show` into the research workspace.
Medium logic (V1/V2 W29-type defined, CSM = PO1f + CSM-Leu) and the caps helper follow Task B's
run_w29_defined.py / run_task_b.py.
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
from holistic_overlay import apply_overlay  # noqa: E402,F401

RESEARCH = Path(os.environ["IYALI26_RESEARCH_ROOT"]).resolve()
ART = RESEARCH / "artifacts/tasks" / HERE.name
IN = ART / "inputs"
PINS = {
    "model": (IN / "E5_coq9_resp_pkg_r132_20261012.xml", "a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515"),
    "medium": (IN / "po1f_csm_leu_w29uptake_20261012.csv", "6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45"),
    "caps": (IN / "w29_13c_upper_caps_20261012.json", "df1e33a8aa72c2d1b83129cb6a0e0d5e6eb283e3aec8feeecd844bcc2dc4257b"),
    "overlay": (HERE / "holistic_overlay.py", "bf6c7ce54d8338c365d98c82ca345ff2bd8c6055b38d402e8144a16b4abf30c4"),
    "profile": (RESEARCH / "state/strain_profiles/po1f_sd_leu.json", "35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383"),
    "experimental": (RESEARCH / "state/essentiality/repository/consensus_essential_genes.csv",
                     "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b"),
    "mtdna_scope": (REPO / "model/expected/mtdna_genes_out_of_scope.json", "16b3e9b4418f24821b8775936fd377ae21f77afa8e8af2b8949a246f62c44996"),
    "w29_fluxes": (RESEARCH / "reference/w29_worland2024/glucose_central_fluxes.csv",
                   "844615d0edfd8ffb033c1bce15ed746e9bb8ea02cbbe75016771b5c5fd08389f"),
}
PARAMS = dict(Threads=1, TimeLimit=60.0, Presolve=0, FeasibilityTol=1e-7, OptimalityTol=1e-7, IntFeasTol=1e-7)
UPTAKE = 2.414
MU = 0.2663
AMINO_ACID_ROWS = ("R1003", "R1202", "R1204", "R1215", "R1217", "R1220", "R1222", "R1223", "R1231", "R1232",
                   "R1233", "R1234")
PDH_FLOOR = 0.866629  # W29 glucose R19 (PDH) 95 % lower bound 86.6629 per 100 glucose


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_pins():
    for key, (path, digest) in PINS.items():
        assert sha(path) == digest, f"Input identity changed: {key}"
    return {k: {"path": str(p), "sha256": d} for k, (p, d) in PINS.items()}


class Counter_:
    lps = 0


def load(medium):
    """medium: 'V1', 'V2' (W29-type defined, uptake fixed 2.414) or 'CSM' (PO1f + CSM-Leu, uptake <= 2.414)."""
    sim = load_effective_simulation_context(model_path=PINS["model"][0], media_path=PINS["medium"][0],
                                            strain_profile_path=PINS["profile"][0])
    model = sim.model
    model.solver = "gurobi"
    configure_solver(model, {"solver": PARAMS})
    assert {r.id: c for r, c in linear_reaction_coefficients(model).items()} == {"biomass_C": 1.0}
    assert math.isclose(-model.reactions.R1070.lower_bound, UPTAKE)
    if medium in ("V1", "V2"):
        current = dict(model.medium)
        for rid in AMINO_ACID_ROWS:
            assert rid in current, rid
            current.pop(rid)
        if medium == "V2":
            current.pop("R1354")
        else:
            current["R1354"] = 1000.0
        model.medium = current
        r612 = model.reactions.R612
        assert tuple(r612.bounds) == (0.0, 0.0), "PO1f URA3 overlay not found"
        r612.bounds = (0.0, 1000.0)
        model.reactions.R1070.bounds = (-UPTAKE, -UPTAKE)
    return sim, model


def growth(model):
    Counter_.lps += 1
    solution = model.optimize()
    value = solution.objective_value
    ok = solution.status == "optimal" and value is not None and math.isfinite(value) and value >= 0
    return (float(value) if ok else None), solution.status


def fva(model, reactions, fraction=1.0):
    present = [r for r in reactions if r in model.reactions]
    Counter_.lps += 1 + 2 * len(present)
    try:
        frame = flux_variability_analysis(model, reaction_list=present, fraction_of_optimum=fraction, processes=1)
        return {rid: [float(frame.loc[rid, "minimum"]), float(frame.loc[rid, "maximum"])] for rid in present}
    except Exception as exc:  # recorded, not hidden
        return {"error": f"{type(exc).__name__}: {exc}"}


def pfba_vector(model):
    Counter_.lps += 2
    try:
        return pfba(model).fluxes.to_dict()
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def r825_reversible(model):
    """In-memory diagnostic only: R825 carnitine/acetylcarnitine antiport allowed in both directions."""
    r = model.reactions.R825
    assert tuple(r.bounds) == (0.0, 1000.0) and r.gene_reaction_rule == "YALI1C03415g"
    r.lower_bound = -1000.0


def pdh_floor(model, fraction=PDH_FLOOR):
    """PDH R693 >= fraction x actual glucose uptake (ASSUMPTION, diagnostic)."""
    expr = model.reactions.R693.flux_expression - fraction * (-model.reactions.R1070.flux_expression)
    model.add_cons_vars([model.problem.Constraint(expr, lb=0, name="pdh_floor")])


def add_caps(model, caps, gcv, substrate="glucose"):
    """Task B caps helper (PC, net malate import, one glycine-cleavage form); ASSUMPTION, diagnostic."""
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


def mfa_rows():
    return {r["ID"]: r for r in csv.DictReader(PINS["w29_fluxes"][0].open())}


def load_caps():
    return json.loads(PINS["caps"][0].read_text())
