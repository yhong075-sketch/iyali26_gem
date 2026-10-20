"""Shared helpers for the orphan-curation scripts (run in the cobra env as subprocesses)."""
import os, json
os.environ.setdefault("HOME", os.path.abspath("cache"))
os.makedirs(os.environ["HOME"], exist_ok=True)


def load_model(path):
    import cobra
    M = cobra.io.read_sbml_model(path)
    M.solver = "glpk"
    M.solver.configuration.presolve = True
    return M


def read_condition(medium_csv, profile_json=None):
    import pandas as pd
    med = pd.read_csv(medium_csv, dtype={"exchange": str})
    medium = {r.exchange.strip(): float(r.uptake) for r in med.itertuples() if float(r.uptake) > 0}
    prof = json.load(open(profile_json)) if profile_json else {}
    medium.update({k: float(v) for k, v in prof.get("medium", {}).get("runtime_uptake_overrides", {}).items()})
    closed = {op["reaction_id"]: (op["set_bounds"]["lower_bound"], op["set_bounds"]["upper_bound"])
              for op in prof.get("operations", []) if op.get("type") == "disable_reaction"}
    return medium, closed


def setup(M, medium, closed):
    M.medium = {k: v for k, v in medium.items() if k in M.reactions}
    for rid, b in closed.items():
        if rid in M.reactions:
            M.reactions.get_by_id(rid).bounds = b


def growth(M):
    v = M.slim_optimize()
    return 0.0 if (v != v or v < 1e-9) else v
