"""The authorized LPs for the lipoylation candidate (GLPK, SD-Leu + PO1f).

Usage: python solves_lipoylation.py MEDIUM.csv PROFILE.json BASELINE.xml CANDIDATE.xml OUT.json
25 LPs (one more was spent on the in-memory prototype, see REPORT.md):
- WT growth on the baseline and the candidate (2).
- Growth not required, candidate: max R1405 (mtFAS octanoyl-ACP step) (1).
- Growth not required, candidate: max R2218, min and max R2216 (the only entries into the R2216-R2241
  export block) (3).
- Single-gene KO on the candidate for 14 genes (14), and on the baseline for the 3 genes that already
  had reactions there (YALI1D30618g, YALI1E20192g, YALI1E00574g) (3).
- Alpha sensitivity on the candidate: biomass lipoyl terms scaled to alpha 1e-6 and 1e-4 (2).
A non-optimal status is recorded as such, not as growth 0.
"""
import json
import os
import sys

os.environ.setdefault("HOME", os.path.abspath("cache"))
import pandas as pd
from cobra.io import read_sbml_model

KO_CANDIDATE = ["YALI1F34729g", "YALI1E13472g", "YALI1C04991g", "YALI1E00574g", "YALI1D30618g", "YALI1E20192g",
                "YALI1D18037g", "YALI1D32594g", "YALI1A20089g", "YALI1C27933g", "YALI1E22262g", "YALI1F38317g",
                "YALI1C26939g", "YALI1F37498g"]
KO_BASELINE = ["YALI1D30618g", "YALI1E20192g", "YALI1E00574g"]
LIPOYL_TERMS = {"m586[C_mi]": -1.0, "m2062[C_mi]": 1.0, "m2065[C_mi]": -1.0, "m2064[C_mi]": 1.0}


def load(path, medium_csv, profile_json):
    model = read_sbml_model(path)
    model.solver = "glpk"
    model.solver.configuration.presolve = True
    med = pd.read_csv(medium_csv, dtype={"exchange": str})
    medium = {r.exchange.strip(): float(r.uptake) for r in med.itertuples() if float(r.uptake) > 0}
    prof = json.load(open(profile_json))
    medium.update({k: float(v) for k, v in prof.get("medium", {}).get("runtime_uptake_overrides", {}).items()})
    model.medium = {k: v for k, v in medium.items() if k in model.reactions}
    for op in prof.get("operations", []):
        if op.get("type") == "disable_reaction":
            model.reactions.get_by_id(op["reaction_id"]).bounds = (op["set_bounds"]["lower_bound"],
                                                                   op["set_bounds"]["upper_bound"])
    return model


def solve(model):
    s = model.optimize()
    return {"status": s.status, "value": s.objective_value if s.status == "optimal" else None}


def ko(model, gid, wt):
    with model:
        gene = model.genes.get_by_id(gid)
        gene.knock_out()
        off = sorted(r.id for r in gene.reactions if r.bounds == (0, 0))
        s = solve(model)
    ratio = s["value"] / wt if (s["value"] is not None and wt) else None
    return {**s, "ko_wt_ratio": ratio, "reactions_off": off}


if __name__ == "__main__":
    medium_csv, profile_json, base_path, cand_path, out = sys.argv[1:6]
    res, n = {}, 0
    base = load(base_path, medium_csv, profile_json)
    cand = load(cand_path, medium_csv, profile_json)
    res["wt_baseline"] = solve(base); n += 1
    res["wt_candidate"] = solve(cand); n += 1
    wt_b, wt_c = res["wt_baseline"]["value"], res["wt_candidate"]["value"]
    probes = {}
    for rid, sense in (("R1405", "max"), ("R2218", "max"), ("R2216", "min"), ("R2216", "max")):
        with cand:
            cand.objective = rid
            cand.objective_direction = sense
            probes[f"{rid}_{sense}"] = solve(cand); n += 1
    res["growth_not_required_probes"] = probes
    res["ko_candidate"] = {}
    for gid in KO_CANDIDATE:
        res["ko_candidate"][gid] = ko(cand, gid, wt_c); n += 1
    res["ko_baseline"] = {}
    for gid in KO_BASELINE:
        res["ko_baseline"][gid] = ko(base, gid, wt_b); n += 1
    res["alpha_sensitivity"] = {"1e-05": wt_c}
    biomass = cand.reactions.biomass_C
    for alpha in (1e-6, 1e-4):
        with cand:
            biomass.add_metabolites({cand.metabolites.get_by_id(m): s * alpha - biomass.metabolites[cand.metabolites.get_by_id(m)]
                                     for m, s in LIPOYL_TERMS.items()}, combine=True)
            res["alpha_sensitivity"][f"{alpha:.0e}"] = solve(cand)["value"]; n += 1
    res["lp_count"] = n
    res["inputs"] = {"baseline": base_path, "candidate": cand_path, "medium": medium_csv, "profile": profile_json,
                     "solver": "glpk (presolve)"}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))
