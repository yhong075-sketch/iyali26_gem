"""The authorized LPs for the acyl-ACP chemistry candidate (GLPK, SD-Leu + PO1f).

Usage: python solves_acyl_acp.py MEDIUM.csv PROFILE.json BASELINE.xml CANDIDATE.xml OUT.json
12 LPs:
- WT growth on the baseline (lipoylation v2) and the candidate (2).
- Growth not required, candidate: max R2230 (C14 -> C16 elongation), min R2225, max R2231, max R2216,
  max R1405 (mtFAS octanoyl-ACP step) (5). The transports are written cytosol -> mitochondrion, so a
  positive flux is import: min R2225 probes myristoyl-ACP export, max R2231 / max R2216 probe palmitoyl-
  / octanoyl-ACP import. (Labels corrected after the audit; in the baseline cycle R2225 runs +36.39
  (import) and R2231 -36.39 (export), so these probes do not test the cycle; the audit found min R2225 =
  3.5e-15 on the baseline, max R2231 / R2216 were not run there. Closure rests on max R2230 and the static proof in REPORT.md section 5.)
- Growth not required, baseline: max R2230, for contrast with the audit's 36.389 (1).
- Single-gene KO on the candidate: YALI1C27933g, YALI1D32594g, YALI1B19844g (FAS1), YALI1B25427g (FAS2) (4).
A non-optimal status is recorded as such, not as growth 0.
"""
import json
import os
import sys

os.environ.setdefault("HOME", os.path.abspath("cache"))
import pandas as pd
from cobra.io import read_sbml_model

KO = ["YALI1C27933g", "YALI1D32594g", "YALI1B19844g", "YALI1B25427g"]


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


def probe(model, rid, sense):
    with model:
        model.objective = rid
        model.objective_direction = sense
        return solve(model)


if __name__ == "__main__":
    medium_csv, profile_json, base_path, cand_path, out = sys.argv[1:6]
    base, cand = load(base_path, medium_csv, profile_json), load(cand_path, medium_csv, profile_json)
    res, n = {}, 0
    res["wt_baseline"] = solve(base); n += 1
    res["wt_candidate"] = solve(cand); n += 1
    res["candidate_probes"] = {}
    for rid, sense in (("R2230", "max"), ("R2225", "min"), ("R2231", "max"), ("R2216", "max"), ("R1405", "max")):
        res["candidate_probes"][f"{rid}_{sense}"] = probe(cand, rid, sense); n += 1
    res["baseline_probes"] = {"R2230_max": probe(base, "R2230", "max")}; n += 1
    wt = res["wt_candidate"]["value"]
    res["ko_candidate"] = {}
    for gid in KO:
        with cand:
            gene = cand.genes.get_by_id(gid)
            gene.knock_out()
            s = solve(cand)
        res["ko_candidate"][gid] = {**s, "ko_wt_ratio": (s["value"] / wt) if (s["value"] is not None and wt) else None,
                                    "reactions": sorted(r.id for r in gene.reactions)}
        n += 1
    res["lp_count"] = n
    res["inputs"] = {"baseline": base_path, "candidate": cand_path, "medium": medium_csv, "profile": profile_json,
                     "solver": "glpk (presolve)"}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res, indent=1))
