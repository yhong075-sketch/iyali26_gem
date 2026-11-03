"""Fluxes for the five-family status figure (read-only): e1e7c868 as built and each family's recommended variant as an
in-memory what-if (whatif.py, the scouts' builders), each applied alone. None was implemented.

Usage: IYALI26_RESEARCH_ROOT=... python flux.py OUT.json
Model and condition are the archived copies the scouts used (e1e7c868.xml, cond_taskA_po1f.json: Task A medium,
R612 closed; no plasmid LEU2 complement of R45), loaded exactly as the scouts loaded them. Gurobi, model.tolerance 1e-9.
Stops unless both SHA-256 match, WT of every variant reproduces the archived value within 1e-9, and the LP count stays
within 1,500. Per variant: WT, pFBA, FVA at 0 % of the optimum over that family's mapped reactions (capacity), and the
energy check (every exchange (0, 1000), xMAINTENANCE lb 0, maximise R_NTP1). As built: the same over all sections.
"""
import hashlib
import json
import logging
import os
import sys
import warnings
from pathlib import Path

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba

sys.path.insert(0, str(Path(__file__).parent))
from whatif import ARCHIVE, VARIANTS  # noqa: E402

warnings.simplefilter("ignore")
logging.getLogger("cobra").setLevel(logging.ERROR)
MODEL, COND = ARCHIVE / "e1e7c868.xml", ARCHIVE / "cond_taskA_po1f.json"
MODEL_SHA = "e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897"
COND_SHA = "db78d70531a447535147c7e462d7d9add9589ec38345d15ec7f894ee29023a63"
BUDGET = 1500
SECTIONS = {
    "sphingolipid": ["R722", "R37", "R1406", "R1407", "R1408", "R1409", "R1410", "R1411", "R1412", "R1413", "R1414",
                     "R1415", "R1416", "R1417", "R1418", "R1419", "R1420", "R1421", "R1422", "R1423", "R1424", "R1425",
                     "R1426", "R1427", "R1428", "R1429", "R1430", "R1431", "R1432", "R1433", "R1508", "R1509", "R1510",
                     "R664", "R198", "R200", "R2136", "R199", "R201", "R206", "R207", "R208", "R209",
                     "R1815", "R1816", "R1817", "R1818", "R1819", "R1820", "R1821", "R1822", "R1823", "R1824",
                     "R1825", "R1826", "R1827", "R1828", "R1829", "R1830", "R1831", "R1832", "R1833", "R1834",
                     "R2143", "R2148", "R2151", "R2154", "R2157", "R2161"],
    "gpi": ["R_GPIGNT", "R371", "R113", "R372", "R373", "R374", "R375", "R376", "R_DM_GPI"],
    "cardiolipin": ["R1791", "R197", "R1836", "R1722", "R1721", "R1723", "R1724"],
    "nglycan": ["R2108", "R2107", "R2106", "R2105", "R173", "R141", "R289", "R850", "R1137", "R602", "R586", "R2181",
                "R2177", "R287"],
    "mito_translation": ["R155", "R158", "R166", "R338", "R393", "R494", "R523", "R532", "R548", "R626", "R760", "R772",
                         "R776", "R799", "R546", "R1251", "R1218", "R122", "R123"],
}
FAMILY_OF = {"sph_V2": "sphingolipid", "gpi_b1": "gpi", "cl_a": "cardiolipin", "ng_C2": "nglycan", "ng_B": "nglycan",
             "mt_rec": "mito_translation"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


for path, want in ((MODEL, MODEL_SHA), (COND, COND_SHA)):
    if sha(path) != want:
        sys.exit(f"STOP: {path} SHA-256 differs from {want}")
J = lambda rel: json.load(open(ARCHIVE / rel))
SPH, NGV = J("sphingolipid/verify_lp/verify2.json"), J("nglycan/verify_lp/verify_c2.json")
EXPECTED = {  # archived WT values (LP lens where it ran the variant; scout for b1)
    "as_built": SPH["baseline"]["wt_growth"],
    "sph_V2": SPH["v2"]["wt_growth"],
    "gpi_b1": J("gpi/s2_variants.json")["variants"]["b1"]["wt"]["growth"],
    "cl_a": J("cardiolipin/verify_lp/verify_cardiolipin.json")["variant_a"]["wt"]["growth"],
    "ng_C2": NGV["C2"]["wt_growth"],
    "ng_B": NGV["B"]["wt_growth"],
    "mt_rec": J("mito_translation/verify_lp/verify2.json")["recommended"]["wt"]["growth"],
}
out = sys.argv[1]
base = cobra.io.read_sbml_model(str(MODEL))
base.solver = "gurobi"
base.tolerance = 1e-9
cond = json.loads(COND.read_text())
base.medium = cond["medium"]
for rid, (lb, ub) in cond["bounds"].items():
    base.reactions.get_by_id(rid).bounds = (lb, ub)
LP = {"n": 0}


def count(n):
    LP["n"] += n
    if LP["n"] > BUDGET:
        sys.exit(f"STOP: LP budget {BUDGET} exceeded ({LP['n']})")


pair = lambda df, rid: [float(df.loc[rid, "minimum"]), float(df.loc[rid, "maximum"])]
base_rx = {r.id: (tuple(r.bounds), {m.id: c for m, c in r.metabolites.items()}) for r in base.reactions}
result = {"model": str(MODEL), "model_sha256": MODEL_SHA, "condition": str(COND), "condition_sha256": COND_SHA,
          "condition_note": "Task A medium (cond_taskA_po1f.json) and R612 (0, 0); no plasmid LEU2 complement of R45",
          "solver": f"gurobi {'.'.join(map(str, __import__('gurobipy').gurobi.version()))}, cobra {cobra.__version__}, "
                    "model.tolerance 1e-9", "sections": SECTIONS, "family_of": FAMILY_OF, "expected_wt": EXPECTED,
          "variants": {}}
reactions, metabolites = {}, {}
for label in ["as_built"] + list(VARIANTS):
    m = base.copy()
    applied = VARIANTS[label](m) if label in VARIANTS else {"biomass_term": {}}
    added = sorted(r.id for r in m.reactions if r.id not in base_rx)
    changed = sorted(r.id for r in m.reactions if r.id in base_rx and
                     (tuple(r.bounds), {x.id: c for x, c in r.metabolites.items()}) != base_rx[r.id] and r.id != "biomass_C")
    fam = FAMILY_OF.get(label)
    here = [rid for sec, ids in SECTIONS.items() if fam in (None, sec) for rid in ids] + added
    sol = m.optimize()
    count(1)
    wt = sol.objective_value if sol.status == "optimal" else None
    if wt is None or abs(wt - EXPECTED[label]) > 1e-9:
        sys.exit(f"STOP: {label} WT {sol.status} {wt} differs from the archived {EXPECTED[label]}")
    p = pfba(m)
    count(2)
    cap = flux_variability_analysis(m, reaction_list=here, fraction_of_optimum=0.0, processes=1)
    count(1 + 2 * len(here))
    with m:
        for r in m.exchanges:
            r.bounds = (0.0, 1000.0)
        m.reactions.xMAINTENANCE.lower_bound = 0.0
        m.objective = "R_NTP1"
        m.objective_direction = "max"
        e = m.optimize()
        count(1)
    result["variants"][label] = {
        "family": fam, "applied": applied, "added": added, "changed": changed, "wt": wt, "wt_status": sol.status,
        "present": here, "bounds": {rid: list(m.reactions.get_by_id(rid).bounds) for rid in here},
        "pfba": {"status": p.status, "fluxes": {rid: float(p.fluxes[rid]) for rid in here}},
        "capacity_taskA": {rid: pair(cap, rid) for rid in here},
        "egc_R_NTP1_no_uptake": {"status": e.status, "value": e.objective_value if e.status == "optimal" else None}}
    for r in (m.reactions.get_by_id(rid) for rid in here):
        key = r.id if label == "as_built" or r.id in added else None
        if key and key not in reactions:
            reactions[key] = {"name": r.name, "lb": r.lower_bound, "ub": r.upper_bound, "gpr": r.gene_reaction_rule,
                              "genes": sorted(g.id for g in r.genes), "subsystem": r.subsystem,
                              "mets": {x.id: c for x, c in r.metabolites.items()}, "origin": label}
        for x in r.metabolites:
            metabolites.setdefault(x.id, {"name": x.name, "comp": x.compartment, "degree": {}, "reactions": {}})
            metabolites[x.id]["degree"][label] = len(x.reactions)
            metabolites[x.id]["reactions"][label] = sorted(y.id for y in x.reactions)
    print(f"{label}: WT {wt!r}, mapped {len(here)}, added {len(added)}, changed {changed}, "
          f"EGC {e.status} {result['variants'][label]['egc_R_NTP1_no_uptake']['value']}, LPs so far {LP['n']}", flush=True)
result["reactions"] = reactions
result["metabolites"] = metabolites
result["lp_count"] = LP["n"]
result["lp_count_rule"] = "per variant: WT 1 + pFBA 2 + FVA (1 + 2 per mapped reaction) + energy check 1"
json.dump(result, open(out, "w"), indent=1)
print("LPs", LP["n"], "->", out)
