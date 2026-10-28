"""Diagnose why FAS / OLE1 / acyl-CoA ligase knockouts are viable once R1726 is open (read-only, about 60 LPs).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/kennedy_closures_20261027/scripts/acyl_leak_diagnosis.py CANDIDATE.xml MEDIUM.csv PROFILE.json OUT.json
Task A medium + PO1f, Gurobi 1e-9, checked solves. For each of the four genes that changed class
(YALI1B19844g, YALI1B25427g, YALI1C07638g, YALI1D22124g): the knockout on the candidate, then the same knockout
with one suspect reaction closed at a time (the R1726 cycle members, the PDAT / PE-DAG acyltransferases, the
phospholipases, the 0.01-scaled PA phosphatase R2145, the lipid-particle lipases). A suspect whose closure
makes the knockout lethal again lies on the rescue route. Also records the pFBA lipid fluxes of the FAS
knockout and the mass balance of every reaction that carries flux there and touches a lipid species.
"""
import json
import re
import sys
import warnings

from cobra.flux_analysis import pfba

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

GENES = ("YALI1B19844g", "YALI1B25427g", "YALI1C07638g", "YALI1D22124g")
SUSPECTS = ("R1726", "R1842", "R1743", "R1806", "R1729", "R1728", "R1730", "R1731", "R1732", "R1733", "R2145",
            "R1763", "R1765", "R1771", "R354", "R1776", "R1778")
LIPID = re.compile(r"phosphatid|glyceride|diacylglycerol|acyl|POOL|choline|lipid|fatty|-CoA|oleate|palmit|stear", re.I)
model_path, medium, profile, out = sys.argv[1:5]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
configure(model)
stoichiometry = matrix(model)
result = {"wt": solve(model, stoichiometry, set()), "knockouts": {}}
for gene in GENES:
    rows = {}
    with model:
        model.genes.get_by_id(gene).knock_out()
        knocked = {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}
        rows["ko_only"] = solve(model, stoichiometry, knocked)
        for rid in SUSPECTS:
            with model:
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
                rows[f"close_{rid}"] = solve(model, stoichiometry, knocked | {rid})
        if gene == "YALI1B19844g":
            fluxes = pfba(model).fluxes
            lipid = []
            for r in model.reactions:
                v = float(fluxes[r.id])
                if abs(v) > 1e-9 and LIPID.search(r.name + " " + " ".join(m.name for m in r.metabolites)):
                    try:
                        imbalance = {k: float(x) for k, x in r.check_mass_balance().items()}
                    except ValueError as error:  # a species without a formula
                        imbalance = {"unknown": str(error)}
                    lipid.append({"reaction": r.id, "name": r.name[:50], "flux": v, "gpr": r.gene_reaction_rule[:60],
                                  "imbalance": imbalance})
            result["fas_ko_pfba_lipid_fluxes"] = sorted(lipid, key=lambda row: -abs(row["flux"]))
    result["knockouts"][gene] = rows
    print(gene, {k: (None if v["growth"] is None else round(v["growth"], 7)) for k, v in rows.items()}, file=sys.stderr, flush=True)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
