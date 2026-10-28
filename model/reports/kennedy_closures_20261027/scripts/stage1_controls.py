"""Stage 1 addendum: targeted controls for the two opened routes (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    stage1_controls.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json
Task A medium + PO1f, Gurobi 1e-9, checked solves. Variants closed / R1725 open / R1726 open / both open, each
with: WT; PSD double knockout (YALI1D27274g + YALI1D04432g); PSD double knockout with SPHPL, R730 and R663
(the sphingoid-base phosphate lyase reactions) also closed; PLD reactions R1743 + R1842 closed; the PLD gene
YALI1E22633g knocked out. Also lists every reaction named or annotated as diacylglycerol kinase (EC 2.7.1.107 /
2.7.1.174), to explain the DAG -> PA substitution of the R1726 cycle.
"""
import json
import sys
import warnings

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context

VARIANTS = {"closed": (), "R1725_open": ("R1725",), "R1726_open": ("R1726",), "both_open": ("R1725", "R1726")}
CONTROLS = {"wt": ((), ()),
            "psd_double_ko": (("YALI1D27274g", "YALI1D04432g"), ()),
            "psd_double_ko_lyase_closed": (("YALI1D27274g", "YALI1D04432g"), ("SPHPL", "R730", "R663")),
            "pld_reactions_closed": ((), ("R1743", "R1842")),
            "pld_gene_ko": (("YALI1E22633g",), ())}
model_path, medium, profile, out = sys.argv[1:5]


def context():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=model_path, media_path=medium, strain_profile_path=profile)
    configure(ctx.model)
    return ctx.model


result = {"variants": {}}
model = context()
result["dag_kinase_like"] = [
    {"reaction": r.id, "name": r.name, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule,
     "ec": r.annotation.get("ec-code"), "equation": r.build_reaction_string(True)[:200]}
    for r in model.reactions
    if "diacylglycerol kinase" in r.name.lower() or "dag kinase" in r.name.lower()
    or any(ec in str(r.annotation.get("ec-code", "")) for ec in ("2.7.1.107", "2.7.1.174"))]
for label, ids in VARIANTS.items():
    model = context()
    for rid in ids:
        model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)
    stoichiometry = matrix(model)
    rows = {}
    for name, (genes, reactions) in CONTROLS.items():
        with model:
            for gene in genes:
                model.genes.get_by_id(gene).knock_out()
            for rid in reactions:
                model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
            knocked = {r.id for g in genes for r in model.genes.get_by_id(g).reactions if r.bounds == (0, 0)} | set(reactions)
            rows[name] = solve(model, stoichiometry, knocked)
    result["variants"][label] = rows
    print(label, {k: (v["status"], v["growth"]) for k, v in rows.items()}, file=sys.stderr, flush=True)
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for row in result["dag_kinase_like"]:
    print(row)
