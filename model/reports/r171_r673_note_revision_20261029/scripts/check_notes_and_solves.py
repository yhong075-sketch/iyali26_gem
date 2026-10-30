"""Check the notes-only revision: changed note keys and values, then 3 sanity solves (read-only).

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ../model/reports/r171_r673_note_revision_20261029/scripts/check_notes_and_solves.py OLD.xml NEW.xml MEDIUM PROFILE OUT.json
Compares every note key of R171 and R673 between the two models and checks the new values against the curated data;
then on NEW: checked WT, and checked knockouts of YALI1F04433g and YALI1E05929g (Gurobi 1e-9, reset per solve).
"""
import json
import sys
import warnings
from html import unescape
from pathlib import Path

from cobra.io import read_sbml_model

from checked_solve import configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL

old_path, new_path, medium, profile, out = sys.argv[1:6]
old, new = read_sbml_model(Path(old_path)), read_sbml_model(Path(new_path))
curated = {"R171": json.loads(MODEL.curation_file("r171_atp_synthase_gpr.json").read_text())["reaction"]["notes"],
           "R673": json.loads(MODEL.curation_file("r673_prolyl_split.json").read_text())["reaction"]["notes"]}
result = {"notes": {}}
for rid, notes in curated.items():
    a = {k: unescape(v) for k, v in old.reactions.get_by_id(rid).notes.items()}
    b = {k: unescape(v) for k, v in new.reactions.get_by_id(rid).notes.items()}
    result["notes"][rid] = {"changed_keys": sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k)),
                            "curated_keys_equal_data": all(b.get(k) == v for k, v in notes.items())}
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    ctx = load_effective_simulation_context(model_path=new_path, media_path=medium, strain_profile_path=profile)
model = ctx.model
configure(model)
wt = solve(model, matrix(model), ())
result["solves"] = {"wt": wt}
for gene in ("YALI1F04433g", "YALI1E05929g"):
    with model:
        model.genes.get_by_id(gene).knock_out()
        knocked = tuple(r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0))
        ko = solve(model, matrix(model), knocked)
        result["solves"][gene] = {**ko, "closed": list(knocked), "ratio": ko["growth"] / wt["growth"] if ko["passed"] else None}
result["model_sha256"] = ctx.canonical_model_sha256
result["condition"] = {"medium": medium, "medium_sha256": ctx.medium_sha256, "strain_profile": profile,
                       "strain_profile_sha256": ctx.strain_profile_sha256,
                       "solver": "Gurobi via checked_solve.configure (FeasibilityTol / OptimalityTol 1e-9, reset per solve)"}
result["lp_count"] = 3
Path(out).write_text(json.dumps(result, indent=1) + "\n")
print(json.dumps({"notes": result["notes"], "wt": wt["growth"],
                  **{g: (v["growth"], v["ratio"], v["closed"], v["passed"]) for g, v in result["solves"].items() if g != "wt"}}, indent=1))
