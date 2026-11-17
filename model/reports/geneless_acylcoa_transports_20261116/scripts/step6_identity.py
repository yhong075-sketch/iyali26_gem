"""Step 6, before the sweep: is 7ccbb568 the same optimisation problem as 353cefdc, the model of the baseline sweep?

No LP. Both models are loaded through load_effective_simulation_context (Task A + PO1f), and these are compared:
reaction ids and order, stoichiometry, effective bounds, gene rules, the objective, the gene ids and the runtime
exclusions. The baseline (sweep684_candidate of folinic_acid_module_20261112) is used only if all of them match.
Usage: python step6_identity.py  ->  tables/step6_identity.json
"""
import json
import sys
import warnings
from pathlib import Path

from common import MEDIUM, MODEL, MODEL_SHA, PROFILE, REPO, RR, check_inputs, sha256

BASE = RR / "artifacts/tasks/folinic_acid_module_20261112/candidate684_seed1.xml"
BASE_SHA = "353cefdc16ef1dc510eae49c168fc999254c2652c3994287af45589621a90cd9"
SWEEP = RR / "artifacts/tasks/folinic_acid_module_20261112/sweep684_candidate.tsv"
SWEEP_SHA = "49b4f27be05beffe32b206b75e7fe571b69d07af6d13a48b568f8e39f07e91db"

sys.path.insert(0, str(REPO / "platform"))
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context  # noqa: E402


def problem(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=path, media_path=MEDIUM, strain_profile_path=PROFILE)
    m = ctx.model
    return {
        "sha256": ctx.canonical_model_sha256,
        "reaction_order": [r.id for r in m.reactions],
        "reactions": {r.id: {"stoichiometry": sorted((k.id, v) for k, v in r.metabolites.items()),
                             "bounds": list(r.bounds), "rule": r.gene_reaction_rule} for r in m.reactions},
        "metabolite_order": [x.id for x in m.metabolites],
        "objective": str(m.objective.expression), "direction": m.objective.direction,
        "genes": sorted(g.id for g in m.genes),
        "excluded_runtime_genes": sorted(ctx.excluded_runtime_genes),
    }


check_inputs()
for path, want in [(BASE, BASE_SHA), (SWEEP, SWEEP_SHA)]:
    if sha256(path) != want:
        sys.exit(f"hash mismatch for {path}")
a, b = problem(MODEL), problem(BASE)
if a["sha256"] != MODEL_SHA or b["sha256"] != BASE_SHA:
    sys.exit("context hash mismatch")
diff = {k: a[k] == b[k] for k in a if k != "sha256"}
differing = sorted(rid for rid in a["reactions"] if a["reactions"][rid] != b["reactions"].get(rid))
out = {"model": MODEL_SHA, "baseline_model": BASE_SHA, "baseline_sweep_tsv": SWEEP_SHA, "same": diff,
       "reactions_differing": differing, "n_reactions": len(a["reactions"]), "n_genes": len(a["genes"]),
       "identical_problem": all(diff.values())}
dest = Path(__file__).resolve().parents[1] / "tables/step6_identity.json"
dest.write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out, indent=1))
