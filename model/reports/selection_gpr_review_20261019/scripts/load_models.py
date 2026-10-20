"""Load source, metadata.xml and candidates; dump reaction fields + notes to JSON (static, no LP)."""
import json, sys, cobra, logging
logging.getLogger("cobra").setLevel(logging.ERROR)
W = "/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/orphan-reaction-newest-model-b8ddc8"
R = "/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/artifacts/tasks"
paths = {
    "source": f"{W}/model/source/iyali26.xml",
    "metadata": f"{R}/model_metadata_trna_push_20260910/staged-checkout/artifacts/coq9_pipeline_integration_20260909/accepted/metadata.xml",
    "a56f": f"{W}/model/candidates/E5_coq9_resp_pkg_r132_20261012.xml",
    "c989": f"{W}/model/candidates/E5_coq9_resp_pkg_r132_gpr_restore_v2_20261019.xml",
    "trna": f"{W}/model/candidates/model_metadata_trna.xml",
}
out = {}
for k, p in paths.items():
    m = cobra.io.read_sbml_model(p)
    out[k] = {
        "reactions": {r.id: {"stoich": {mm.id: c for mm, c in r.metabolites.items()}, "bounds": list(r.bounds),
                             "gpr": r.gene_reaction_rule, "name": r.name, "notes": {a: str(b) for a, b in r.notes.items()},
                             "ec": r.annotation.get("ec-code")} for r in m.reactions},
        "mets": {mm.id: {"formula": mm.formula, "charge": mm.charge, "compartment": mm.compartment, "name": mm.name} for mm in m.metabolites},
        "genes": {g.id: {"name": g.name, "annotation": {a: (b if isinstance(b, (str, list)) else str(b)) for a, b in g.annotation.items()}} for g in m.genes},
    }
    print(k, len(m.reactions), len(m.metabolites), len(m.genes), flush=True)
json.dump(out, open(sys.argv[1], "w"))
