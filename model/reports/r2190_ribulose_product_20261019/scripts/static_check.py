"""Static after-state checks of the new candidate: R2190 balance, species usage, build record (no LP)."""
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "platform"))
from cobra.io import read_sbml_model
from scripts.gem_annotate.coq9 import exact_residual

path = REPO / "model/candidates/E5_coq9_resp_pkg_r132_r2190_ribulose_20261019.xml"
model = read_sbml_model(str(path))
build = json.loads(path.with_suffix(".build.json").read_text())
r = model.reactions.R2190
usage = lambda mid: sorted(x.id for x in model.metabolites.get_by_id(mid).reactions)
record = {
    "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    "build_output_sha256": build["output"]["sha256"], "build_head": build["head"],
    "build_dirty_status": build["dirty_status"], "build_argv": build["argv"],
    "build_step_record": build["reaction_selection"].get("r2190_ribulose_product"),
    "requested_build_complete": build["requested_build_complete"],
    "source_unchanged_during_build": build["source_unchanged_during_build"],
    "R2190": {"equation": r.reaction, "named": r.build_reaction_string(use_metabolite_names=True),
              "bounds": list(r.bounds), "gpr": r.gene_reaction_rule, "exact_residual": exact_residual(r)},
    "m2061_reactions": usage("m2061[C_cy]"), "m454_reactions": usage("m454[C_cy]"),
    "m1823_reactions": usage("m1823[C_cy]"),
    "residuals_of_m454_reactions": {rid: exact_residual(model.reactions.get_by_id(rid)) for rid in usage("m454[C_cy]")},
}
out = Path(__file__).resolve().parents[1] / "static_check.json"
out.write_text(json.dumps(record, indent=1, default=str) + "\n")
print(json.dumps(record, indent=1, default=str))
