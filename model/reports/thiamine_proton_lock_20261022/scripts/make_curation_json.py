"""Write model/curation/chemistry/thiamine_proton_lock.json from before.json plus the reviewed evidence text.

Usage (from the worktree root): python make_curation_json.py BEFORE.json OUT.json
"""
import json
import sys
from xml.sax.saxutils import unescape

before_path, out = sys.argv[1:3]
before = json.load(open(before_path))
assert before["sha256"] == "ef2265c1fd6c3cd12eebacc8c03901091f90289a6f77a7be66d2ede9d7379440"

PH73 = ("Checked on 2026-10-22 with every partner at its MetaNetX 4.5 pH 7.3 form, and recomputed by the "
        "thiamine phosphate chemistry audit: the reaction balances as written, H+ included. ")
EQUATIONS = {
    "R752": "thiamine(1+) + ATP(4-) = thiamine diphosphate(2-) + AMP(2-) + H+",
    "R408": "4-methyl-5-(2-hydroxyethyl)thiazole + ATP(4-) = HET-P(2-) + ADP(3-) + H+",
    "R750": "thiamine diphosphate(2-) + 2 H2O = thiamine(1+) + 2 hydrogenphosphate(2-) + H+",
}

reactions = {}
for rid, equation in EQUATIONS.items():
    row = before["targets"][rid]
    assert row["balance"] == {"charge": 1.0, "H": 1.0}, (rid, row["balance"])
    reactions[rid] = {
        "name": row["name"], "bounds": row["bounds"], "gpr": row["gpr"], "stoichiometry": row["stoichiometry"],
        # The residual the lock protects: exactly the pH 7.3 H+ under the interim KEGG forms.
        "residual": row["balance"],
        "notes_before": {k: unescape(v, {"&apos;": "'", "&quot;": '"'}) for k, v in row["notes"].items()},
        "evidence": PH73 + f"pH 7.3 equation: {equation}.",
        "evidence_tier": ("verified on 2026-10-22 by in-memory balance checks; the proton gate dry run on ef2265c1 "
                          "removes this H+ because the residual is exactly one H+ under the KEGG forms"),
        "rationale": ("The H+ is correct for the pH 7.3 convention the pipeline is moving to. Under the current "
                      "mixed forms, KEGG fully protonated thiamine phosphates and HET-P with neutral ATP, ADP, AMP "
                      "and phosphate, the reaction is off by exactly that H+; it balances only when all of these "
                      "move to pH 7.3 together. The lock keeps the proton/water gate, from the 2026-10-22 "
                      "microspecies.py on, from removing it, also on a model read back from SBML."),
        "previous_state": "no proton lock; the gate would remove the H+",
    }

spec = {
    "schema_version": 1,
    "curation_id": "THIAMINE-PROTON-LOCK-20261022",
    "status": "user_authorized_proton_water_lock",
    "authorization": ("2026-10-22 user: Approved: implement the H+ lock plan as written: gate honours a string note "
                      "lock_proton_water_stoichiometry = true, new opt-in --thiamine-proton-lock on top of ef2265c1 "
                      "for R752/R408/R750, default-build and control checks, no LPs. No commit or push."),
    "baseline_model_sha256": before["sha256"],
    "report": "model/reports/thiamine_proton_lock_20261022/",
    "lock_note": {"lock_proton_water_stoichiometry": "true"},
    "reactions": reactions,
}
# libSBML drops a note holding < or > and escapes ' and " on the next read, so new note texts avoid them.
for text in [spec["authorization"]] + [row[k] for row in reactions.values()
                                       for k in ("evidence", "evidence_tier", "rationale", "previous_state")]:
    assert not set(text) & set("<>&'\""), text
with open(out, "w") as handle:
    json.dump(spec, handle, indent=1)
    handle.write("\n")
print("wrote", out, sorted(reactions))
