"""Mass/charge residual of every reaction the two steps touch, in the built candidate (read-only).

Usage: python -I balance_check.py CANDIDATE.xml DIFF.json LIPOYLATION.json ACYL_ACP.json OUT.json
Every reaction in DIFF.json is checked against the curated after-state: balanced, or exactly the residual
recorded in acyl_acp_chemistry.json (residuals_after) or lipoylation_module.json
(preexisting_residuals_after_fill, replaced by the acyl-ACP record where that step corrected the reaction).
biomass_C is a pseudo-reaction and is reported, not judged.
"""
import json
import sys
import warnings

from cobra.io import read_sbml_model

model_path, diff_path, lipo_path, acyl_path, out = sys.argv[1:6]
lipo, acyl = json.load(open(lipo_path)), json.load(open(acyl_path))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(model_path)
touched = sorted({c["item"] for c in json.load(open(diff_path)) if c["section"] == "reactions"})
expected = {rid: {} for rid in lipo["balanced_after"]}
expected.update(lipo["preexisting_residuals_after_fill"])
expected.update({rid: r for rid, r in acyl["residuals_after"].items()})


def residual(r):
    return {k: v for k, v in r.check_mass_balance().items() if abs(v) > 1e-9}


def same(a, b):
    return b is not None and set(a) == set(b) and all(abs(a[k] - b[k]) <= 1e-9 for k in b)


rows = {}
for rid in touched:
    if rid == "biomass_C":  # contains formula-less species; not a balanceable reaction
        rows[rid] = {"residual": None, "expected": None, "verdict": "pseudo-reaction"}
        continue
    reaction = model.reactions.get_by_id(rid)
    if len(reaction.metabolites) == 1 and rid in lipo["new_reactions"]:  # R2290, the labelled 5'-dA sink
        rows[rid] = {"residual": residual(reaction), "expected": None, "verdict": "one-sided sink by design"}
        continue
    actual = residual(reaction)
    rows[rid] = {"residual": actual, "expected": expected.get(rid),
                 "verdict": "pseudo-reaction" if rid == "biomass_C" else
                 ("balanced" if not actual and same(actual, expected.get(rid, {})) else
                  "curated residual" if same(actual, expected.get(rid)) else "MISMATCH")}
with open(out, "w") as handle:
    json.dump(rows, handle, indent=1)
    handle.write("\n")
for rid, row in rows.items():
    print(rid, row["verdict"], row["residual"] if row["verdict"] != "pseudo-reaction" else "")
print("mismatches:", [rid for rid, row in rows.items() if row["verdict"] == "MISMATCH"])
