"""Count how the model balances protons in cytosolic NAD/NADH and NADPH/NADP reactions (read-only, no solves).

Usage (from platform/): python ../model/reports/dolichol_route_20261028/scripts/proton_convention_scan.py MODEL.xml OUT.json
For every reaction with both members of a pair in C_cy and more than two species: the H+ m10[C_cy] coefficient
signed by the oxidised cofactor's coefficient, and the H / charge imbalance from check_mass_balance.
"""
import collections
import json
import sys
import warnings

from cobra.io import read_sbml_model

model_path, out = sys.argv[1:3]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(model_path)
PAIRS = {"NAD": ("m122[C_cy]", "m123[C_cy]"), "NADP": ("m41[C_cy]", "m40[C_cy]")}
result = {"species": {mid: [model.metabolites.get_by_id(mid).formula, model.metabolites.get_by_id(mid).charge]
                      for pair in PAIRS.values() for mid in pair}, "pairs": {}}
for label, (oxidised, reduced) in PAIRS.items():
    counts, examples = collections.Counter(), {}
    for r in model.reactions:
        s = {m.id: c for m, c in r.metabolites.items()}
        if oxidised in s and reduced in s and len(s) > 2:
            # H+ coefficient per oxidised cofactor formed (both read in the direction that forms it).
            h = s.get("m10[C_cy]", 0.0) * (-1 if s[oxidised] < 0 else 1) / abs(s[oxidised])
            bal = r.check_mass_balance()
            key = f"H+ per oxidised cofactor formed {h:+g}; imbalance " + (", ".join(f"{k} {v:+g}" for k, v in sorted(bal.items()) if k in ("H", "charge")) or "none")
            counts[key] += 1
            examples.setdefault(key, []).append(r.id)
    result["pairs"][label] = {k: {"count": v, "examples": examples[k][:5]} for k, v in counts.most_common()}
with open(out, "w") as handle:
    json.dump(result, handle, indent=1)
    handle.write("\n")
for label, rows in result["pairs"].items():
    for k, v in list(rows.items())[:4]:
        print(label, k, v["count"], v["examples"])
