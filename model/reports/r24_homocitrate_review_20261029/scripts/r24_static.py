"""Static facts about R24/R395 and their metabolites in two models (read-only): stoichiometry, balance,
charges, and every reaction touching the pathway intermediates.

Usage: /opt/miniconda3/bin/python r24_static.py SOURCE.xml CANDIDATE.xml
"""
import sys
import warnings
from pathlib import Path

import cobra

warnings.simplefilter("ignore")
for path in sys.argv[1:]:
    m = cobra.io.read_sbml_model(Path(path))
    print("#####", Path(path).name)
    for rid in ("R1170", "R24", "R395", "R398", "R396"):
        r = m.reactions.get_by_id(rid)
        print(f"{rid} {r.bounds} {r.build_reaction_string()} | rule {r.gene_reaction_rule}"
              f" | imbalance {r.check_mass_balance()}")
    for mid in ("m75[C_mi]", "m76[C_mi]", "m624[C_mi]", "m628[C_nu]"):
        x = m.metabolites.get_by_id(mid)
        print(f"  {mid} {x.name} {x.formula} charge {x.charge} reactions {sorted(r.id for r in x.reactions)}")
