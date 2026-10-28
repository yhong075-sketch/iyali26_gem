"""Mass and charge balance of the two touched reactions, before and after (read-only).

Usage: python balance.py BASELINE.xml NEW.xml
biomass_C is a pseudo-reaction; R1342's residual comes from the formula-less TPP and the mis-annotated
TMP (C6H4Cl3NO). The point is that neither residual changes except by the TPP term itself.
"""
import sys
import warnings

from cobra.io import read_sbml_model

for path in sys.argv[1:3]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = read_sbml_model(path)
    print(path.split("/")[-1])
    for rid in ("biomass_C", "R1342"):
        reaction = model.reactions.get_by_id(rid)
        try:
            balance = reaction.check_mass_balance()
        except ValueError as error:  # e.g. chitin H2O(C8H13NO5)n in biomass_C
            balance = f"not computable: {error}"
        print(f"  {rid:9s} {reaction.bounds} terms={len(reaction.metabolites)} balance {balance}")
