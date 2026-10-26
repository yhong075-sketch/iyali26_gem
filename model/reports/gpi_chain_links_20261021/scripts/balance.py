"""Mass and charge balance of the five added reactions in a built model (read-only).

Usage: python balance.py MODEL.xml
"""
import sys
import warnings

from cobra.io import read_sbml_model

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = read_sbml_model(sys.argv[1])
for rid in ["R_GPIGNT", "R_PEt_em_er", "R_DAGt_er_em", "R_ACt_er_cy", "R_DM_GPI"]:
    reaction = model.reactions.get_by_id(rid)
    balance = "boundary" if reaction.boundary else reaction.check_mass_balance()
    print(f"{rid:13s} {reaction.bounds} {reaction.build_reaction_string()} | balance {balance}")
