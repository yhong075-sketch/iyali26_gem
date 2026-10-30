"""Reproduce the stale knock_out of a reactionless gene in cobra (read-only).

Usage: python -I repro_functional.py MODEL.xml   (MODEL = e1e7c868, where YALI1E17935g has no reaction)
"""
import sys
from pathlib import Path

from cobra.io import read_sbml_model

A1, A2 = "YALI1D11984g", "YALI1E17935g"
m = read_sbml_model(Path(sys.argv[1]))
print("fresh read: ACO2 linked to model:", m.genes.get_by_id(A2)._model is m,
      "| ACO1 linked to model:", m.genes.get_by_id(A1)._model is m)


def ko(genes):
    with m:
        for g in genes:
            m.genes.get_by_id(g).knock_out()
        closed = sorted({r.id for g in genes for r in m.genes.get_by_id(g).reactions if r.bounds == (0, 0)})
    return closed


print("fresh ACO2 functional", m.genes.get_by_id(A2).functional)
print("ACO2 KO (no reactions):", ko((A2,)), "-> ACO2 functional after context exit:", m.genes.get_by_id(A2).functional)
with m:
    m.reactions.R24.gene_reaction_rule = f"{A1} or {A2}"
    print("same model, rule C; ACO1 KO closes:", ko((A1,)))
m2 = read_sbml_model(Path(sys.argv[1]))
with m2:
    m2.reactions.R24.gene_reaction_rule = f"{A1} or {A2}"
    with m2:
        m2.genes.get_by_id(A1).knock_out()
        print("fresh model, rule C; ACO1 KO R24 bounds:", m2.reactions.R24.bounds)
