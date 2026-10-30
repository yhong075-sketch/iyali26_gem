"""Colour scales and scenario titles shared by render.py (gap-census palette, adapted from the UM11 figure).

Scenarios: pFBA at WT (identical in both models) and pFBA with the alpha subunit YALI1F04433g knocked out, before
(d8d1bac7) and after (ca43ae0c) the R171 rule change.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12}]
# pFBA |v| at the optimum; the energy fluxes are ~1 to 20, the proton leak ~150
PFBA_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#a7d7a0", "size": 12},
              {"type": "value", "value": 1, "color": "#4fae5b", "size": 18},
              {"type": "value", "value": 10, "color": "#1b6b2a", "size": 26},
              {"type": "value", "value": 100, "color": "#0b3d16", "size": 32}]
RXN_SCALES = {"pfba_wt": PFBA_SCALE, "pfba_ko_before": PFBA_SCALE, "pfba_ko_after": PFBA_SCALE}
TITLES = {
    "pfba_wt": "WT, pFBA, Task A + PO1f: growth 0.3117 in both models (the rule change does not touch the WT optimum). ATP synthase R171 carries 17.1",
    "pfba_ko_before": "BEFORE (d8d1bac7): alpha subunit YALI1F04433g knocked out. R171 stays open through the single-gene OR branches g / k: growth 0.3117 = WT, fluxes = WT",
    "pfba_ko_after": "AFTER (ca43ae0c): alpha subunit YALI1F04433g knocked out. R171 closed by the 12-gene AND: growth 0.0342 (KO/WT 0.1097); proton leak R1384 dissipates the gradient",
}
