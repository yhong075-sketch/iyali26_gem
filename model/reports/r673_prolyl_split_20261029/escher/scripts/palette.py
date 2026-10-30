"""Colour scales and scenario titles shared by render.py (gap-census palette, adapted from the UM11 figure).

Scenarios: pFBA at WT (identical in both models) and pFBA with YALI1E05929g knocked out, before (ca43ae0c) and after
the R673 rule change. Copied from the R171 figure; titles changed.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12}]
# pFBA |v| at the optimum; the prolyl-tRNA fluxes are ~1e-2
PFBA_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#a7d7a0", "size": 12},
              {"type": "value", "value": 1e-3, "color": "#4fae5b", "size": 18},
              {"type": "value", "value": 1e-1, "color": "#1b6b2a", "size": 26},
              {"type": "value", "value": 10, "color": "#0b3d16", "size": 32}]
RXN_SCALES = {"pfba_wt": PFBA_SCALE, "pfba_ko_before": PFBA_SCALE, "pfba_ko_after": PFBA_SCALE}
TITLES = {
    "pfba_wt": "WT, pFBA, Task A + PO1f: growth 0.3117 in both models (the rule change does not touch the WT optimum)",
    "pfba_ko_before": "BEFORE (ca43ae0c): YALI1E05929g knocked out. R673 stays open through the OR partner YALI1D19110g: growth = WT",
    "pfba_ko_after": "AFTER: YALI1E05929g knocked out. R673 closed (the only cytosolic prolyl-tRNA synthetase step): growth 0",
}
