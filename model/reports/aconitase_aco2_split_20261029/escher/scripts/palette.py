"""Colour scales and scenario titles shared by render.py (gap-census palette, adapted from the UM11 figure).

Scenarios: pFBA at WT (same growth in both models) and pFBA with YALI1D11984g knocked out, before (059316ef) and after
the aconitase rule change. Copied from the R673 figure; titles and comment changed.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12}]
# pFBA |v| at the optimum; the TCA fluxes are ~1, the residual knockout fluxes ~1e-2
PFBA_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#a7d7a0", "size": 12},
              {"type": "value", "value": 1e-3, "color": "#4fae5b", "size": 18},
              {"type": "value", "value": 1e-1, "color": "#1b6b2a", "size": 26},
              {"type": "value", "value": 10, "color": "#0b3d16", "size": 32}]
RXN_SCALES = {"pfba_wt": PFBA_SCALE, "pfba_ko_before": PFBA_SCALE, "pfba_ko_after": PFBA_SCALE}
TITLES = {
    "pfba_wt": "WT, pFBA, Task A + PO1f: growth 0.3117 in both models (the rule change does not touch the WT optimum)",
    "pfba_ko_before": "BEFORE (059316ef): YALI1D11984g knocked out. R238/R216/R1556 stay open through the OR partner YALI1E17935g; only R239 closes: growth = WT",
    "pfba_ko_after": "AFTER: YALI1D11984g knocked out. All four aconitase steps closed: growth 0.00761 (KO/WT 0.0244)",
}
