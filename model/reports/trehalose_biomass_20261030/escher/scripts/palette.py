"""Colour scales and scenario titles shared by render.py (gap-census palette, adapted from the UM11 figure).

Scenarios: pFBA at WT and with YALI1E17616g (TPS1) knocked out, each before (e1e7c868) and after the trehalose
biomass term was removed. Adapted from the aconitase figure; scenarios, titles and comment changed.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12}]
# pFBA |v| at the optimum; glycolytic entry ~2, the trehalose branch ~2e-3
PFBA_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#a7d7a0", "size": 12},
              {"type": "value", "value": 1e-3, "color": "#4fae5b", "size": 18},
              {"type": "value", "value": 1e-1, "color": "#1b6b2a", "size": 26},
              {"type": "value", "value": 10, "color": "#0b3d16", "size": 32}]
RXN_SCALES = {"pfba_wt_before": PFBA_SCALE, "pfba_wt_after": PFBA_SCALE, "pfba_ko_before": PFBA_SCALE,
              "pfba_ko_after": PFBA_SCALE}
TITLES = {
    "pfba_wt_before": "BEFORE (e1e7c868), WT: growth 0.31166; R143 and R767 carry 0.00214 into the biomass_C trehalose term",
    "pfba_wt_after": "AFTER, WT: growth 0.31230 (+0.20 %); no trehalose term, the branch carries no flux",
    "pfba_ko_before": "BEFORE, YALI1E17616g (TPS1) knocked out: R143 closed, trehalose unreachable, growth 0",
    "pfba_ko_after": "AFTER, YALI1E17616g (TPS1) knocked out: R143 closed, growth 0.31230 = WT",
}
# Side-by-side PNGs written by render.py: (output name, left scenario, right scenario)
COMPARES = [("compare_wt", "pfba_wt_before", "pfba_wt_after"), ("compare_ko", "pfba_ko_before", "pfba_ko_after")]
