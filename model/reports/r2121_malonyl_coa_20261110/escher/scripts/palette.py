"""Colour scales and scenario titles shared by render.py and build_html.py (gap-census palette).

Adapted from the UM7 nuclear-diffusion figure: pFBA before / after and an on-off status scenario.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12},
             {"type": "value", "value": 1, "color": "#898781", "size": 12}]  # Escher needs two points; all nodes are 0
FLUX_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#9cc2ee", "size": 12},
              {"type": "value", "value": 1e-4, "color": "#5d9be3", "size": 15},
              {"type": "value", "value": 1e-2, "color": "#2a78d6", "size": 19},
              {"type": "value", "value": 0.1, "color": "#0d366b", "size": 24},
              {"type": "value", "value": 1.0, "color": "#6250d6", "size": 28}]
# status: 1 flux before only, 2 after only, 3 both, 4 neither, 5 closed after
STATUS_SCALE = [{"type": "value", "value": 1, "color": "#2a78d6", "size": 18},
                {"type": "value", "value": 2, "color": "#6250d6", "size": 18},
                {"type": "value", "value": 3, "color": "#9cc2ee", "size": 12},
                {"type": "value", "value": 4, "color": "#b5b5b5", "size": 8},
                {"type": "value", "value": 5, "color": "#3d3d3d", "size": 16}]
RXN_SCALES = {"before": FLUX_SCALE, "after": FLUX_SCALE, "status": STATUS_SCALE}
TITLES = {
    "before": "BEFORE: e584151f, pFBA flux at the optimum on Task A + PO1f. R2121 makes the cytosolic malonyl-CoA; TAG is broken down, beta-oxidised and re-elongated",
    "after": "AFTER: R2121 closed, pFBA flux at the optimum on Task A + PO1f. Acetyl-CoA carboxylase R88 makes the malonyl-CoA; the TAG cycle stops",
    "status": "STATUS: blue = flux before only, purple = after only, light blue = both, grey = neither, black = closed after (R2121)",
}
