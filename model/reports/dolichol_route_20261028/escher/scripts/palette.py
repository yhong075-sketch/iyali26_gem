"""Colour scales and scenario titles shared by render.py and build_html.py (gap-census palette)."""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 2, "color": "#eb6834", "size": 20},
             {"type": "value", "value": 3, "color": "#e87ba4", "size": 20}]
FLUX_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#9cc2ee", "size": 12},
              {"type": "value", "value": 1e-5, "color": "#5d9be3", "size": 15},
              {"type": "value", "value": 1e-4, "color": "#2a78d6", "size": 19},
              {"type": "value", "value": 1e-2, "color": "#0d366b", "size": 24},
              {"type": "value", "value": 0.2, "color": "#6250d6", "size": 28}]
# census (open condition): 1 freed, 2 blocked in both, 3 closed (0, 0) after, 4 new and unblocked, 5 unblocked in both
CENSUS_SCALE = [{"type": "value", "value": 1, "color": "#2a78d6", "size": 18},
                {"type": "value", "value": 2, "color": "#eb6834", "size": 14},
                {"type": "value", "value": 3, "color": "#3d3d3d", "size": 14},
                {"type": "value", "value": 4, "color": "#6250d6", "size": 18},
                {"type": "value", "value": 5, "color": "#9cc2ee", "size": 12}]
RXN_SCALES = {"before": FLUX_SCALE, "after": FLUX_SCALE, "census": CENSUS_SCALE}
TITLES = {
    "before": "BEFORE: ce9e1112 (ce9e1112ee940531...), pFBA on Task A + PO1f. Grey = no flux; R_DOLSYN not in this model; biomass_C (light grey) lacks the Dol-P term",
    "after": "AFTER: 9d655c2b (9d655c2bdbf82721...), pFBA on Task A + PO1f. Blue = flux (darker = larger); R2178 closed (0, 0)",
    "census": "CENSUS (open condition): blue = freed, orange = blocked in both, black = closed (0, 0), purple = new and unblocked, light blue = unblocked in both",
}
