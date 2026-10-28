"""Colour scales and scenario titles shared by render.py and build_html.py (gap-census palette)."""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 2, "color": "#eb6834", "size": 20},
             {"type": "value", "value": 3, "color": "#e87ba4", "size": 20}]
FLUX_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#9cc2ee", "size": 12},
              {"type": "value", "value": 1e-5, "color": "#5d9be3", "size": 15},
              {"type": "value", "value": 1e-4, "color": "#2a78d6", "size": 19},
              {"type": "value", "value": 1e-3, "color": "#0d366b", "size": 24},
              {"type": "value", "value": 0.1, "color": "#6250d6", "size": 28}]
CENSUS_SCALE = [{"type": "value", "value": 1, "color": "#eda100", "size": 16},
                {"type": "value", "value": 2, "color": "#e87ba4", "size": 16},
                {"type": "value", "value": 3, "color": "#2a78d6", "size": 16}]
RXN_SCALES = {"before": FLUX_SCALE, "after": FLUX_SCALE, "census": CENSUS_SCALE}
TITLES = {
    "before": "BEFORE: a1237821 (a12378210c7c6672...), pFBA on Task A + PO1f. Grey = no flux; no freed reaction carries flux. biomass_C (light grey) lacks these six terms here",
    "after": "AFTER: 00e8e397 (00e8e39775b64cfe...), pFBA on Task A + PO1f. Blue = flux (darker = larger)",
    "census": "CENSUS: where each reaction was freed (amber = open only, pink = Task A only, blue = both)",
}
