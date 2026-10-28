"""Colour scales and scenario titles shared by render.py and build_html.py (gap-census palette).

Adapted from the dolichol-route figure: capacity scenarios (FVA at 0 % on Task A) and a status scenario replace the
census one, because pFBA carries no flux through the UM7 reactions at the optimum in either model.
"""
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
# capacity = max |v| over FVA at 0 % of the optimum on Task A (how much flux the reaction can carry at all)
CAP_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
             {"type": "value", "value": 1e-6, "color": "#9cc2ee", "size": 12},
             {"type": "value", "value": 1e-2, "color": "#5d9be3", "size": 16},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 100, "color": "#0d366b", "size": 26}]
# status (Task A, whole-model FVA at 0 %): 1 freed, 2 blocked in both, 3 closed (0, 0) in both, 4 new and unblocked,
# 5 unblocked in both
STATUS_SCALE = [{"type": "value", "value": 1, "color": "#2a78d6", "size": 18},
                {"type": "value", "value": 2, "color": "#eb6834", "size": 14},
                {"type": "value", "value": 3, "color": "#3d3d3d", "size": 14},
                {"type": "value", "value": 4, "color": "#6250d6", "size": 18},
                {"type": "value", "value": 5, "color": "#9cc2ee", "size": 12}]
RXN_SCALES = {"capacity_before": CAP_SCALE, "capacity_after": CAP_SCALE, "before": FLUX_SCALE, "after": FLUX_SCALE,
              "status": STATUS_SCALE}
TITLES = {
    "capacity_before": "CAPACITY BEFORE: 9d655c2b, max |flux| any feasible state can carry (FVA at 0 %), Task A + PO1f. Grey = blocked",
    "capacity_after": "CAPACITY AFTER: UM7 candidate, max |flux| (FVA at 0 %), Task A + PO1f. Blue = can carry flux; R_ADPRtn and R_PPItn are new",
    "before": "BEFORE: 9d655c2b, pFBA flux at the optimum on Task A + PO1f. Grey = no flux",
    "after": "AFTER: UM7 candidate, pFBA flux at the optimum on Task A + PO1f (the optimum does not use the new links)",
    "status": "STATUS (Task A, whole-model FVA): blue = freed, orange = blocked in both, black = closed (0, 0), purple = new and unblocked, light blue = unblocked in both",
}
