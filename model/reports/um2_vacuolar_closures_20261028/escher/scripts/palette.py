"""Colour scales and scenario titles shared by render.py and build_html.py (gap-census palette).

Adapted from the UM7 figure: the model is not changed (Option A0), so the scenarios are the as-built status and
capacity plus two in-memory what-if options that were not implemented, and the no-uptake loop of the census test.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 2, "color": "#eb6834", "size": 20}]
# capacity = max |v| over FVA at 0 % of the optimum on Task A (how much flux the reaction can carry at all)
CAP_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
             {"type": "value", "value": 1e-6, "color": "#9cc2ee", "size": 12},
             {"type": "value", "value": 1e-2, "color": "#5d9be3", "size": 16},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 100, "color": "#0d366b", "size": 26}]
# loop = max |v| with no uptake at all (every exchange (0, 1000), biomass 0): any non-zero value is an internal cycle
LOOP_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#f0a3a3", "size": 14},
              {"type": "value", "value": 100, "color": "#d03b3b", "size": 22},
              {"type": "value", "value": 1000, "color": "#8a1414", "size": 28}]
# status as built (Task A, whole-model FVA at 0 %): 3 closed (0, 0), 2 open but blocked, 5 can carry flux
STATUS_SCALE = [{"type": "value", "value": 2, "color": "#eb6834", "size": 14},
                {"type": "value", "value": 3, "color": "#3d3d3d", "size": 14},
                {"type": "value", "value": 5, "color": "#9cc2ee", "size": 12}]
RXN_SCALES = {"status": STATUS_SCALE, "capacity_as_built": CAP_SCALE, "capacity_whatif_cprime": CAP_SCALE,
              "loop_whatif_open17": LOOP_SCALE}
TITLES = {
    "status": "AS BUILT (decision A0, unchanged): black = closed (0, 0) by inheritance, orange = open but blocked, light blue = can carry flux. Task A + PO1f",
    "capacity_as_built": "AS BUILT: capacity = max |flux| any feasible state can carry (FVA at 0 %), Task A + PO1f. Every UM2 reaction is 0",
    "capacity_whatif_cprime": "WHAT-IF C' (NOT implemented): pump, water, R1163 as a vac->cyt leak and the 11 AVT/ATG22 transports opened in memory. Capacity on Task A: futile ATP cycles only, pFBA 0",
    "loop_whatif_open17": "WHAT-IF 'open all 17' with R1163 as written (NOT implemented): flux with NO uptake at all. Red = internal loop (H+ in by R1163, out with the amino-acid shuttles)",
}
