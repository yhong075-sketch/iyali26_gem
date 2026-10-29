"""Colour scales and scenario titles shared by render.py and build_html.py (gap-census palette).

Adapted from the UM2 figure: the model is not changed (Option B), so the scenarios are the as-built status, pFBA and
capacity, plus the capacity and pFBA of the in-memory Option A what-if (thioredoxin links), which was not implemented.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 2, "color": "#eb6834", "size": 20},
             {"type": "value", "value": 3, "color": "#a50f15", "size": 30}]
# capacity = max |v| over FVA at 0 % of the optimum on Task A (how much flux the reaction can carry at all)
CAP_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
             {"type": "value", "value": 1e-6, "color": "#9cc2ee", "size": 12},
             {"type": "value", "value": 1e-2, "color": "#5d9be3", "size": 16},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 100, "color": "#0d366b", "size": 26}]
# pFBA flux at the optimum (|v|); the RNR fluxes are ~1e-3 to 3e-2
PFBA_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#a7d7a0", "size": 12},
              {"type": "value", "value": 1e-3, "color": "#4fae5b", "size": 18},
              {"type": "value", "value": 1e-1, "color": "#1b6b2a", "size": 26}]
# status as built (Task A, whole-model FVA at 0 %): 2 blocked, 3 closed (0, 0), 4 what-if only, 5 can carry flux
STATUS_SCALE = [{"type": "value", "value": 2, "color": "#eb6834", "size": 14},
                {"type": "value", "value": 3, "color": "#3d3d3d", "size": 14},
                {"type": "value", "value": 4, "color": "#8e5cc2", "size": 14},
                {"type": "value", "value": 5, "color": "#9cc2ee", "size": 12}]
RXN_SCALES = {"status": STATUS_SCALE, "pfba_as_built": PFBA_SCALE, "capacity_as_built": CAP_SCALE,
              "capacity_whatif_A": CAP_SCALE, "pfba_whatif_A": PFBA_SCALE}
TITLES = {
    "status": "AS BUILT (decision B, unchanged): orange = blocked (the 11 UM11 reactions), light blue = can carry flux, purple = what-if link NOT in the model. Task A + PO1f",
    "pfba_as_built": "AS BUILT: pFBA flux at the optimum, Task A + PO1f. Deoxynucleotides are made in the cytosol (R702, R710, R1913 at this optimum); every nuclear copy is 0",
    "capacity_as_built": "AS BUILT: capacity = max |flux| any feasible state can carry (FVA at 0 %), Task A + PO1f. Every UM11 reaction is 0",
    "capacity_whatif_A": "WHAT-IF A (NOT implemented): reduced and oxidized thioredoxin linked nucleus <-> cytosol in memory. All 11 gain capacity; nothing else in the model changes",
    "pfba_whatif_A": "WHAT-IF A (NOT implemented): pFBA flux at the optimum. Every nuclear copy is still 0: they add capacity, never growth",
}
