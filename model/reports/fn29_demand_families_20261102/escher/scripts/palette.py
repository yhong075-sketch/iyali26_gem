"""Colour scales and scenario titles shared by render.py and build_html.py (gap-census palette).

Adapted from the UM11 figure: the scenarios are the as-built status, pFBA and capacity, plus the capacity and pFBA of
each family's recommended in-memory what-if (each applied alone; none implemented) and the N-glycan B-only pFBA.
"""
MET_SCALE = [{"type": "value", "value": 0, "color": "#898781", "size": 12},
             {"type": "value", "value": 1, "color": "#eb6834", "size": 18},
             {"type": "value", "value": 3, "color": "#a50f15", "size": 30},
             {"type": "value", "value": 4, "color": "#8e5cc2", "size": 20}]
# capacity = max |v| over FVA at 0 % of the optimum on Task A (how much flux the reaction can carry at all)
CAP_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
             {"type": "value", "value": 1e-6, "color": "#9cc2ee", "size": 12},
             {"type": "value", "value": 1e-2, "color": "#5d9be3", "size": 16},
             {"type": "value", "value": 1, "color": "#2a78d6", "size": 20},
             {"type": "value", "value": 100, "color": "#0d366b", "size": 26}]
# pFBA flux at the optimum (|v|); the family fluxes are ~1e-5 to 2e-1
PFBA_SCALE = [{"type": "value", "value": 0, "color": "#b5b5b5", "size": 7},
              {"type": "value", "value": 1e-6, "color": "#a7d7a0", "size": 12},
              {"type": "value", "value": 1e-4, "color": "#4fae5b", "size": 18},
              {"type": "value", "value": 1e-1, "color": "#1b6b2a", "size": 26}]
# status as built (Task A): 2 blocked, 3 closed (0, 0), 4 what-if only, 5 can carry flux (pFBA 0), 6 pFBA flux
STATUS_SCALE = [{"type": "value", "value": 2, "color": "#eb6834", "size": 14},
                {"type": "value", "value": 3, "color": "#3d3d3d", "size": 14},
                {"type": "value", "value": 4, "color": "#8e5cc2", "size": 14},
                {"type": "value", "value": 5, "color": "#9cc2ee", "size": 12},
                {"type": "value", "value": 6, "color": "#4fae5b", "size": 16}]
RXN_SCALES = {"status": STATUS_SCALE, "pfba_as_built": PFBA_SCALE, "capacity_as_built": CAP_SCALE,
              "capacity_whatif": CAP_SCALE, "pfba_whatif": PFBA_SCALE, "pfba_whatif_nglycan_B": PFBA_SCALE}
TITLES = {
    "status": "AS BUILT (e1e7c868): orange = blocked, dark grey = closed (0,0), light blue = can carry flux but 0 at the pFBA optimum, green = carries pFBA flux, purple = what-if reaction NOT in the model. Task A + PO1f condition file",
    "pfba_as_built": "AS BUILT: pFBA flux at the optimum. Every family reaction is 0 except R289 / R287 (mannan via the PMT lump) and the mitochondrial Met / Ile exits R1251 / R1218",
    "capacity_as_built": "AS BUILT: capacity = max |flux| any feasible state can carry (FVA at 0 %). Blocked reactions are 0; the GPI chain and much of the sphingolipid set can carry flux but nothing demands it",
    "capacity_whatif": "WHAT-IF (NOT implemented): capacity under each family's recommended variant, each applied ALONE (sphingolipid V2, GPI b1, cardiolipin a, N-glycan C2, mitochondrial translation a2 + R1251 / R1218 reversible)",
    "pfba_whatif": "WHAT-IF (NOT implemented): pFBA flux at the optimum under each family's recommended variant, each applied ALONE; the demanded route now carries flux in every section",
    "pfba_whatif_nglycan_B": "WHAT-IF N-glycan B only (NOT implemented): R289 1:1, Golgi mannan route, Golgi Pi exit, R1137 closed; no LLO link or N-glycan term. Other sections: no data",
}
