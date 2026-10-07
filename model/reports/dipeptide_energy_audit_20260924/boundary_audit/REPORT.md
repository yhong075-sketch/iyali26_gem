# Independent B1 boundary audit

Checked on 2026-09-24 by `/root/energy_boundary_audit`, within the current workspace only. This audit reads the fixed XML, historical complete-run effective state, loader source and configurations; it does not call the solver or modify a scientific model. `check.py` is the executable static check; complete identities and reaction lists are in `facts.json`.

## Closure recommendation

The fixed XML SHA256 remains `aad701126d12d113816fda4b872333b614b4469ee1b6d8ab8c419231c89e965f`. Historical effective state has 2,314 reactions, 1,877 material constraints, and exactly two nonnegative solver variables per reaction. All 1,877 rows are homogeneous steady-state constraints; no extra/custom constraint, nonzero row bound or SBML boundary-condition species was found. This statement is about the inspected fixed/historical state, not an independently rerun loader.

Close both directions of all 183 one-sided reactions, independent of compartment, ID prefix or current medium: 182 external one-metabolite exchanges plus internal biomass drain R1373. There are no empty reactions or additional multispecies one-sided reactions in this state. Explicitly close all six producers of the shared biomass species `m1401[C_cy]`: R1372, R1387, R1710, xBIOMASS, **newBiom**, and biomass_C. Finding newBiom requires tracing biomass-species adjacency, not just matching `biomass` in IDs/names. The resulting explicit list has 189 reactions.

Strictly, when R1373 is closed, its complete biomass-species row already forces all six producer fluxes to zero because they all have positive coefficients and nonnegative lower bounds. Explicit closure therefore clarifies the diagnostic without silently relaxing or deleting any material row. Each of the 20 TRNA_BIOMASS residue-production reactions is also forced to zero by its own complete two-entry row when biomass_C is closed. Likewise xLIPID and xAMINOACID are forced to zero by their product pools when xBIOMASS is closed. They are not unrecorded free sources.

Eight reversible xPOOL_AC/FA columns are two-sided internal mixing reactions, with real constituent metabolites on one side and an aggregate on the other. They should be enumerated and their chemical incompleteness recorded, not automatically classified as external inputs merely because of their names. Their generic/absent aggregate formula prevents full elemental verification. Closing them can be an explicitly labelled additional conservative test; it must not be confused with identifying the cause of an internal ATP defect. All eight, xLIPID and xAMINOACID have zero flux in the old maintenance-max witness inspected here.

No DIAG source/drain from the earlier supply scenarios exists in the baseline effective reaction set. When constructing new diagnostic copies, assert the absence of those temporary columns, and explicitly close any additional artificial source or drain rather than relying on `model.medium = {}`. Do not close internal reactions merely because a stored formula is incomplete or a chemical residual is nonzero: that would pre-remove potential mechanisms before the witness audit.

## Maintenance, compulsory fluxes and zero feasibility

The only positive lower bound is xMAINTENANCE `[7.8625,1000]`; no reaction has negative upper bound. Keep xMAINTENANCE and set its diagnostic bounds to `[0,1000]`. Its exact stored forward chemistry is:

`ATP[C_cy] + H2O[C_cy] -> ADP[C_cy] + phosphate[C_cy]`.

The stored formulae are C10H16N5O13P3, H2O, C10H15N5O10P2 and H3O4P; all stored charges are 0. The maintenance column is element- and charge-balanced under that convention. Do not add a proton-bearing surrogate dissipation column: the preceding memote surrogate had a different proton convention. Formula/charge consistency here does not assert physiological predominance of neutral ATP.

After the stated closure and maintenance-LB reset, all reaction bounds contain 0 and all constraints are homogeneous. The all-zero flux vector is therefore an exact algebraic feasible point of the inspected state. A new-loader zero-feasibility failure would require examination of state differences or implementation, not a biological-death interpretation. All internal reaction directions and all material constraints must remain unchanged.

## Historical witness scope

The previous complete run contains a saved `closed_local/baseline/xMAINTENANCE/max` witness with maintenance 1000 and 11 nonzero columns: R211, R419, R603 (reverse), R632, R694, R778, xMAINTENANCE, R2013, R2092, R_NTP7 (reverse), and R_PGAM1_PhosHydro. None of the proposed additional biomass or pool closures excludes that saved witness. This is a read-only observation, not a new optimized witness and not proof that any one listed reaction is the cause or that the support is minimal. Chitin's stored formula `H2O(C8H13NO5)n` is a generic polymer formula; it is not a fully specified molecular formula. The upcoming fixed-unit witness and directional blocking tests must retain this distinction.
