# Maximum growth at maximum explicit TAG output

The corrected lipid-unlump candidate returns maximum TAG mass output **0.606736134771869 g TAG/gDW/h**. Fixing that output, maximizing growth returns **0.14650760568105112 h^-1**, numerically equal to the preserved growth floor (0.14650760568106092 h^-1). The associated doubling time is approximately **4.73 h**. This is lexicographic optimization: TAG output has priority, growth is optimized second. It is not simultaneous maximization of two objectives or maximization of a product/substrate ratio.

The fixed input XML is the charge-corrected candidate (SHA-256 50f086e5d12d80c09f3c3fc989e7198dbe7e01109f247b46fdeb82cb53974433), built by the correction subsequently pushed as commit 6868777. All 343 explicit TAGs were verified neutral before optimization. The TAG mass weights, freely varying composition, SD-Leu/PO1f effective bounds and absolute growth floor match the earlier experiment. Previous effective bounds were transferred by reaction ID after exact stoichiometric equality checks. A first preflight stopped before any solve because the old runtime configuration changes R45 to its plasmid pseudo-gene; that exact recorded overlay was then restored. The historical context fingerprint is retained as historical provenance, not assigned to the corrected model.

## Results and validation

| Stage | Result | Saved-point assessment |
|---|---|---|
| Maximize TAG output | 0.606736134771869 g/gDW/h | Passes 1e-7 feasibility checks |
| Maximize growth with TAG fixed at the above maximum | 0.14650760568105112 h^-1 | Passes 1e-7 checks; TAG equality discrepancy 2.56e-13 |
| Minimize split flux on the maximum TAG face | Solver objective 526.3350493187479 | Rejected: bound violation 1.5004471833076028e-7 exceeds 1e-7 |
| Maximize growth with a new 1.25-fold pFBA cap | Not executed | Stopped after preceding validation failure |

For both accepted stages, the maximum mass residual is 8.18e-8 and maximum split-bound violation is 9.55e-8. These approach the prescribed tolerance; optimum values are solver-reported, tolerance-limited results, not exact-arithmetic certificates. `verify.py` independently recomputed net mass/bound residuals from saved JSON and NPZ without solving and confirmed the pass/fail decisions. No values were clipped.

The maximum-growth solution's glucose uptake is 10 mmol/gDW/h, and its only TAG demand above 1e-7 is tripalmitin at about 0.751544 mmol/gDW/h. This is one returned optimal mixture; it does not establish unique composition or physiological lipid composition. Extra intracellular TAG demand is an accumulation proxy, not a modeled secretion mechanism. Units assume the model's mmol/gDW/h convention.

The optimum growth remaining at the lower bound indicates that, within this constrained model and numerical precision, there is no additional growth available while retaining the maximum TAG output. It supports a growth/product competition interpretation in this model; it is not an experimental growth prediction. The chemical correction leaves the FBA problem unchanged, and the small difference from the older reported TAG maximum is not evidence of increased production capability.

A cap restricts the feasible region and cannot improve the exact mathematical maximum growth. Nevertheless, no new capped maximum-growth solve was completed here, so capped feasibility and its optimum are not reported as tested results. No tolerance was relaxed, no TAG target was lowered and no CFSA sampling was run.

## Execution and evidence

Three backend optimizations were used against the four-call/300-second budget. Primary requested optimization completed; supplemental cap comparison stopped. `completed_result.json` retains status `stopped` rather than hiding the supplemental failure. The earlier `result.json` and `preflight_*` files preserve the zero-solve preflight. Full inputs, software, code hashes, runtime bounds, objective constraints and split/net solutions are retained in the accompanying JSON/NPZ files. Existing models, historical outputs and activation gates were preserved.
