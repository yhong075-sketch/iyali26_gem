# Independent source audit — panel C mathematical expansion

Audited at 2026-09-16 21:24:53 UTC. Scope: the proposed CoQ9 cofactor-dilution explanation only. Read-only source inspection and exact arithmetic; no model optimization, trajectory rerun, scientific input changes, or LibreOffice use. This audit record is the only file written.

## Decision and coverage

The seven qualified claims below are supported. The equations describe a simplified two-state pool and the separate runtime sensitivity formulation, not the saved static screen in panel D or experimentally calibrated physiology.

`total 7 | audited 7 | supported 7 | unresolved 0 | contradicted 0 | unchecked 0`

This count concerns the qualified mathematical/source claims only. Biological pool calibration and transfer to real-cell growth remain outside the supported scope.

| ID | Claim | Verdict and evidence |
|---|---|---|
| C1 | With no pool drain, the two schematic rows are `v_syn - v1 + v2 = 0` and `v1 - v2 = 0`; hence `v_syn = 0`, while `v1 = v2 > 0` can satisfy these rows. | supported as algebra of the schematic. Define v1 as Q9 reduction and v2 as Q9H2 oxidation. All other model constraints still apply; this is not proof of unbounded flux or whole-network growth. |
| C2 | The runtime artificial source and dilution drain both act on oxidized Q9. | supported by the exact historical helper: `Q9_METABOLITE_ID = m468[C_mi]`; source coefficient +1 and drain coefficient -1. The matching XML names/formulas identify this as ubiquinone-9, C54H82O4, whereas m471 is ubiquinol-9, C54H84O4. |
| C3 | With the runtime extension, the two schematic rows are `v_syn + v_res - v1 + v2 - v_dil = 0` and `v1 - v2 = 0`; their sum gives `v_syn + v_res = v_dil = alpha*mu`. | supported. The helper imposes `v_dil - alpha*v_biomass = 0` with equal zero constraint bounds. Independent XML row summation again gives only native R385:+1 in the total Q9/Q9H2 pool; the runtime helper adds source:+1 and drain:-1. |
| C4 | If alpha > 0 and both net sources are closed, mu = 0 in this formulation. With synthesis closed but an available artificial reserve, reserve release can temporarily satisfy the imposed demand. | supported, conditional. Use **can**, not **must**: other model constraints can restrict or prevent growth, and the bound need not be reached. Zero model growth is not a cell-death measurement. |
| C5 | With synthesis off, a constant positive alpha, fixed volume, and no other reserve replenishment or loss, `dR/dt = -v_res*B`, `dB/dt = mu*B` imply `R + alpha*(B-B0) = R0` and `B <= B0 + R0/alpha`. | supported as the idealized reserve budget. R is the explicitly tracked **artificial reserve concentration**, not the intracellular Q9 pool. The implementation uses finite steps and clipping; see the numerical qualification below. |
| C6 | B is gDW/L, R is mmol/L, alpha is mmol/gDW, mu is h^-1, and each v is mmol/gDW/h. Alpha = 1e-4 with illustrative mu = 0.3 gives 3e-5 mmol/gDW/h. | supported by dimensions and exact arithmetic. The manifest uses alpha = 1e-4 and labels it `sensitivity_only_not_calibrated`. Mu = 0.3 is an illustrative input, not a saved run result. |
| C7 | Panel C's runtime dilution/reserve mechanism is distinct from panel D's saved static screen. | supported; retained from the preceding source audit. Runtime/reference and static-screen model identities remain unchanged and were rehashed in this pass. The preceding audit opened the saved screen driver, context loader and KO function and verified their hashes against the screen manifest; no runtime CoQ9 source/drain is added in that path. |

## Numerical qualification: ideal budget versus implementation

The matching runtime runner uses `v_res <= R_n/(B_n*dt)`, with an initially disabled helper source opened by the runner. A helper source alone does not constitute a finite reserve.

Before numerical clamping, the implemented steps are

```text
B_(n+1) = B_n * (1 + mu_n*dt)
R_(n+1) = R_n - v_res,n * B_n * dt
```

For `v_syn = 0` and `v_res,n = alpha*mu_n`, these steps preserve `R + alpha*(B-B0)` exactly in exact arithmetic. The actual code caps withdrawals by the stored reserve, applies a nonnegative lower bound, and snaps sufficiently small remaining reserves to zero. Therefore the displayed equality should be labelled **Ideal constant-alpha reserve budget**, not asserted as an exact floating-point invariant of every saved step. The existing depletion report separately preserves raw withdrawals, inventory adjustments, residuals and small biomass-cap overage. No numerical trajectory was rerun for this audit.

The equality assumes no additional reserve loss; loss can make the biomass upper bound stricter. Other metabolic or nutrient constraints may stop growth before the bound is reached. The finite-reserve relationship is conditional on synthesis remaining closed; do not apply it when R385 continues supplying the pool.

## Recommended concise English wording

- **Steady-state recycling alone** — “These two pool-balance rows allow recycling without net synthesis. Other model constraints still apply.”
- **Add growth-coupled replenishment** — “The runtime extension drains oxidized Q9 at `v_dil = alpha*mu`.”
- **Close both sources** — “For alpha > 0, `v_syn = v_res = 0` forces mu = 0 in this formulation.”
- **Finite artificial reserve** — “With synthesis off, an available reserve can temporarily meet the imposed demand.”
- Above the integrated equation: **Ideal constant-alpha reserve budget**.
- Footnote: “R is the tracked artificial reserve; alpha is uncalibrated. Other constraints may limit growth earlier.”
- Numerical example: **Illustrative calculation**, not a measurement or run output.
- Panel separation: “Runtime sensitivity extension; separate from the saved static screen in panel D.”

## Reopened sources and locators

Runtime helper, source/drain placement and linear constraint:

`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/scripts/gem_annotate/coq9_dilution.py`, lines 14–19, 69–105, 109–118.

Runtime inventory implementation:

`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/scripts/gem_annotate/quinone_dfba_essentiality.py`, lines 569–590, 631–646, 656.

Recorded scope, parameters, topology and saved evidence:

- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/coq9_supply_mechanism_20260908/run_manifest.json`, `parameters`, `inputs`, `scope`, `calibration_status`.
- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/coq9_reserve_depletion_20260908/run_manifest.json`, `inputs`, `scope`, `calibration_status`.
- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/coq9_reserve_depletion_20260908/REPORT.md`, opening questions 3–5, results table, and clipping explanation.
- `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/model.xml`, stored Q9/Q9H2 species identities and complete sum of their stoichiometric rows.

Previously audited saved static-screen path retained without recomputing results:

`/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/artifacts/screen_test_metadata_trna_20260910/run_manifest.json` and its source-hash-bound `run_screen.py`, context loader and single-gene deletion implementation.

## Exact identities rechecked in this pass

| Object | SHA256 |
|---|---|
| Runtime/reference XML | `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee` |
| Runtime dilution helper | `f223faa225c9c748ecae3d40c80d1b903e90e11fa782b32cd7946619dc8d4ff4` |
| Runtime dFBA runner | `d3810a404f6e2c93802f2efd14686f9d3c7b0336d3caa18ecc90cebbeb67344a` |
| Saved static-screen XML | `d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0` |

Both runtime manifests record the same helper and runner hashes observed here. No claim of recovering the entire historical environment is made.
