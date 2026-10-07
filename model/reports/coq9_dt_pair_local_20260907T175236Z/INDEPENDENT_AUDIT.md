# Independent read-only audit

Status: **PASS with explicit provenance limitations**.

The independent reviewer did not run optimization and did not edit the bundle. The final saved snapshot was checked after postprocessing was frozen.

## Verified

- All eight trajectories follow the preregistered serial order and have unique condition keys.
- Input hashes match the manifest: model `bc2aac…1bee`, medium `ed176d…e85`, profile `353078…383`, runner `d3810a…344a`, and reused r3 metric code `0bd322…a066`.
- The model file hash is unchanged before/after execution, and all in-memory model restoration checks passed.
- All eight solver readbacks are `Threads=1`, `Seed=0`, `Method=1`, `FeasibilityTol=1e-9`; the recorded additional values are mutually consistent.
- The 2,810 trajectory rows have exact time steps, continuous recorded states, nonnegative finite states, and four valid zero-duration terminal-only `infeasible` rows. Invalid terminal fluxes remain `NaN`.
- Same-dt metrics and unified-common-window metrics independently recalculate to the stored values.
- All 19 files that preceded this audit note were fully inventoried and passed SHA256 verification; this note is added to the same ledger below.

## Independently confirmed results

- `finite_batch`: `T*=24 h`; maximum endpoint doubling delta `2.39808e-14`; maximum dimensionless ratio delta `0.000745411`; maximum event-time delta `0.125 h`. The strict diagnostic fails because the event delta exceeds `0.0625 h` and the WT uracil exact-zero event is present at only one dt.
- `po1f_nonlimiting`: `T*=5.1875 h`; maximum endpoint doubling delta `0.1770026324`; maximum dimensionless ratio delta `5.21583e-13`; maximum event-time delta `0.125 h`. The strict diagnostic fails the doubling and event-time thresholds.
- Finite-batch KO uracil exact-zero occurs at `10.375/10.25 h`; finite-batch WT is `NA/3.40625 h`; nonlimiting glucose exact-zero is `5.3125/5.1875 h` for both WT and KO.
- Finite-batch KO initial-growth ratio is `0.3288438295`, 24 h net-increment ratio is numerically 1, and biomass-AUC ratios are `0.7451716691/0.7459170802`.
- Every trajectory keeps the artificial Q9 reserve at `1e-6 mmol/L`; all source-positive interval counts are zero. Derived reserve-enabled second-solve attempts are `328,219,0,0,659,441,0,0`. These are runner-policy inferences, not native backend LP counts.

## Scope and provenance limitations

- Native Gurobi backend optimization counts were not instrumented before execution and remain unknown. The 4,461 stored runner-level calls are pFBA/fallback calls, not native LP counts.
- The summed `simulate_gene` wall time is `357.326139 s`; it excludes fresh-context model loading.
- Effective-context provenance is stored once because every loaded context matched; only the model file was rehashed after execution.
- The driver and plots were corrected during read-only postprocessing. Their final hashes prove the delivered implementation, not the exact original orchestration-script bytes at solve start. No trajectory was rerun.
- The figures were visually and structurally checked. PO1f glucose-zero/solver termination events are described in the report/table rather than marked directly on the curves.
- All conclusions remain `runtime_only` and `sensitivity_only_not_calibrated`.

Gene identity wording checked:

- **YALI1A14736g — no established gene name — uncharacterized protein (model/GPR assignment only; R305).**
- **YALI1A21711g — no established gene name — uncharacterized protein (model/GPR assignment only; R2062).**
