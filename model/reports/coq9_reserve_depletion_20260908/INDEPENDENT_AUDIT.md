# Independent read-only audit

## Verdict

**PASS for the bounded runtime implementation question after reporting corrections.** Two independent read-only reviews ran no solver and edited no files. They audited all 7 final ledger claims: 7 supported within their stated scope, 0 contradicted, 0 unchecked. This is not biological validation.

## Recomputed result checks

- Structure: 32/32 advanced intervals, steps 0–31, exact 0–1 h grid at `dt=0.03125 h`; time, biomass, Q9-reserve and glucose continuity residuals are 0.
- Budget inputs: `alpha=1e-4 mmol/gDW`, `B0=0.01 gDW/L`, `Q0=1.0000000000000002e-6 mmol/L`, initial source cap `0.0032 mmol/gDW/h`, encoded biomass cap `0.020000000000000004 gDW/L`.
- R385 is exactly 0 in all selected intervals. Max absolute stored source-demand residual is `3.5910402256e-15 mmol/gDW/h`; dilution-demand residual is 0. Reconstructed from the stored TSV, max absolute Q9-, Q9H2-, and total-Q-row residuals are `3.03708e-15`, `3e-15`, and `3.59104e-15 mmol/gDW/h`.
- Euler biomass-update and stored reserve-update residuals are 0. Max local and cumulative budget residuals are `2.0998740856e-18` and `3.9700434238e-18 mmol/L`.
- Raw source integral is `1.0000000000012337e-6 mmol/L`; actual inventory deduction is `1.0000000000000002e-6 mmol/L`. At step 15 (`0.46875–0.5 h`), raw pre-cap inventory is `-1.2333196759e-18 mmol/L`; the inventory cap yields saved zero. There is no positive `<=1e-12 mmol/L` snap loss.
- Final biomass is `0.0200000000000397 gDW/L`, an overage of `3.9697412024e-14 gDW/L` above the encoded cap, within the prespecified `1e-12` check tolerance but not mathematical exactness.
- Steps 16–31 start after reserve depletion and have selected solver status `optimal` with raw biomass and source fluxes exactly 0. Raw and clipped values agree throughout; no hidden tiny positive or negative growth was found.
- Glucose remains `110.92822956102589 mmol/L`; no glucose-zero event or infeasible interval occurs. Uracil inventory is intentionally NaN in `po1f_nonlimiting`. Twelve other finite medium pools have no stored dynamic inventories, so broader exclusion of every nutrient limitation is not supported.
- Solver accounting: 64/64 pFBA calls and 128/128 backend optimizations are optimal; 0 fallback, 0 exception; wall time `3.542641874984838 s`.
- Numerical observations retained without changing tolerance: step 4 R558 bound violation `1.7813790807e-10`, below `FeasibilityTol=1e-9`; depletion-row source-cap overage `2.0164940681e-15 mmol/gDW/h`.

## Provenance audit

- Reused package checksum ledger independently passes 12/12 entries. Old conditions 04 and 05 match their cited mode, background, R385/source settings, solver status and growth class; neither was rerun.
- Current manifest input hashes match 10/10 files; effective context matches the reused package 7/7; solver readback covers 9/9 frozen settings.
- Persistent canonical model SHA remains `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`. The manifest records in-memory restoration as true, but its before/after fingerprints were not saved separately, so that restoration assertion is not independently reconstructable after the run.
- Code identity records the direct runner/helpers but not every transitive local import; complete historical-environment replay remains limited by that omission.
- The final report and ledger explicitly retain `runtime_only / sensitivity_only_not_calibrated`, the old R558 exception, and the old STOP/time-step conclusions. No model, GPR or project-state claim is authorized by this audit.

## Decision

The stored evidence supports one narrow statement: **with R385 closed, the existing artificial finite-reserve implementation consumes its budget and is sufficient to force solver-optimal zero growth after the saved reserve reaches zero in this chosen runtime trajectory.** It does not establish a measured CoQ9 pool, exclusive biological causation, cell death, formal GPR support, or model calibration. No extension run is warranted for this task.
