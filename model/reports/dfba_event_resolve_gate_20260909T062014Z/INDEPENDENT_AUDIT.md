# Independent post-run audit

Two independent read-only reviewers separately audited the completed outputs.
Neither reviewer edited files or ran a GEM solve.

## Coverage

- Frozen execution identities: `3/3` matched.
- Output hashes recorded in the manifest: `3/3` matched.
- Trajectory segments: `86/86` inspected.
- Inventory records: `1204/1204` inspected (`86 x 14`).
- Depletion events and same-time re-solves: `5/5` linked.
- Backend budget: `172/256`; wall budget: `12.804310625011567/600 s`.

## Audit conclusion

**Execution integrity PASS.** The trajectory has 85 optimal advanced segments
and one final zero-duration `infeasible` record. Time, biomass, and every
inventory are continuous. All raw balance residuals and clips satisfy the
frozen tolerances. Persisted numeric NaNs occur only where fluxes are undefined
in the permitted terminal record. Q9 source usage is zero throughout the
advanced trajectory. Output hashes match the manifest, and the model file,
complete disposable model state, and solver hook are restored or unchanged.

**Scientific gate FAIL.** The glucose-zero event occurred at
`5.048080683449264 h`, which is within `0.0625 h` of the stored quarter-step
event, but the candidate terminated before `T*=5.109375 h`. Therefore the
candidate value at `T*` is unavailable without forbidden extrapolation, and the
doubling-gap criterion is false. This is a prespecified numerical/scientific
failure, not an infrastructure failure or evidence of biological death.

The independent conclusions agree with the manifest and with a separate
read-only full-row recomputation of continuity, exposure balances, event links,
budgets, hashes, and restoration fields.

