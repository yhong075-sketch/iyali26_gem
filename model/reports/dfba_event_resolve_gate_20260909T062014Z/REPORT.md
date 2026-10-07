# WT event-and-resolve real-GEM gate

## Decision

- Execution integrity: **PASS**.
- Prespecified scientific gate: **FAIL**.
- Scope: `runtime_only` and `sensitivity_only_not_calibrated`.
- Production adoption: **not authorized**.

This was the single authorized WT trajectory. It was not retried. The frozen
condition was PO1f SD-Leu, `po1f_nonlimiting`, alpha `1e-4 mmol/gDW`, pool
multiplier `1`, outer step `0.0625 h`, and a `6 h` horizon. The candidate used
physical instantaneous medium bounds, one exact biomass exposure for biomass
and all 14 inventories, and an immediate same-time re-solve after each
depletion event.

## Prespecified result

The candidate reached an `infeasible` terminal re-solve at
`5.048080683449264 h`, immediately after the glucose inventory reached zero.
The common comparison time was `T*=5.109375 h`. Because the candidate ended
before `T*`, biomass and doublings at `T*` remain null; no extrapolation was
performed, and the doubling-gap criterion could not be evaluated.

The glucose-event subcriterion passed: its absolute difference from the stored
quarter-step event time was `0.06129431655073603 h`, within the prespecified
`0.0625 h` tolerance. That subcriterion alone is insufficient, so the combined
gate remains FAIL.

The five recorded early events were R1003, R1215, R1217, R1202, and R1070
(glucose). The first four same-time re-solves remained optimal; the glucose
event produced the permitted zero-duration `infeasible` terminal record.
`Infeasible` is a solver state and is not biological death evidence.

## Execution checks

- Exactly one GEM trajectory, with 85 advanced segments and one terminal
  record; zero trajectory retries.
- All 86 segments contain exactly 14 inventory records (`1204/1204`).
- All five events have an immediate same-time re-solve.
- Time, biomass, and inventory continuity passed.
- Maximum biomass-balance residual: `1.66e-15`.
- Maximum raw inventory-balance residual: `6.989200884710556e-15 mmol/L`.
- No material clamp; Q9 source flux was exactly zero in all 85 optimal segments.
- Backend entries: `172/256`; wall time: `12.804310625011567/600 s`.
- Model fingerprint, full disposable model state, solver hook, and model file
  were restored or unchanged.

## Interpretation and stop decision

The event-and-resolve candidate executed cleanly, but it failed the frozen
comparison gate because it did not reach `T*`. This result does not establish
time-step convergence and does not authorize refinement, production wiring,
parameter changes, or changes to `model.xml`, GPRs, bounds, stoichiometry,
curated data, or the FN dossier. Under the preregistration, this candidate path
is rejected or remains unresolved unless a separately authorized study is
defined.

## Frozen identities

- Driver SHA-256: `c9c8a7f4de9f9c44329e15fd5d20e13a08287cdbafce51e1b2a7634a438eb9fd`
- Preregistration SHA-256: `03b3f72641fea8132150fb61c33e489d5fe0d76eee987e0395a993d009e87fba`
- Candidate SHA-256: `596d6908f11dab0cda876d4703cf7bf31f0dcee0089de1b5ae89767666aeb440`
- Model SHA-256: `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`

