# N2 — consistent-exposure dFBA diagnosis and isolated candidate

## Decision at a glance

The stored WT trajectory is **not time-step converged**: at the common physical time \(T^*=5.109375\) h, successive doubling differences are 0.225522 and 0.109728, both above the declared 0.01 line. The differences shrink, but not enough to call convergence.

The current production runner is nevertheless internally consistent explicit Euler: biomass, finite nutrients, and artificial Q9 reserve all use the same left-endpoint exposure \(B_n\Delta t\). The numerical weakness is first-order exposure and event timing, not an update that integrates biomass one way and nutrients another way.

This package adds one default-off, pure-standard-library event kernel and no production wiring. It is **toy-only, not GEM-verified**. Scope remains `runtime_only` and `sensitivity_only_not_calibrated`.

Actual budget: zero GEM solves; three deterministic toy examples; no HPCC, grid, scan, model, GPR, bounds, labels, benchmark, or production-runner change.

## Reused evidence and common condition

The analysis reuses three checksum-verified WT trajectories at \(\Delta t=0.0625\), 0.03125, and 0.015625 h. No matching prior N2 package with this objective, inputs, conditions, and method was available in the inspected artifact set, so this is a first execution rather than a rerun.

The common condition is `po1f_nonlimiting`, no added runner KO within the enabled `po1f_sd_leu_accrispr_v1` overlay, \(\alpha=10^{-4}\) mmol/gDW, pool multiplier 1, initial biomass 0.01 gDW/L, initial glucose 111 mmol/L, `runtime_topology=false`, and no runtime GPR scenario. Q9 source flux is zero in every stored interval; the artificial Q9 reserve remains 1.0000000000000002e-6 mmol/L and never depletes. Thus Q9 is not the event causing this WT comparison.

| dt (h) | terminal time/status | biomass at \(T^*\) (gDW/L) | doublings at \(T^*\) | glucose at \(T^*\) (mmol/L) |
|---:|---|---:|---:|---:|
| 0.0625 | 5.3125 h / infeasible | 12.893902045456 | 10.332473213197 | 23.043547309470 |
| 0.03125 | 5.1875 h / infeasible | 15.075557240933 | 10.557995614408 | 8.149716906919 |
| 0.015625 | 5.109375 h / infeasible | 16.266899421744 | 10.667723574692 | 0 |

Endpoint values at their different terminal times are not used as a convergence comparison. `infeasible` is a solver record, not evidence of biological death.

## What can and cannot be separated

Before the final constrained interval, all three trajectories use the same \(\mu=1.464804644951635\ \mathrm{h^{-1}}\) and glucose uptake 10 mmol gDW\(^{-1}\) h\(^{-1}\). Under the stored explicit-Euler rule,

\[
B(T^*)=B_0(1+\mu\Delta t)^n(1+f\mu\Delta t).
\]

For 0.0625 and 0.03125 h this reconstructs the saved common-window biomass within \(1.3\times10^{-14}\) gDW/L. Their 0.225522401210515-doubling gap is therefore Euler biomass approximation for this metric and window.

For 0.03125 to 0.015625 h, the additive decomposition is:

- explicit-Euler continuation: +0.117764542066295 doublings;
- quarter-step final glucose-limited/event correction: -0.008036581782056 doublings;
- observed net: +0.109727960284239 doublings.

Each track has exactly one reduced-growth final interval, but its cap and location depend on \(\Delta t\):

| dt (h) | final constrained interval (h) | glucose cap/uptake | growth rate (h⁻¹) |
|---:|---|---:|---:|
| 0.0625 | 5.25–5.3125 | 4.019587321722 | 0.563426506265 |
| 0.03125 | 5.15625–5.1875 | 2.045547764650 | 0.265781143412 |
| 0.015625 | 5.09375–5.109375 | 7.583797088062 | 1.101144555476 |

The terminal-time shifts (0.125 then 0.078125 h) cannot be uniquely assigned to one mechanism. The horizon-dependent inventory cap \(C/(B_n\Delta t)\), left-endpoint exposure, zero clamp, solver active set, and terminal infeasibility are coupled. Twelve other finite inventories are not emitted, and a dedicated raw pre-clamp glucose field was not stored. Reconstructed final raw glucose is only roundoff-sized (-4.44e-16, 0, and 0 mmol/L), but that does not prove glucose was the exclusive limiting constraint.

The separate reserve-depletion package confirms its own Q9 budget and zero-snap behavior under an R385-closed condition. It does not establish convergence or Q9 causality here.

## Isolated numerical candidate

`event_exposure_integrator.py` implements exactly one method, `constant_specific_rate_event_exposure_v1`. For prescribed constant specific growth \(\mu\) and withdrawals \(q_i\), it uses

\[
E(\tau)=
\begin{cases}
B_0\tau,&\mu=0,\\
B_0\operatorname{expm1}(\mu\tau)/\mu,&\mu>0,
\end{cases}
\]

for biomass and every inventory. It finds the earliest \(C_i/q_i\) exposure, advances all states to that same event, and requires a caller to re-solve afterward. It does not post-hoc clip solver fluxes.

The three passing examples are: smooth growth with nutrient and Q9 withdrawals plus an internal split-step identity; an H-Q9-1-like reserve event that doubles biomass from 0.01 to 0.02 gDW/L at 0.4731134451 h; and zero-growth maintenance that exhausts glucose at 0.25 h without changing biomass.

This is exact only for fixed rates over a segment. In a GEM, growth and uptake are endogenous, making the exact inventory constraint nonlinear. The candidate is therefore neither a drop-in FBA cap nor evidence that a real-GEM trajectory will converge.

## Exactly one future real-GEM comparison — not run

The only next numerical decision is whether to authorize **one event-and-resolve WT trajectory**; do not authorize a grid, KO, or production replacement yet.

Freeze the existing model, medium, profile, runner context, solver settings, `po1f_nonlimiting`, no added KO within the PO1f overlay, \(\alpha=10^{-4}\), pool multiplier 1, initial biomass 0.01 gDW/L, \(dt=0.0625\) h, and a 0–6 h horizon. Reuse the stored Euler 0.0625 and 0.015625 h trajectories as comparators. Record all 14 finite inventory starts, raw ends, saved ends, exposure, and event/re-solve log.

Budget: one new WT trajectory, at most 96 nominal intervals, at most 256 backend optimization calls, one CPU thread, 600 s wall time, no retry. Stop on any identity mismatch, unrecorded pool, substantive pre-clamp deficit, unexplained pre-depletion nonoptimal status, or budget overrun.

Acceptance is conjunctive:

1. initial solver status/growth/flux state matches the stored context before integration;
2. every biomass and inventory balance uses the same exposure, raw residuals stay within documented tolerance, and no substantive clamp occurs;
3. an inventory event is recorded and followed by a fresh solve rather than a hidden overshoot;
4. at the fixed common time \(T^*=5.109375\) h, the candidate 0.0625 h absolute doubling gap from the stored 0.015625 h diagnostic reference is strictly below the old 0.335250361495 gap; the candidate's first advanced interval with positive starting glucose and zero ending glucose must occur within 0.0625 h of the stored quarter-step glucose-zero diagnostic at 5.109375 h.

Passing would only justify a separately authorized refinement comparison; it would not establish convergence or permit production adoption. Failing rejects this candidate path without changing the model or scientific assumptions.

## Reproduction

From the repository root:

```bash
/usr/bin/python3 -B artifacts/dfba_consistent_integrator_n2_20260908/recompute_existing_diagnostics.py
/usr/bin/python3 -B artifacts/dfba_consistent_integrator_n2_20260908/event_exposure_integrator.py
```

Both scripts use only the Python standard library. The first verifies all three trajectory hashes before reading them and performs zero solver calls. The second runs the three deterministic examples. Exact input and environment identities are in `run_manifest.json`; file identities are in `SHA256SUMS`.

No named gene is a biological subject of this WT-only N2 analysis. “WT” is operationally defined above and must not be read as a claim that the overlay is an unmodified wild-type genome.
