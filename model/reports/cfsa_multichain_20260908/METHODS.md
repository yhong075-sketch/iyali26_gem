# Bounded multichain diagnostic methods

## Material Passport

Task: implement numerical diagnostics for the planned same-polytope CFSA pilot; no solver calls or model sampling in this implementation subtask. Inputs are arrays ordered **chains × retained draws × variables**. The function `diagnose(data)` returns `rhat`, `ess_bulk`, and `ess_tail` arrays of length variables, plus `lag1` with shape chains × variables. NumPy 2.4.3 and SciPy 1.18.0 were used for the implementation check. This document describes methods, not results from the planned CFSA chains.

## Definitions and primary implementation

For each variable, split each chain into equal first/last halves; an odd central draw is omitted. Pool average ranks across split chains and transform using z = Φ⁻¹[(rank − 3/8)/(S + 1/4)]. Compute R̂ = √[((N−1)W/N + B/N)/W], with N draws per split chain, W the mean sample variance within chains, and B = N × variance of chain means. Report the maximum of this rank-normalized R̂ and the corresponding R̂ after folding split values about their pooled median.

Bulk ESS applies the multichain ESS estimator to split rank-normalized values. Tail ESS is the minimum ESS of indicators x ≤ q₀.₀₅ and x ≤ q₀.₉₅, with linear quantiles calculated on the original pooled values before splitting. FFT autocovariances use the draw count as denominator. The estimator combines within- and between-chain variance, Geyer's initial positive and monotone paired sequences, the final positive even-lag correction, and the ArviZ lower bound 1/log₁₀(S) on integrated autocorrelation time. These definitions follow the pinned [ArviZ v0.22.0 diagnostic functions](https://github.com/arviz-devs/arviz/blob/v0.22.0/arviz/stats/diagnostics.py) and [autocovariance/quantile helpers](https://github.com/arviz-devs/arviz/blob/v0.22.0/arviz/stats/stats_utils.py).

Lag-1 is the Pearson correlation between adjacent retained values within each original chain, calculated before splitting. It is descriptive and has no independent pass threshold here.

## Frozen pilot criterion and limits

The coordinating task specifies three chains × 1,000 retained draws, seeds 20260921/20260922/20260923, the same feasible polytope and no extra burn-in removal. It requires **R̂ < 1.01, bulk ESS ≥ 400 and tail ESS ≥ 400 jointly** for every tested varying reaction and five key observables. Those are task criteria; this helper neither selects variables nor changes thresholds. Report individual and joint passing fractions, undefined counts, maximum R̂ and minimum bulk/tail ESS with identities. Medians alone are insufficient. Any nonfinite diagnostic is nonpassing.

The [Vehtari et al. primary paper, version 5](https://arxiv.org/html/1903.08008v5), recommends rank/fold diagnostics, R̂ < 1.01, sufficiently large ESS and at least four chains by default (§2); the planned three-chain run is a bounded pilot limitation. Six split chains do not constitute six independent initializations. Distinct seeds do not themselves establish dispersed starts. Satisfying marginal diagnostic thresholds cannot prove convergence, global exploration, biological accuracy or engineering-target validity, and bulk ESS is not a direct ESS for every downstream functional.

## Deliberate adapter behavior

- `_rhat` and `_ess` were extracted from the pinned ArviZ source, preserving the arithmetic and Geyer loops while replacing optional Numba variance dispatch with NumPy. The rank, split, autocovariance and quantile wrappers use installed NumPy/SciPy equivalents. No ArviZ/xarray installation was added.
- Globally constant or nonfinite variables return NaN throughout. ESS for a constant transformed/indicator input returns NaN rather than ArviZ's nominal draw count. `np.maximum`/`np.minimum` propagate undefined folded or tail components. These fail-closed departures prevent an undefined diagnostic being treated as a pass.
- Scientific fixed/near-constant classification belongs to the frozen coordinating analysis. Do not let a low within-chain range hide different constants across chains. Excluded items must retain their identities and reasons; exclusions are not passing diagnostics.
- The input boundary requires at least two chains and eight draws. Odd-chain splitting, thresholds on the original draws and lag-1 shape are explicit. No burn-in is silently discarded beyond the middle draw required for equal odd-length halves.

## Verification and provenance

Run `.venv/bin/python artifacts/cfsa_multichain_20260908/mcmc_diagnostics.py` from the project root. The small self-check passed for IID normal draws, a shifted chain, stationary-initialized AR(1) draws with ρ=0.95, constant inputs and a nonfinite input. An independent focused check matched 50 finite values exactly against extracted official ArviZ functions across normal/tied inputs and 8/9/100/101/1,000 draws (maximum absolute difference 0); the conservative NaN departures above were identified and documented. The scope and evidence are in `focused_code_review.md`. These are synthetic software checks, not new CFSA samples or a complete validation of MCMC estimators; no broad or biological audit is claimed.

Official release snapshots and the Apache-2.0 license are preserved under `vendor/`; `method_sources.json` records exact source URLs, version, fetch date and full SHA-256 values, plus the implementation identity. Network retrieval initially failed inside the restricted sandbox, then the authorized scoped download succeeded. The saved release is selected for reproducibility, not asserted to be the latest release. The paper's identity and relevant methods/recommendations were checked through arXiv; implementation details were checked against the official source.
