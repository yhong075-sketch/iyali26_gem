# N1 — CFSA sampling candidate evaluation

## Decision

**Do not continue this candidate in its present form.** The shared warmup-pair hit-and-run pilot improved several diagnostics and reduced solver setup work, but it did not meet the preregistered criterion and made two global observables worse. These results do not establish better sampling from the target distribution.

## Frozen problem and target

The executed input was `scenario_2.json.gz` (SHA-256 `3296863c5fa0be686098e55de78aaf8dffd660b869ce89e716ba518ec3dab61d`). The sampled split-variable region used:

- growth ≥ 0.14650760568106092;
- mixed-lipid pseudopool demand ≥ 23.74006264072655;
- \(\sum_i(f_i+r_i)\le 2330.701999560738\);
- the same remaining reaction bounds, medium, stoichiometric equalities, and solver tolerances as the preserved C baseline.

The complete sampler matrices and variable order matched the C baseline exactly. The cap was restored from the preserved executed solver problem, rather than inferred from a summary JSON. Net flux was calculated only as \(v_i=f_i-r_i\). The nominal target remained uniform on the complete **split-variable** feasible region; no uniformity claim is made for its projected net-flux distribution.

## Candidate and target-preservation argument

At each transition, the candidate selects an ordered pair of distinct frozen FVA warmup points and uses

\[
d=w_i-w_j.
\]

Both points satisfy the same equalities, so \(Ad=0\); the observed maximum equality residual for a checked warmup direction was \(4.26\times10^{-14}\). Uniform ordered-pair selection gives \(q(d)=q(-d)\), and the direction law is independent of the current state. Uniform sampling along the feasible chord therefore retains the uniform split-space distribution as the intended invariant target. All hard reactions and constraints remain present.

COBRApy invoked its numerical fallback 6 times in 300,000 transitions. That fallback resets from the sampler center and is outside the clean detailed-balance argument. It was rare, but it limits the claim to preservation of the **nominal** target; this run does not prove exact stationarity or exact uniformity.

## Fair comparison

Both methods used 3 independent seeded chains, 1,000 retained draws per chain, thinning 100, and therefore 300,000 transitions and 3,000 retained draws. Bound-fixed and observed near-constant reactions were excluded from convergence pass counts.

| Metric | Native OptGP baseline | Candidate |
|---|---:|---:|
| Feasible retained draws | 3,000 / 3,000 | 3,000 / 3,000 |
| Bound-fixed reactions | 69 | 69 |
| Observed near-constant reactions | 876 | 876 |
| Varying reactions evaluated | 1,369 | 1,369 |
| R-hat < 1.01 | 64 | 188 |
| Bulk ESS ≥ 400 | 0 | 14 |
| Tail ESS ≥ 400 | 22 | 130 |
| All three criteria | 0 | 13 |
| Backend optimizations | 13,470 | 4,490 |
| Hit-and-run transitions | 300,000 | 300,000 |

The backend-call reduction comes from generating one shared frozen warmup set instead of regenerating the same-domain setup for each chain. It is a setup-cost result, separate from mixing quality. Historical wall times (126.71 s baseline; 58.82 s candidate through sampling and the failed JSON write) were not collected as a matched timing experiment, so no speed multiplier is claimed.

| Key observable | Baseline R-hat / bulk ESS / tail ESS | Candidate R-hat / bulk ESS / tail ESS | Candidate pass |
|---|---:|---:|:---:|
| Growth | 1.0039 / 305.0 / 204.1 | 1.0105 / 303.6 / 446.5 | No |
| Mixed-lipid demand | 1.0113 / 180.2 / 139.5 | 1.0044 / 496.4 / 868.9 | Yes |
| Glucose uptake | 1.0254 / 92.1 / 95.5 | 1.0925 / 30.8 / 18.7 | No |
| Net L1 | 1.0664 / 35.7 / 95.6 | 1.1565 / 12.8 / 42.2 | No |
| Forward + reverse total | 1.0093 / 221.2 / 193.9 | 1.0013 / 533.2 / 418.4 | Yes |

The candidate passed all criteria for 13/1,369 varying reactions and 2/5 key observables. The preregistered overall criterion therefore failed. The mixed-lipid and split-total improvements are real observations from this bounded pilot, while the poorer glucose-uptake and net-L1 diagnostics rule out a general improvement claim.

## Execution record and limitation

The one authorized run used 4,490 backend optimizations and 58.82 seconds, within the limits of 16,000 optimizations and 1,200 seconds. Sampling completed before a NumPy integer failed JSON serialization in the combined runner. The exact executed code was preserved as `run_candidate_executed.py` with SHA-256 `70753beda64f2120c71491e2b3d69761b2c7f79b8cc59278104a83d4313c1405`; it matches the code hash recorded during execution. The three already-written feasible chains were then diagnosed offline with zero additional optimization or sampling. The runnable source fixes only that serialization cast.

This is a computational diagnostic result for one constrained CFSA region. It is not a TAG-yield result, a physiological validation, or evidence for gene or model changes.

## Only next decision

Close this candidate as a negative/mixed pilot. Any further N1 work should preregister a different single kernel or geometry treatment rather than lengthening these chains.
