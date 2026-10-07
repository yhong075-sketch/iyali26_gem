# CoQ9 candidate: chemistry, growth dependency and closed-input ATP checks

Checked on 2026-10-05 (America/Los_Angeles; run started 2026-10-06 00:56 UTC). **CoQ-dependent growth works, and the specified closed-input ATP tests pass; R305 local chemistry fails.** These are eight new optimizations of the saved XML files, not reused growth results.

The comparison uses [original E5](../atp_candidate_repair_20260924/candidates/E5.xml) and the [selected α = 10⁻⁴ candidate](../coq_biomass_candidate_20261005/E5_coq9_biomass_alpha_1e-4_validated.xml). The candidate adds 0.0001 mmol CoQ9/gDW to `biomass_C`, an uncalibrated assumption. It retains E5's R305 chemistry. The filename “validated” refers to its construction/export checks, not acceptance of all model chemistry.

## 1. Local chemistry — fails

Both files encode the same R305 equation:

\[
Q_9H_2+2\,cyt\,c_{ox}+1.5H^+_{mi}
\rightarrow Q_9+2\,cyt\,c_{red}+1.5H^+_{cy}.
\]

Exact arithmetic using the stored molecular formulas and charges gives products-minus-reactants residuals of **−2 H atoms and −2 charge units**; all other elemental residuals are zero. This is a stoichiometric defect, not solver noise.

The existing [curation record](../../data/coq9_curation.json) contains a separate Q-cycle proposal: mitochondrial H⁺ coefficient **−2**, cytosolic H⁺ coefficient **+4**. Substituting those two coefficients gives zero elemental and charge residuals. This substitution was checked arithmetically only; it was **not applied to either model or used in the optimizations**. The present reaction note explicitly says the Q-cycle candidate was not applied.

Balance alone requires the product-side proton coefficient to exceed the reactant-side coefficient by two; it does not uniquely establish 2/4 or validate a mechanism. The proposal also uses `C_cy` as a P-side/intermembrane proxy. Its existing mechanism rationale is in the [prior CoQ report, §3.1](../coq9_pipeline_20260909/input/REPORT_zh.md); no new literature verification was performed here. Thus the current candidate cannot pass the requested mechanism/chemistry criterion.

## 2. Growth and CoQ dependency — passes for the current model

Static PO1f/SD-Leu conditions, maximum `biomass_C`, glucose uptake limit 10 and the existing nonlimiting uracil allowance were used. NGAM remained 7.8625 mmol/gDW/h throughout these growth tests. These conditions differ from the earlier finite-uracil dFBA runs.

| Model / perturbation | Maximum growth, h⁻¹ | R385, mmol/gDW/h | Solver status |
|---|---:|---:|---|
| E5 baseline | 1.429177892279185 | 0 | optimal |
| E5, R385 disabled | 1.429177892279190 | 0 | optimal |
| α = 10⁻⁴ candidate | 1.428916296065595 | 0.00014289162960656 | optimal |
| Candidate, R385 disabled | **0** | **0** | **optimal** |

Summing the complete Q9 and Q9H₂ metabolite rows identifies only R385 as net input and candidate biomass as net removal:

\[
v_{R385}-10^{-4}\mu=0.
\]

The saved candidate solution has zero residual for this relation. Setting R385 to [0,0] therefore forces μ = 0, confirmed by an optimal, feasible zero-growth solution with normal maintenance retained. This is reaction dependency under the encoded assumptions; it neither establishes gene essentiality nor validates the native CoQ pathway. CoQ9 is a required biomass constituent here, not the sole growth substrate.

## 3. Closed-input ATP consistency — passes within the tested scope

For each model, all single-sided sources, sinks and exchanges were closed, including oxygen, water and proton exchanges; all six configured biomass reactions were closed. NGAM's lower bound was relaxed from 7.8625 to 0. Zero flux satisfies the resulting variable and constraint bounds, so a forced-maintenance infeasibility cannot explain the zero ATP result.

Each test maximized a balanced ATP hydrolysis drain, which measures the network's ability to regenerate ATP continuously. Cytosolic ATP used the existing maintenance reaction; mitochondrial ATP used a temporary balanced ATP + H₂O → ADP + phosphate drain in the model's stored neutral chemical representation.

| Maximum ATP dissipation, mmol/gDW/h | E5 baseline | α = 10⁻⁴ candidate |
|---|---:|---:|
| Cytosol | 0, optimal | 0, optimal |
| Mitochondrion | 0, optimal | 0, optimal |

No positive substrate-free ATP generation was detected at the 10⁻⁷ tolerance. No anomaly-directed follow-up solves were needed. These results apply to these drains and the stated fully closed boundary; they do not prove global thermodynamic consistency or cure R305's local imbalance.

## Execution record and disposition

Eight optimizations completed; no model, GPR, medium file, default build or existing result was modified. No new XML was exported. Input and code hashes remained unchanged during execution. Full identities, actual media/strain overlays, solver settings, chemistry and closure edits are in [manifest.json](manifest.json); raw results are in [growth.tsv](growth.tsv), [energy.tsv](energy.tsv) and their per-case JSON files. The [solver ledger](solver_budget.json) and [archived check script](run.py) retain the bounded execution record. The script uses original manifest paths and refuses overwriting this run; it is not a portable rerun interface. The historical dirty environment was recorded, not reconstructed. The prior CoQ mechanism report linked above is a local historical reference; the tracked curation record retains the proposed coefficients.

One independent read-only check covered the run code, recorded hashes, the two result tables and all eight saved optimization problems/results. It confirmed the pool identity, R385 knockout outcome, all 183 single-sided closures per model, zero-feasible bounds, the four ATP maxima and the chemistry failure. No additional optimization, literature search or separate audit artifact was produced.

**Disposition:** the biomass coupling is operational; this candidate is not an R305 chemistry repair. Any later mechanism correction needs its own candidate and checks using the changed equation.
