# R305-only correction on the CoQ9 growth-coupled candidate

Executed 2026-10-05 America/Los_Angeles (started 2026-10-06 01:07 UTC). **The requested correction was applied, exported and reloaded. R305 balances exactly; CoQ-dependent growth remains intact; both specified closed-input ATP tests return optimal zero. Exactly four new candidate optimizations were performed.**

## Model and edit

The [pre-edit control](../coq_biomass_candidate_20261005/E5_coq9_biomass_alpha_1e-4_validated.xml) was preserved byte-for-byte. Before editing, its biomass requirement was confirmed as `biomass_C: m468[C_mi] = −0.0001`, with total mature CoQ pool row `R385 − 0.0001 × biomass_C = 0`.

The new, separate [R305 candidate XML](E5_coq9_alpha_1e-4_R305_qcycle.xml) contains these final signed coefficients:

| R305 metabolite | Control | Candidate |
|---|---:|---:|
| m28[C_mi], H⁺ | −1.5 | **−2** |
| m10[C_cy], H⁺ | +1.5 | **+4** |

Reloaded equation:

\[
Q_9H_2+2\,cyt\,c_{ox}+2H^+_{mi}
\longrightarrow Q_9+2\,cyt\,c_{red}+4H^+_{cy}.
\]

In model IDs:

```text
2 m28[C_mi] + 2 m2[C_mi] + m471[C_mi]
  -> 4 m10[C_cy] + 2 m3[C_mi] + m468[C_mi]
```

The [stoichiometric diff](stoichiometric_diff.tsv) contains exactly those two entries. The only additional XML change is R305's `coq9_proton_state` note, now recording the applied candidate and its limitations. Full parsed-XML comparison and reloaded model/notes comparison verified preservation of all other coefficients, GPRs, bounds, formulas, charges, annotations and compartments. The α value remains an uncalibrated assumption; `C_cy` remains cytoplasm under the existing compartment approximation, not an experimentally established intermembrane-space compartment.

Before any optimization, exact arithmetic with the stored formulas/charges found **C, H, Fe, N, O and S residuals all zero; charge residual zero**. The previous R305 H and charge residuals were both −2.

## Four candidate results and matched controls

The comparison is the same CoQ-coupled model with old versus revised R305. Control results below are **reused existing results**, not new control solves. Model SHA, medium, strain overlay, solver settings and each complete reconstructed optimization problem—including its constraint matrix—matched the four saved `candidate_*` problems from the [previous checks](../coq_candidate_checks_20261005/REPORT.md) exactly.

| Quantity | Pre-edit control, reused | R305 candidate, new |
|---|---:|---:|
| Normal maximum growth, h⁻¹ | 1.4289162961 | **1.4445584333** |
| R305 flux at normal growth | 28.2360050951 | **27.1788473760** |
| R385 flux at normal growth | 0.0001428916296 | **0.0001444558433** |
| R385 − 10⁻⁴ × growth | 0 | **0** |
| Maximum growth with R385 disabled, h⁻¹ | 0 | **0** |
| Closed-input cytosolic maximum ATP dissipation | 0 | **0** |
| Closed-input mitochondrial maximum ATP dissipation | 0 | **0** |

All four new solver statuses are **optimal**. Reaction fluxes and ATP dissipation are in mmol/gDW/h. Zero values are solved optima, not infeasible cases, failures or substituted missing values. The normal-growth optimum increases by **1.094685%** under this isolated edit. R305 fluxes are values in the saved optimal solutions; no FVA or flux-uniqueness claim is made.

Both growth tests retain the original static PO1f/SD-Leu conditions and NGAM lower bound 7.8625. R385 disabled means its bounds are [0,0]. The candidate's zero-growth solution remains feasible with maintenance retained; its saved R305 flux is 4.0733327268, so zero growth is not a statement that every metabolic flux is zero.

The ATP protocol is unchanged: close all 183 single-sided reactions, close the six configured biomass reactions, and relax only the compulsory NGAM lower bound to zero. Maximize the existing balanced cytosolic ATP drain or the same temporary balanced mitochondrial ATP + H₂O → ADP + phosphate drain. All closure edits match the control protocol; zero flux is feasible. No additional constraints were relaxed and no follow-up optimizations were run. ATP results apply to these two compartments, drains and fully closed boundary at tolerance 10⁻⁷; they do not establish global thermodynamic or native biological validity.

## Evidence and scope

The [curation specification](../../data/coq_r305_candidate.json), [executed script](run.py), [manifest](manifest.json), [results table](results.tsv) and [four-call ledger](solver_budget.json) record the source/candidate identities, full code hashes, medium/strain context, raw statuses, equations and budget. Per-case JSON files retain the full optimization problems and flux witnesses. Original inputs and protected files were unchanged at completion. Historical dirty state is recorded, not reconstructed.

A bounded independent read-only audit checked the complete XML pair, exact chemical residuals, saved old/new problems, numerical witnesses, closure and solver ledger: **6 claims | 6 audited | 6 supported | 0 unresolved | 0 contradicted | 0 unchecked**. It confirmed that the only old/new optimization differences are the intended two R305 entries and corresponding matrix coefficients. The audit ran no optimizations or literature searches.

This delivers the requested isolated model candidate and tests. Default builds and earlier models remain unchanged; no new GPR, mechanism acceptance, commit or push is implied.
