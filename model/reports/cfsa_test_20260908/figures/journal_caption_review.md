# Proposed journal-caption review

Reviewed 2026-09-08T07:43:11Z. Scope: consistency of the proposed three-panel figure and caption with saved `sample_diagnostics.json`, `NUMERICAL_AUDIT.md` and `FIGURE_CAPTIONS.md`; applicable instructions and the relevant CFSA project-state entry were read. The review used existing verified evidence, with no solver calls, resampling, raw-sample recomputation, external source audit or source changes. The proposed final figure was not visually inspected. This is manuscript wording for human review, not publication acceptance or biological validation.

## Proposed legend

Comparative flux sampling of the provisional iYali26 reference under SD-Leu/PO1f conditions. Growing, slow-growing and producing scenarios each contain 500 retained samples from one chain (thinning interval, 100; no additional burn-in removal). All scenarios retain a minimum growth of 10% of the original wild-type maximum; the growing and producing scenarios use 90% optimality constraints. **A,** Growth rate versus mixed-lipid pseudopool demand flux for all 1,500 samples; diamonds denote arithmetic means and the horizontal axis is logarithmic. **B,** Sum of absolute net reaction fluxes per sample. Boxes show the interquartile range, central lines the medians and whiskers the full range; the vertical axis is logarithmic. Dashed segments mark the scenario-specific pFBA forward-plus-reverse flux caps, which also bound the net-flux sums; no such cap applies to the producing scenario. **C,** Empirical cumulative distributions of lag-1 Pearson autocorrelation across 1,369 varying reactions per scenario (sample peak-to-peak range >10⁻⁷). The lipid-demand constraint imposes the high-demand producing state. Lipid demand and total flux use model reaction scales and do not represent a TAG yield, measured secretion or enzyme cost. Saved net-flux samples passed the previously performed feasibility checks at a tolerance of 10⁻⁷; original forward/reverse trajectories were not saved. The draws are strongly correlated, and adequate mixing, convergence and biological target validity remain unestablished.

## Substantive labeling risks

- Label panel A's horizontal axis **Mixed-lipid demand (model flux scale)**, not TAG yield or productivity. Use **Growth rate (h⁻¹)** vertically. The scenario names identify imposed constraints, not observed phenotypes or successful interventions.
- Label panel B **Total absolute net flux (model flux scale)** or **Σ|vᵢ| (model flux scale)**. Its dashed caps originate in the forward-plus-reverse formulation; the saved net sum is not the original split-variable sum, and it is not an enzyme-cost measurement. The absence of a producing-scenario cap prevents interpreting the difference as a controlled biological effect; no specific futile cycle was identified by the saved audit.
- Label panel C **Lag-1 autocorrelation** and **Cumulative fraction of varying reactions**. Its denominator is 1,369 reactions within each scenario, not 500 independent samples. If a 0.9 reference line is retained, add “The 0.9 reference is descriptive” to the legend; it is not a convergence cutoff. High autocorrelation alone does not establish nonstationarity.
- The legend's full-range whiskers and median lines must match the actual rendering. Mean diamonds are specified for panel A; if added to panel B, append “diamonds denote arithmetic means” to its description. Do not describe the 500 draws as replicates.

## Evidence scope and identity

This review does not create a new source-audit or numerical-validation claim. Existing coverage remains: source review 6 total / 6 audited / 5 supported / 1 unresolved; numerical review 7 total / 7 audited / 6 supported / 1 unresolved. Both previously reported 0 contradicted and 0 unchecked. No new acceptance decision follows from caption preparation.

| Saved input | SHA-256 |
|---|---|
| `sample_diagnostics.json` | `28253f66d58fb4021cf3540fb1346fe93e9abd2f10b2819873852010e7f28472` |
| `NUMERICAL_AUDIT.md` | `7327f5edf467997c04545f94ccecf808c8768cfa7690b72ed881791326386705` |
| `figures/FIGURE_CAPTIONS.md` | `5076c21b6b2fe3661be7a335b408c63c5eddf54f5a435078f0f7f6e9c1db414e` |
