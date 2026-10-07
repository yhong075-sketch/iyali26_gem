# Independent numerical audit

Audited 86 saved solver outcomes, including 74 optimal primals, 5 full route witnesses and 6785 ledger entries. All checks passed. No optimization was run.

Maximum native constraint violation: 9.10772e-09. Independent full net S*v residual (compensated summation): 9.10772e-09. Maximum net reaction-bound violation: 4.44402e-08; maximum native split-variable bound violation: 5.21167e-08. Maximum net/split difference: 0. Native sparse coefficients reproduce the stored net stoichiometry and the sole permitted growth-floor constraint. Source and hydrolysis are fixed to the same positive q, paired routes share the same growth floor, and complete Gly-Pro pool accounting has no unexplained source or export.

All ledger terms and weighted route sums were independently recomputed. Every adjacent reaction for each ledger metabolite is represented. Closed energy tests retain all mass-balance rows and have no diagnostic source or growth floor; verdicts agree with actual statuses and objectives. The saved code and explicit input identities match their recorded hashes. A deliberately corrupted primal was rejected.

An independent ElementTree comparison of all four exported XML files confirms that original species, gene products, groups, objectives and metadata are unchanged, all original reactions retain their full XML definition apart from the declared R2039 capacity/isolation bound, and the only new reaction is the balanced, GPR-free cytosolic hypothesis. Historical chemical and energy notes are preserved with the original XML elements.

Primary FBA and growth-optimal pFBA differ in 96, 98, 90, 104 reaction net fluxes above 1e-7, for O_CY/C_ONLY, O_CY/V_ONLY, O_VA/C_ONLY and O_VA/V_ONLY respectively. Their recorded biomass values match. Whole-cell optimal fluxes are therefore not unique within the recorded numerical precision; target-pump FVA at matched q/growth is a separate result. Detailed differences and source hashes are in optimal_flux_multiplicity.json.

Scope: this validates the recorded optimization problems and primal evidence. It does not independently establish optimality or infeasibility through dual/Farkas certificates, and it does not verify native protein activity, localization, pH, membrane potential, or gene essentiality.
