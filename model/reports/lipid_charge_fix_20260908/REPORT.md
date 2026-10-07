# Linoleoyl charge correction integrated into the lipid-unlump builder

The default `build_candidate` path now corrects all five linoleoyl-CoA copies before expanding lipid states. It retains the existing neutral CoA convention: C39H66N7O17P3S, charge 0, verified directly against https://www.ebi.ac.uk/chebi/CHEBI:15530 (formula and net-charge fields). Conflicting stearoyl identity/structure annotations are replaced with verified neutral linoleoyl identifiers. Original annotations and chemistry remain preserved in the frozen source and previous candidate.

The user explicitly authorized this bug fix and pipeline integration. Scope is the strict sn lipid builder and its curated input, regression tests, and a new candidate artifact. This is not the separate global biochemical-pH CoA migration. No formal model, old result, GPR, culture condition, or activation gate was changed. No commit or push was performed.

## Validation performed

- 14 regression tests passed in 37.484 seconds, including the existing growth, pFBA, route dependency, source preservation, and deterministic SBML round-trip tests.
- Independent formula oracle checked every TAG against glycerol + three fatty acids - three water molecules. All 343 TAG species have the expected formula and charge 0, including 108 single-linoleoyl, 18 double-linoleoyl and 1 triple-linoleoyl combinations.
- All 1,134 generated reactions are mass/charge balanced after SBML export and read-back.
- New and old candidate have identical reaction identities, stoichiometric coefficients, bounds, objectives and GPRs. All metabolite formulas are unchanged. This metadata correction therefore leaves the FBA optimization problem unchanged; it does not establish a new yield or resolve the earlier CFSA warmup failure.
- 341 metabolites have charge and/or annotation changes: five CoA copies and 336 propagated lipid states. The original 127 incorrectly charged TAGs are included in these 336 states.
- Among affected pre-existing reactions, R1870/R1871/R1872 now balance. R1869/R1873 retain pre-existing ±1 H/charge residuals, and the generic pool reaction R2165 retains elemental imbalance. The biomass reaction cannot be chemically assessed completely because m401[C_cy] has no usable elemental formula. No previously balanced, checkable affected reaction became unbalanced. These remaining chemical issues are not silently repaired by this scoped patch.

## Pipeline integration

The correction is in `scripts/lp_sn12_candidate.py`, backed by schema version 5 of `data/lipid_unlump_sn_core_curation.json`. The source fingerprint is still checked before mutation. The old legacy tuple remains recorded; the separate correction target is explicit. The builder rejects altered correction contracts and non-neutral generated TAGs. Tests are in `tests/test_lp_sn12_candidate.py`.

The actual source used was `../iyali26_gem_integration/model.xml`, whose SHA matches the builder's frozen source contract. The current directory's `model.xml` is a different input and is deliberately not substituted. The generated deliverable is `model_lipid_unlump_charge_fixed.xml` in this report's directory. `verification.json` records full source/code/output hashes, changes and residuals; `build_report.json` retains the existing broader activation status. `verify.py` repeats the saved old/new artifact checks without optimization or sampling. Chemical source verification was performed directly; no independent external reviewer was used.

Executed validation commands (repository root):

```sh
IYALI26_SOURCE_MODEL='/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_integration/model.xml' .venv/bin/python -m unittest discover -s tests -p test_lp_sn12_candidate.py
.venv/bin/python scripts/lp_sn12_candidate.py '../iyali26_gem_integration/model.xml' --report artifacts/lipid_charge_fix_20260908/build_report.json --candidate-sbml artifacts/lipid_charge_fix_20260908/model_lipid_unlump_charge_fixed.xml
.venv/bin/python artifacts/lipid_charge_fix_20260908/verify.py
```

The validation scope was the existing lipid regression suite, one build/report execution and read-only artifact comparisons; no CFSA resampling, global CoA migration or additional research matrix was run.
