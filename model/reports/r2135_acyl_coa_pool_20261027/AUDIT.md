# Independent audit: R2135 acyl-CoA pool remap and R2145 1:1 rescale (2026-10-27)

**Auditor:** a fresh read-only general-purpose agent with no conversation history. It was given the
deliverables, the logs and the claims, and was not given the reasoning.

**Scope:** this folder, the curation file, the step, the test and the `--r2135-acyl-coa-pool` lines.

**Method:** it recomputed every claim itself, using shasum, a raw XML diff, external sources, and Gurobi 1e-9
LPs on in-memory models. It also ran the focused pytest. It changed no file.

**Overall verdict:** "the change is what was authorized and nothing more". No blocking defect was found.

| Claim | Verdict | Auditor's own evidence (summary) |
| --- | --- | --- |
| C1 SHA-256 of baseline / candidate | supported | `00e8e397…`, `e9d35857…` |
| C2 diff limited to R2135 / R2145 stoichiometry + notes | supported | Raw XML diff: m68 → m567 in R2135, 0.01 → 1 on m1993 / m1994 in R2145, two new notes blocks; nothing else at byte level. |
| C3 control = baseline; seeds identical | partially supported | Bytes verified for all builds, and the final-round `.build.json` records the final code and data SHAs. The **PYTHONHASHSEED values are not recorded** in `.build.json` or the logs; only the file names carry them. |
| C4 identities (MNXM1136895, acetyl-CoA, rxn09551) | supported | Opened MetaNetX, KEGG REST (C00040, C00024, R02251, R01496, yli:2911239, yli:2910950), UniProt Q6C7L9, PMC9505966 and the ModelSEED Solr record. Additional facts: ModelSEED lists rxn09551 as alias of BiGG DAGPYP_SC (iMM904 / iND750) and flags it mass-imbalanced; R2145's MNXR144651 differs from ModelSEED's alias MNXR136076. |
| C5 WT 0.3120049 → 0.3116824 | supported | Recomputed; Pi exchange −0.02589 → −0.08455. |
| C6 no-uptake ATP and biomass 0 | supported | All 185 boundary reactions closed to import, xMAINTENANCE relaxed: 0.0, optimal, both models. |
| C7 FAS knockout lethal on the candidate with and without R1726, not lethal on the baseline with R1726 | supported | 0.0 / 0.0 / 0.312056; knocked flux 0. |
| C8 no class change in the two-model sweep | supported | `ko_compare.tsv` re-read, 0 changes at all cutoffs; spot checks YALI1D10264g, YALI1D35593g, YALI1F01981g reproduce. |
| C9 no class change in the stacked R1726 sweep | supported | Stacked TSV vs candidate columns: 0 differences, 1,081 genes. |
| C10 re-census, R1343 / R1344 GLPK fringe | supported | Gurobi 1e-9 FVA on both models: R1343, R1344, R970 max 1.47e-4. **R970 is still GLPK-blocked on the candidate**, as on the baseline. |
| C11 pytest equals before + 4 new tests | supported (logs) | `pytest_before.log` predates the code edits; `pytest_after.log` postdates the final code. The full suite was not re-run. |
| C12 fails closed, idempotent | supported | Code read; focused test 4 passed, 29 subtests. Low observations: `before_stoichiometry` is pinned only in data; the cli.py diff also edits the cofactor help text. |
| C13 limits stated honestly | partially supported | See the gaps below. |

**Gaps the auditor raised and what was done:**
1. AUDIT.md was linked but missing: this file.
2. The seeds are not recorded. REPORT now says so. The builds were run as `PYTHONHASHSEED=0 …` and
   `PYTHONHASHSEED=1 …` on the command line, but the build record does not store the variable.
3. "At the candidate's optimum, TAG comes through R1771" holds at the pFBA point only. A plain FBA optimum
   used R2135 → R2168 at 7.3e-5 instead (alternate optima). REPORT reworded.
4. "The 0.01 comes from rxn09551" is provenance by annotation (inferred). REPORT now says so, and mentions
   the BiGG DAGPYP_SC alias and ModelSEED's own imbalance flag.
5. The code review was not on file: [logs/code_review.md](logs/code_review.md) now records its findings and
   fixes.
6. Positive coverage and the unresolved exclusion were not stated. REPORT now states them: 326 in-model
   positives, and 16 unresolved excluded from the confusion counts.
7. The budget overrun was not disclosed: 6 builds instead of 3, because the builds were repeated after the
   review fixes, and about 70 built-file LPs instead of about 20. REPORT now discloses it.
8. Criterion 9.6 relied on an earlier Gurobi fringe check; the auditor's own check on `e9d35857` closes it
   (C10).

The auditor's scratch scripts are in the session scratchpad. They are not part of the deliverables.
