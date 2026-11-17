# Independent audit: methylcitrate cycle on 7ccbb568, 2026-11-16

## How the audit was run

- **Auditor.** One fresh, read-only general-purpose agent, given no history and no expected verdicts. Its brief named
  TASK.md, REPORT.md (first draft, pytest section still pending), the data, code, logs, outputs, both model files,
  the large outputs and the other worktree the step was copied from.
- **Rules.** No writes, no build, no LP. Reading SBML with cobra, hashing, diffing and running the step in memory
  were allowed.
- **What it recomputed:**
  - every SHA, and byte identity with `cmp`: control vs 7ccbb568, seed 1 vs seed 0, the copied step and data vs the
    other worktree, the coq9 TSVs;
  - the two `.build.json` files, key by key;
  - its own whole-model cobra diff and a plain XML diff, which is byte-identical to
    [logs/xml_line_diff_vs_7ccbb568.txt](logs/xml_line_diff_vs_7ccbb568.txt);
  - the LP identity of 353cefdc and 7ccbb568 (bounds, stoichiometry, GPR, objective);
  - the matrix CoA-moiety scan on both models;
  - the step in memory on 7ccbb568 compared with the built file, and a rerun on the built file;
  - four extra guard probes;
  - the new test (8 passed, 43 subtests, 70 s);
  - every number in REPORT §4–§5 from the JSON and TSV outputs, including its own full-column sweep comparison and
    the scores;
  - the figure texts.
- **Web.** 12 read-only requests: KEGG REST, the Rhea TSV, Europe PMC (abstracts and two full texts) and UniProt
  counts.
- **Scope.** The full pytest suites were still running, so test criterion 7 was not audited. The user's
  authorizations and the cbf2ef95 state of the other line were outside its reach.

## Supported

48 claims were checked:
- 43 were supported, some with the qualifications noted below;
- 3 were partially supported (claims 25, 34, 35);
- 1 was overstated (claim 33);
- 1 was unsupported (claim 36).

The 5 that fell short are P1–P5 below. The supported ones include:
- every SHA and byte-identity claim;
- the argv: the baseline's 47 arguments plus the flag;
- the build record;
- that only the targets changed;
- R_OBDHm's definition and balance (C46H69N14O33P5S on both sides, charge 0);
- R2075's state;
- the guard (fails closed, all or nothing, rerun "already_correct");
- the wiring: opt-in, last;
- every WT, FVA, energy-check and knockout number;
- the pFBA equal-total-flux claim;
- the forced route (no CoA or acyl-CoA transport);
- the sweep: same gene set, no status, certificate or class change, maximum |ΔKO| 1.39e-14;
- the scores and the 326 positives;
- the figure numbers (spot check);
- the after-copy containing the code changes.

External sources:
- **KEGG** (no yli or sce gene for EC 5.4.99.2; `hsa:4594`), **PMID 11092862** (ICL2 in the matrix), **PMID 3718468**
  (rat PDH and BCKDH), **PMC5991449** (phd1 deletion, no growth on propionate), **Rhea 33223** and the **UniProt**
  counts in R2075's note were all read and match.
- **PMC6987463** matches as a general pathway statement in its introduction, not a *Y. lipolytica* experiment.

## Problems raised and how they were resolved

| # | Problem | Resolution |
| --- | --- | --- |
| P1 | Capacity changes outside the nine were not reported. The WT maximum maintenance flux rose 63.33341 → 63.37226. R693, R412, R741, R83 and R410 maxima rose slightly. "Nothing else in the model changes" was too strong | Verified (sweep summaries: Δ +0.0388). Added to the answer, §4 (table and paragraph) and §9. The wording is now "no call, score or WT growth changes". Energy check 0; the NADH explanation is marked inferred |
| P2 | `YALI1C24124g` was called listed essential; R25's genes are not listed either | Verified: in the sweep TSV, only `YALI1F03803g`, `YALI1F39620g` and `YALI1F20788g` have `experimental_positive=True`. §8 is split into listed (FN, KO/WT 1.0) and not listed (TN), and §9 is reworded |
| P3 | §3 called "no reaction could take CoA moiety out" proven; "34" is the 7ccbb568 count | §3 now separates the proven identity from the observed (LP) and statically verified parts, gives 34 in 7ccbb568 and 35 in ad901ddb (audit re-scan), and says the scan was not re-run in this round's scripts |
| P4 | R963's 7ccbb568 value (2.67) was the forward maximum, while the figure prints 2.94 (largest absolute) | The §4 table now gives min / max for every row |
| P5 | "R834 only through side route B" was stated as fact | Reworded in the answer and §3: the observed route is B, and that it is the only route is inferred (no R834 route solve on ad901ddb) |
| P6 | The status-figure legend did not explain the blue-grey arrows | Legend now reads "blue-grey = can carry flux before and after; light grey = blocked"; figures re-rendered |
| P7 | A session anchor link does not resolve in the repo | Replaced by the session title and its branch |
| P8 | Pytest and AUDIT.md were pending | This file; pytest is in REPORT §6 |
| P9 | `logs/pytest_after_extra_files.txt` lists its own temporary `.z` name, and the after-copy was made at 19:22 | Recorded in REPORT §6. Harmless: all code files in the copy equal the worktree |
| P10 | The LP total excludes LPs run inside the pytest suites | Stated in REPORT §11 |
| P11 | "43 subtests of rejection" (it is 38 rejection + 5 knockout); "statically verified" for tests; relative `--output-model` in §1 | Corrected: "Tested"; the subtest split; the absolute output path |

## Not checked by the audit

- The user's authorizations (TASK §0).
- That cbf2ef95 shares the pinned state. Only the step's success on 7ccbb568 was checked.
- The full pytest comparison.
- Whether the investigation's in-memory LP results hold beyond what this round re-solved.
- The biology of R2112's direction and of side route B.
