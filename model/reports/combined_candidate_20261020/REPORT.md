# Combined candidate of the 2026-10-19 opt-in changes

**Candidate:** `model/candidates/E5_coq9_resp_pkg_r132_combined_20261020.xml`
(SHA-256 `164f235f4fe823181df0cc83082d6058ee9cf59675777adfe93557e9a77583a6`).
It is a56f3649 plus all eight user-approved opt-in changes from 2026-10-19, built in one run. It is
the newest candidate, not a designated reference model. Branch: `claude/combined-candidate-20261019`.
Scope, branches and acceptance criteria: [TASK.md](TASK.md).

## 1. What the candidate contains

| Change | Reactions | Effect |
| --- | --- | --- |
| 4-GPR restore | R612, R2274, R1302, R1303 | URA3, LIP2 and YALI1A08512g GPRs back. R2274 gets its 3 H₂O and R1303 loses the extra H⁺. |
| R671 restore | R671 | GatB YALI1E11370g off prephenate dehydrogenase again |
| EG-GPR-001 (CPA/URA2) | R159, R190 | R159 URA2-only and channelled (equation rewritten). R190 = CPA2 **and** CPA1. |
| Builder fix (code only) | (R765/R766 path) | The builder no longer removes YALI1E07744g, which is the transketolase. No model bytes change from this alone: the selection step re-adds the gene, and the control still reproduces a56f3649. |
| Transketolase fragment | R765, R766 | YALI1D02625g removed; GPR = YALI1E07744g |
| R153 | R153 | Peptidase YALI1D17462g removed; R153 is gene-less (gene kept, unlinked) |
| R2041 | R2041 | D-arabinokinase forward-only, [0, 1000] |
| R2190 | R2190 | Product written as D-ribulose (new species `m2061[C_cy]`) |
| R1964 | R1964 | Provisional GPR `YALI1F03452g or YALI1E06616g` (new gene YALI1E06616g); D-xylulose cross-references |

Each change's evidence, tier and limits are in its own report, linked from `model/STATE.md`
(entries dated 2026-10-19). The orphan inventory and the 33-reaction gene search are merged as
reports only.

## 2. How it was combined

- **Commits made for the plan:** the uncommitted hash-seed fix (great-wu) and the gene-search
  report were committed as their sessions left them. This session's branch was committed too.
- **Merges:** all eight branches were merged one by one into `claude/combined-candidate-20261019`,
  starting from `origin/main` `285eaee`.
- **Conflicts:** all were parallel additions.
  - `STATE.md`, both READMEs, `cli.py` and `main.py` kept both sides.
  - Every flag block stays at the stage its branch chose: R2041 after `apply_r1159_direction`;
    the others after the R132 step, in the order restores, fragment, R153, R2190, R1964.
  - For this session's branch, the code files were taken from the combined side, because
    `dc471d7` holds the same restore code, data and test plus later fixes.
  - STATE.md entries are newest first by date. Within 2026-10-19 they are in merge order.
- **Base has moved.** `origin/main` was `285eaee` when the merges were made. It has since moved to
  `1540a7d` (PR #7, gem_review literature layer), so this branch is 5 commits behind. Merging into
  main later will conflict at the top of `model/STATE.md`.
- **One test fix after the merge (§4):** `ORDER_PROBE` in
  `test_selection_gpr_restore_change_sets.py` now also applies the R765/R766 change set.

## 3. Verification

| Check | Result | Level |
| --- | --- | --- |
| Control build (a56f3649's flags, merged code) | `a56f3649…`, byte-identical | verified |
| Combined build, `PYTHONHASHSEED=0` and `=1` | Both `164f235f…`, byte-identical | verified |
| Build record | Input and sources unchanged, build complete. All 8 steps report `applied`. | verified |
| Union diff ([union_diff.json](union_diff.json), `scripts/union_diff.py`) | 50 field-level changes from a56f3649 = exactly the union of the seven part candidates' diffs. 0 extra, 0 missing, 0 value mismatches, 0 fields changed two ways. | verified |
| Counts | 2306 reactions (=), 1876 species (+1 D-ribulose), 1074 genes (+1 YALI1E06616g) | verified |
| Balance | All 13 touched reactions balance, except R671 (charge +1) and R765/R766 (charge +2). Those are GPR-only changes, and the same residuals are already in a56f3649. | verified |
| pytest | `origin/main`: 165 passed, 11 failed, 14 errors, 31 skipped. Combined: 207 passed, 11 failed, 14 errors, 31 skipped. Per test: identical plus 42 new passes from the merged branches ([logs](logs/)). | verified |

## 4. Deviation: one test failed after the merge (fixed with approval)

The first full test run on the merged tree failed `test_restored_species_order_does_not_depend_on_the_hash_seed`
with `KeyError: 'R765'`.
- That test (great-wu) checks every change set in `CHANGE_SETS`.
- After merging cool-nightingale, `CHANGE_SETS` also holds the R765/R766 fragment change set, but
  the probe applied only the first three change sets.
- Work stopped, and the user approved extending the probe to the R765/R766 change set.
- That change set adds no species, so it cannot make species order depend on the hash seed. It does
  change bytes (the R765/R766 GPRs and notes).
- `ORDER_PROBE` still lists the change sets by hand, so the next change set will need the same
  edit. Deriving the runs from `CHANGE_SETS` would prevent that; not done here.
- Rerun: identical to the baseline plus 42 new passes.

## 5. Solves (20 LPs, GLPK, SD-Leu + PO1f; observed in runs, [solves.json](solves.json))

| Test | a56f3649 | Combined |
| --- | --- | --- |
| WT growth | 1.430806 h⁻¹ | 1.430806 h⁻¹ |
| Loop test (all boundaries closed, forced fluxes relaxed) | R1964 −1000, R2041 −1000, R2042 +1000, R2190 −1000 ([loop_check.txt](../orphan_reactions_20261019/loop_check.txt); re-run by the audit) | **0** for all four |

Single-gene KO on the combined model (KO/WT):

| Gene | Reactions now | KO/WT | Consensus essential list |
| --- | --- | --- | --- |
| YALI1E07744g (transketolase) | R765, R766 | **≈0** (6.9e-10) | listed |
| YALI1C33005g (CPA2) | R190 | **0.1439** | listed |
| YALI1D09420g (CPA1) | R190 | **0.1439** | not listed |
| YALI1D02625g | none | 1.0 | not listed |
| YALI1E11768g (URA2) | R159 | 1.0 | not listed |
| YALI1E11370g (GatB) | none | 1.0 | listed (false negative; it has no reaction in the model) |
| YALI1D17462g | none | 1.0 | — |
| YALI1F03452g | R1964, R2190 | 1.0 | — |
| YALI1E06616g | R1964 | 1.0 | — |
| YALI1A21372g (LIP2) | R1895, R1896, R2274 | 1.0 | — |

"—" means not checked.

**What these mean:**
- **Arabinose loop gone.** The net-zero loop that ran an arabinose kinase backwards to remake ATP
  (R803 → R714 → R2042 → R2041 → R2190 → R1964) carries nothing in the combined model.
  - R2041 is now forward-only, and R2190 makes D-ribulose rather than D-arabinose.
  - R2190 carrying 0 holds by construction. Its new product D-ribulose (`m2061[C_cy]`) is used by no
    other reaction, so R2190 is blocked, as the R2190 report discloses.
  - The test does not say which change alone would close the loop; each part was built separately.
- **Transketolase.** The transketolase gene YALI1E07744g is now the only gene on R765/R766, and
  its KO is lethal.
  - In a56f3649 the OR with the fragment YALI1D02625g kept it non-essential.
  - The gene is on the consensus essential list (`1e887f5a…`, benchmark.md), so this call now
    agrees with the experiment.
  - The R765/R766 branch had not tested this (no LP). This run is its first observation.
- **CPA genes on R190.** KO/WT is 0.1439 for both, matching the selection-fix branch's own
  report. That is essential at the 15 % cutoff and non-essential at 1, 5 and 10 %.
  - YALI1C33005g is listed and YALI1D09420g is not, so at 15 % one call becomes a hit and one a
    false positive.
- **Not validation.** These changes were curated with the essentiality data in view (EG-GPR-001
  is an essentiality patch), so per benchmark.md they are not independent validation. No full
  benchmark was run.

**Solve provenance:**
- **Models:** baseline `a56f3649…`; combined `164f235f…` (the seed-1 file, byte-identical to the
  candidate).
- **Condition:** medium `model/reports/reference_pipeline_restore_20260909/research/state/media/sd_leu.csv`
  (`ed176d26…`) and profile `…/strain_profiles/po1f_sd_leu.json` (`35307853…`).
- **Software:** GLPK with presolve; cobra 0.30.0, optlang 1.8.3, Python 3.13.5.
- **Command:** TASK.md, "Exact commands".
- `solves.json` itself does not record these.

## 6. Independent audit

A fresh read-only reviewer got TASK.md, the merge commits, the deliverables and the logs. It
re-derived every claim:
- a libsbml element-level dump for the union;
- `--remerge-diff` of every merge;
- its own run of the loop test on both models;
- the KOs and the consensus list;
- a full pytest on copies of `285eaee` and of HEAD.

| Claim | Verdict |
| --- | --- |
| 1 Ancestry, merge resolutions, each flag block once and at its stage, nothing lost from b9eb2b6 | supported |
| 2 Control = a56f3649; seeds 0 and 1 byte-identical | supported (the seed values are not stored in the build records) |
| 3 Diff = exact union of the 7 parts | supported (110 attribute and note-level keys, 0 extra, 0 missing, 0 conflicts) |
| 4 Balance except inherited residuals | supported |
| 5 pytest = baseline + 42 new passes; test edit scoped | supported, reproduced |
| 6 20 LPs, WT, loop, KO ratios, consensus list | supported, with wording corrections |
| 7 Every guard `applied`; build record clean | supported |

**Corrections it asked for, all applied here or in TASK.md:**
- the moved `origin/main`;
- missing build commands and solve provenance;
- a broken table in `model/candidates/README.md`;
- the deviation wording ("cannot change … bytes");
- two parts built on `51c6fabd`, not directly on a56f3649;
- the code-only builder fix;
- the a56f3649 loop values (±1000, including R2190);
- R2190 blocked by construction;
- the YALI1E11370g consensus status;
- the raw YALI1E07744g value.

No scientific or model defect was found.

## 7. Limits

- **Benchmark:** no full single-KO screen or four-cutoff benchmark was run on the combined model.
  Only 10 genes were knocked out.
- **Loop test:** the test is a closed-boundary LP check. It shows no cycle through these four
  reactions; it is not a loopless FVA of the whole model.
- **Opt-in flags:** all eight changes stay behind their flags. The default build (`b4ce0974…`) and
  every earlier candidate are unchanged.
- **Gate:** the CoQ9 approval gate and the other orphan groups are untouched.
- **Part candidates:** the seven part candidates and their reports are kept as merged; they
  document each change on its own.

## 8. Corrections to earlier reports in this tree

- [selection_gpr_restore_20261019/REPORT.md](../selection_gpr_restore_20261019/REPORT.md) §6 and
  its STATE.md entry call YALI1E07744g a "known mis-annotation" that the builder removes. The
  later review (`model/reports/selection_gpr_review_20261019/`, merged here) found the 2026-07-06
  removal was a YALI1/YALI0 locus mix-up: YALI1E07744g is the transketolase.
  - The builder fix in `dc471d7` and the fragment change follow that finding.
  - The earlier text is left as written; this note supersedes it.
