# Combined candidate c1f70ac4: R24 option C on the trehalose line

Scope, authorization and starting state are in [TASK.md](TASK.md). Verified in this run, 2026-10-31, unless marked.

## 1. Result

- **New candidate:**
  [E5_…_aco_split_treh_biomass_r24_aco_20261031.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_20261031.xml),
  `c1f70ac4198e443ce2bc089f534c021d4111d8c2e999888e5beae06cc7087540`.
  - It is `164542fb`'s arguments plus `--r24-homocitrate-aconitase`.
  - It is the newest model of the line: R24 option C on top of `164542fb`. That `164542fb` carries every model change from
    PR #34 to #40 is inferred from the candidates README chain (`164542fb` ← `e1e7c868` ← `059316ef` ← …), not re-traced.
  - It is not a designated reference.
- **Combination:** it equals `164542fb` with R24's change, and `f505571d` with the trehalose change, field for field.
- **Scores:** no knockout class and no knockout growth value changes against `164542fb`.
- **15 %, primary rule** (certified-infeasible as growth 0; mtDNA and untargeted genes out of scope):
  **136 / 190 / 43 / 695**, 0 unresolved.

## 2. What was merged

The branch `claude/combined-r24-trehalose-20261031` was cut from `origin/main` @ `74146c9` and received two merges:

| Merge | Branch | Conflicts and resolution |
| --- | --- | --- |
| `0b1a5f8` | `claude/nice-bouman-b161f4` @ `2ea1e31` (R24 option C, review and candidate `f505571d`) | Four files. In each, both sides added an item in the same place. `cli.py`: both flags kept. `main.py`: both parameters and both calls kept. Order: aconitase split, then R24, then trehalose. R24 must follow the aconitase split, because its guard expects YALI1E17935g to have no reactions. Its place before trehalose is a convention: the two guards cover disjoint fields (R24 and three genes; `biomass_C`, `m294` and the objective), so either order satisfies both (inferred from the code, not tested). The `0b1a5f8` commit message gives the incomplete reason ("the trehalose step must stay after every `biomass_C` edit"); this paragraph corrects it. `STATE.md`: main's entries kept, plus the branch's R24 entry and its R24/R395 review entry. `candidates/README.md`: main's rows kept, plus the `f505571d` row |
| `59e1547` | `claude/essential-gene-prediction-e76be3` @ `40212be` (V-ATPase report, no model change) | `STATE.md` only: the V-ATPase entry added above the rest |

Code and curated data that differ from `origin/main` (verified with `git diff --stat`):
- the R24 step `r24_homocitrate_aconitase.py`;
- its test;
- its curated file `model/curation/gpr/r24_homocitrate_aconitase.json`;
- 2 lines in `cli.py` and 8 in `main.py`.

## 3. Checks

| # | Check | Result |
| --- | --- | --- |
| 1 | Focused tests (R24, trehalose, aconitase) on the merged code | 17 passed, 77 subtests ([log](logs/pytest_focused.txt)) |
| 2 | Full pytest compared test by test with `origin/main` | see section 4 |
| 3 | Control build of `164542fb` from its recorded argv on the merged code | `164542fb…` byte for byte |
| 4 | Control build of `f505571d` from its recorded argv on the merged code | `f505571d…` byte for byte |
| 5 | Combined build under PYTHONHASHSEED 0 and 1 | both `c1f70ac4…` |
| 6 | Whole-model diff against `164542fb` | 5 entries: R24 GPR and notes; reaction links of YALI1D11984g, YALI1E03387g, YALI1E17935g ([log](logs/model_diff_vs_164542fb.log)) |
| 7 | Whole-model diff against `f505571d` | 2 entries: `biomass_C` stoichiometry and notes ([log](logs/model_diff_vs_f505571d.log)) |
| 8 | Changed values equal the parents' values | all 7 True: R24 rule `YALI1D11984g or YALI1E17935g` ([log](logs/field_identity.txt)) |
| 9 | Knockout screens of `164542fb` and `c1f70ac4` | WT 0.3122963651 on both; 0 class changes under rule A or B at any cutoff; largest knockout growth change 0.0; 0 points failing the check; 20 certified infeasible on each, none uncertified ([compare.json](compare.json)) |

The build commands are in [build_commands.txt](logs/build_commands.txt), the code state in
[build_head.txt](logs/build_head.txt) (HEAD `59e1547`, no tracked changes; the round folder and the new candidate files were untracked). The hash seed of each build is recorded only in `build_commands.txt`, and the input hashes in
[sweep_inputs_sha256.txt](logs/sweep_inputs_sha256.txt).

**Screen setup:**
- Task A medium + PO1f, Gurobi at 1e-9, model reset per knockout, targeted check.
- `ko_screen_certified.py`, copied unchanged from the infeasible-rule round (SHA-256 `313d5d01…`).
- 1,123 LPs per model, 2,246 in all.

## 4. Full pytest

The baseline ran in a detached worktree of `origin/main` @ `74146c9` at the short path `.claude/worktrees/m74`.
A first attempt from a `git archive` copy in the scratchpad hit cobra's 260-character path limit: relative to the `m74` baseline it has 58 newly failing tests (57 errors and 1 failure). Its other 14 errors also occur in the baseline. It is
kept as [the long-path attempt](logs/pytest_before_main_attempt1_longpath.txt) and was not used.

| Run | Failed | Passed | Errors | Subtests |
| --- | --- | --- | --- | --- |
| `origin/main` `74146c9` ([log](logs/pytest_before_main.txt)) | 19 | 420 | 14 | 1,031 |
| Combined branch ([log](logs/pytest_after_combined.txt)) | 19 | 420 | 20 | 1,008 |

- **The 19 failures are the same IDs on both.**
- **The 6 extra errors are path-length artifacts.** They are all in `test_dolichol_route.py`. This worktree's longer
  path makes the dolichol baseline path 268 characters, over cobra's limit, against 239 in `m74`. Re-run on the
  merged commit `59e1547` from the short-path worktree `.claude/worktrees/c59`, all 6 pass, with 52 subtests
  ([log](logs/pytest_dolichol_shortpath_59e1547.txt)).
- **The counts reconcile exactly:**
  - passed: 420 − 6 dolichol + 6 new R24 tests = 420;
  - subtests: 1,031 − 52 + 29 = 1,008. The 29 R24 subtests come from a separate run of the R24 test file on the merged code: 6 passed, 29 subtests ([log](logs/pytest_r24_only.txt));
  - errors: 14 + 6 = 20.

  See [pytest_compare.txt](logs/pytest_compare.txt). So no test regresses on the merged code.

## 5. Scores (c1f70ac4, identical to 164542fb)

Negatives are the user-defined class (model genes not on the consensus list). mtDNA genes are out of scope. This
is a development reference.

| Cutoff | Primary: rule B, mtDNA + untargeted out of scope | Earlier rule A, mtDNA out of scope |
| --- | --- | --- |
| 1 % | 108 / 218 / 36 / 702 | 97 / 218 / 32 / 702, 16 unresolved |
| 5 % | 115 / 211 / 38 / 700 | 104 / 211 / 34 / 700, 16 unresolved |
| 10 % | 119 / 207 / 43 / 695 | 108 / 207 / 39 / 695, 16 unresolved |
| 15 % | **136 / 190 / 43 / 695** | 125 / 190 / 39 / 695, 16 unresolved |

At 15 % under the primary rule: recall 0.4172, precision 0.7598, MCC 0.4423.

**Coverage:**
- 326 of the 1,612 listed positives are model genes; the other 1,286 are out of scope as listed genes with no model
  counterpart, so they are never counted as FN.
- Out-of-scope model genes: the 16 mtDNA genes and 1 untargeted gene (YALI0F04114g, COX9). All 17 are non-positives,
  so the scope rules remove negatives only.

The infeasible-rule round reported 136 / 190 / 46 / 692 on `e1e7c868`. The 3 fewer FP here are the trehalose genes,
which moved FP → TN in PR #40. This round's screen of `164542fb` uses the same script and inputs as that round (verified). No other summary under
`model/reports` carries the `164542fb` SHA, so it is probably the first under the new rule (the research root was not
scanned).

The two TSVs are identical in every class, status and growth column. One cell differs: the numerical `max_residual`
of YALI1E03387g, 4.7e-14 vs 8.0e-15.

## 6. Meaning and limits

- **Why nothing moves:** R24 option C changed no class on `e1e7c868`, and the trehalose term does not involve R24.
  So the combination was expected to score like `164542fb`, and it does, down to every knockout growth value.
  That the two edits are independent was inferred beforehand; checks 6 to 9 now show it.
- **Proved vs observed:**
  - Byte identity, field identity and the diffs are proved for these files.
  - The scores are observed in this run on one medium and one strain profile.
- **What still holds from the parent rounds:** the R24 option C evidence, and its open questions on the LYS4-like
  YALI1E03387g, belong to the R24 report. The trehalose evidence belongs to the trehalose report.
- **Not done:** no new biology, no FVA, no other medium.
- **Bookkeeping:** two detached worktrees were created for the short-path tests: `.claude/worktrees/m74` (`74146c9`)
  and `.claude/worktrees/c59` (`59e1547`). Both are left in place for the user to remove.

## 7. Independent audit

A fresh read-only reviewer, given only this folder, the logs, the deliverables and the branch history, checked every
claim on 2026-10-31. It found the merges, the byte identity, the diffs, the screens and the pytest reasoning supported.
It re-ran the comparisons with its own code. Its corrections are applied above:
- the reason for the step order;
- the long-path attempt count;
- coverage and out-of-scope reporting;
- the inferred lineage;
- "no tracked changes".

Its scripts are in this session's scratchpad (`audit/`), not in the repository.
