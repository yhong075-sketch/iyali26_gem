# Independent audit of the R24 option C round, 2026-10-30

**Setup.**
- Reviewer: a fresh general-purpose agent with no conversation history, read-only (no edits, builds, solves or state
  changes).
- Given: TASK.md, REPORT.md (first draft), the diffs, logs and deliverables, the review folder and AGENTS.md.
- Below is a summary of its verdicts. The last column of the issues table says what was done about each one.

## Verdicts per claim

| # | Claim | Verdict | Deciding evidence (as reported by the auditor) |
| --- | --- | --- | --- |
| 1 | Baseline `e1e7c868…`, candidate `f505571d…` | supported | `shasum` on both files |
| 2 | Control = baseline; seed 1 = seed 0 | supported | `logs/builds_sha256.txt` (lines unlabelled; the output paths name them) |
| 3 | Build argv = `e1e7c868` argv | supported | `args_e1e7c868.txt` vs the baseline `.build.json` |
| 4 | Exactly 5 diffs; R24 notes only gained keys | supported | `model_diff_vs_e1e7c868.json`; the auditor's raw XML diff: 6 added `<p>` note lines and the R24 GPR, nothing else |
| 5 | Build record: HEAD `199dec2`, step `applied` | supported | `.build.json` `reaction_selection`; `source_sha256` of the step, data, cli, main |
| 6 | R304 libSBML error pre-existing | supported | all build logs here and in the aconitase round |
| 7 | Every check runs before any edit; only R24 edited | supported, with a caveat | the post-edit `RuntimeError` fires after editing but cannot be reached in practice |
| 8 | Each guard is tested | partially supported | no missing-gene case; no gene guard on the applied path; one misleading comment |
| 9 | Sorting of "R239" vs "R24" | supported | both sides use `sorted()` |
| 10 | Context exit restores rule and notes | supported | code and `test_context_exit_restores_the_baseline` |
| 11 | Must run after the aconitase split, which fails closed on `f505571d` | supported (code read) | both gene guards; `aco2_still_linked` case |
| 12 | Knockout table, 5 LPs, all optimal and passed | supported | `candidate_kos.json` |
| 13 | Equals the review's variants C and A | supported | `kos_vs_review_variant_c.txt`; only difference 8.67e-19 |
| 14 | `YALI1E03387g` KO is set by the lysine bound 0.00595 | supported | medium CSV and log |
| 15 | R24 carries no flux in the ACO1-type KO | partially supported | one optimum, no FVA |
| 16 | Erratum: cause of the stale knockout | supported | `repro_functional.log` |
| 17 | Erratum: only variant C's ACO1-type row affected | supported (auditor's reasoning over `variants.json`) | — |
| 18 | No other gene can change class | supported (reasoning plus diff) | — |
| 19 | 15 % score 125 / 190 / 42 / 692 | unverified here | labelled "not recomputed"; matches the README row |
| 20 | pytest counts | supported | log summary lines |
| 21 | Only the 6 new tests plus 29 path errors changed | supported | identical FAILED sets; after-errors ⊂ before-errors; all 43 before-errors are cobra `OSError` |
| 22 | 23 pass in the aconitase run; 6 dolichol errors; 263 vs 254 characters | supported | comparison log; the auditor's own length check |
| 23 | Python / cobra / pytest versions | partially supported | no version header in the pytest logs |
| 24 | Evidence statements match review rows E1–E11 | supported | evidence files |
| 25 | "Converted **only** by aconitase plus homoaconitase … verified" | partially supported | E2: the yeast part is "not shown" |
| 26 | The ACO2-type member is not overstated | supported | tier says transfer from *S. cerevisiae*; if anything understated (E4 genetics) |
| 27 | TASK budget 6 LPs | partially supported | the report says 5 |
| 28 | TASK criterion 5 (only new tests differ) | not met as written | 29 others differ, cause explained |
| 29 | AUDIT.md and STATE.md entry exist | unsupported (at audit time) | — |

## Issues and resolution

| # | Severity | Issue | Resolution |
| --- | --- | --- | --- |
| 1 | medium | Criterion 5 reported as passing although it does not hold literally | REPORT §5 now says "not met as written" and gives the cause of each of the 29 differences. A clean same-path before run was not made. |
| 2 | medium | AUDIT.md and the STATE.md entry are missing | this file; STATE.md entry added |
| 3 | low | REPORT §1 overstates "only … verified" against E2 | REPORT §1 reworded. The curated note keeps its wording, because changing it needs a rebuild; this is recorded in REPORT §7. |
| 4 | low | Test gaps (missing gene; gene guard on the applied path; comment) | 3 subcases added and the comment corrected; focused run 12 passed, 58 subtests (`logs/pytest_focused2.txt`). The full runs used the earlier test file. |
| 5 | low | TASK says 6 LPs, REPORT says 5 | TASK §3 corrected |
| 6 | low | "carries no flux" rests on one optimum | reworded "at this optimum (no FVA was run)" |
| 7 | low | Provenance: unlabelled hash lines, empty log, no version header, compare counts | build lines labelled in REPORT §3; the empty `pytest_before_rerun_copied.txt` (this round's own) removed; `logs/pytest_env.txt` added; compare counts explained in REPORT §5 |
| 8 | info | Scope notes: stale-knockout reach, deleted UniProt IDs, order-sensitive rule guard | added to REPORT §7 |
