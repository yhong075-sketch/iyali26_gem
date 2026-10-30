# R24 homocitrate dehydration to both aconitases (option C), 2026-10-30

**Status: done; see [REPORT.md](REPORT.md).** One `curation-change` round.

## 0. Authorization

- The user, 2026-10-30, in this session, after the review
  [r24_homocitrate_review_20261029](../r24_homocitrate_review_20261029/REPORT.md): "Option C", then "Approve plan".
  The approved plan:
  1. Commit the review on `claude/nice-bouman-b161f4` (done: `ac8921d`).
  2. Merge the local, unpushed branch `refs/heads/claude/aconitase-aco2-split-1c492f` (`a006df1`) into it, so that the
     builder carries `--aconitase-aco2-split` (done: merge `199dec2`).
  3. Implement option C through the pipeline: curated data, a guarded `apply_*` step, the flag
     `--r24-homocitrate-aconitase`, a focused test.
  4. Build a new candidate and a byte-identical control, diff the whole model, run the targeted solves, audit, report.
- **No push.** Option H (the 3 H+) and option N (R24 name and cross-references, m76 cross-references) are not in scope.

## 1. Starting state (verified 2026-10-30)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/nice-bouman-b161f4` @ `199dec2` (merge of `ac8921d` and `a006df1`; base `origin/main` `4794ec3`) |
| Baseline model | `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml`, `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897`. Not a designated reference. |
| Baseline build arguments | the `argv` of that candidate's `.build.json` |
| Medium / profile / list | Task A `po1f_csm_leu_w29uptake_20261012.csv`; PO1f `po1f_sd_leu.json`; `consensus_essential_genes.csv` (SHAs as in the review TASK.md) |

Targets in `e1e7c868` ([before.json](before.json), [logs/before_state.log](logs/before_state.log)):

| Item | Before | After (option C) |
| --- | --- | --- |
| R24 rule (C_mi, homocitrate ⇌ H2O + 3 H+ + m76, (−1000, 1000)) | `YALI1E03387g` | `YALI1D11984g or YALI1E17935g` |
| R24 notes | `PROTEIN_CLASS: 4.2.1.114` | the same plus this round's notes |
| `YALI1E03387g` (LYS4-like) reactions | R24, R395 | R395 |
| `YALI1D11984g` (ACO1-type) reactions | R1556, R216, R238, R239 | R1556, R216, R238, R239, R24 |
| `YALI1E17935g` (ACO2-type) reactions | none | R24 |

Nothing else may change: no name, bound, stoichiometry, formula, charge or annotation.

## 2. Expected effect and acceptance criteria

1. Control build (same arguments, no new flag) is byte-identical to `e1e7c868`.
2. Whole-model diff against `e1e7c868` shows only R24's rule and notes (gene-reaction links follow from the rule).
3. R24 stays mass- and charge-balanced (unchanged stoichiometry).
4. Targeted solves on Task A + PO1f, Gurobi 1e-9, reset per solve: WT and the knockouts of `YALI1E03387g`,
   `YALI1D11984g`, `YALI1E17935g` and `YALI1D11984g`+`YALI1E17935g` match the review's in-memory variant C
   ([variants.json](../r24_homocitrate_review_20261029/variants.json)): WT 0.3116608216; no class change at 1/5/10/15 %.
   No other gene's knockout can change, because only R24's rule changes and it names only these three genes.
5. Full pytest: only the new tests differ from the previous run.
6. Independent read-only audit: no unsupported claim left unaddressed.

## 3. Budget and stop condition

Two builds (candidate, control), plus a seed-1 candidate build; one whole-model diff; 5 LPs on the candidate (WT +
4 knockouts; the first draft said 6, corrected after the audit); two pytest runs; one audit agent.
Stop after REPORT.md and the STATE.md entry. The approved plan did not cover committing this round: commit only on
the user's instruction. Committed as `ff4dd32` on the user's instruction ("Commit on branch", 2026-10-30). No push.
