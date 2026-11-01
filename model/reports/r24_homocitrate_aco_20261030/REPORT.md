# R24 homocitrate dehydration to both aconitases (option C): report, 2026-10-30

**Outcome.** New candidate
[E5_…_aco_split_r24_aco_20261030.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_r24_aco_20261030.xml),
SHA-256 `f505571d3da6d76920f242e334aa74521d0f754c452b65eb12791e89c079e3c6`. It is `e1e7c868` plus
`--r24-homocitrate-aconitase`. R24's rule changes from `YALI1E03387g` to `YALI1D11984g or YALI1E17935g`, and R24 gets
six notes. Nothing else in the model changes. No knockout class changes. Not a designated reference. Scope and
authorization: [TASK.md](TASK.md).

## 1. What changed and why

| Item | Before (`e1e7c868`) | After (`f505571d`) |
| --- | --- | --- |
| R24 rule (C_mi, homocitrate ⇌ H2O + 3 H+ + m76) | `YALI1E03387g` | `YALI1D11984g or YALI1E17935g` |
| R24 notes | `PROTEIN_CLASS: 4.2.1.114` | unchanged, plus 6 `r24_homocitrate_aco_*` notes |
| `YALI1E03387g` (LYS4-like) | R24, R395 | R395 |
| `YALI1D11984g` (ACO1-type) | R1556, R216, R238, R239 | the same plus R24 |
| `YALI1E17935g` (ACO2-type) | no reaction | R24 |

The evidence comes from the review ([r24_homocitrate_review_20261029](../r24_homocitrate_review_20261029/REPORT.md),
rows E1–E11). Evidence was not re-fetched this round.
- **`YALI1E03387g` does not catalyse R24.** Yeast and *A. fumigatus* homoaconitases show no homocitrate activity
  (Fazius 2012, [PMC3556520](https://pmc.ncbi.nlm.nih.gov/articles/PMC3556520/), Table 1; verified from the full
  text in the review). In the coupled assay, homocitrate was converted by aconitase plus homoaconitase. That was
  shown with *A. fumigatus* AcoA; with yeast Aco1p or Aco2p the authors state it but the data are "not shown" (E2).
- **`YALI1E03387g` stays on R395.** Its Swiss-Prot entry [Q6C791](https://rest.uniprot.org/uniprotkb/Q6C791.json)
  lists only RHEA:15485, the R395 step, by similarity. Verified in the review.
- **Either aconitase is enough in yeast.** In *S. cerevisiae*, *aco1* and *aco2* single mutants grow without lysine.
  The double mutant is completely lysine auxotrophic (same paper). Verified in the review. The OR rule transfers
  that pattern to *Y. lipolytica*: inferred.
- **The ACO2-type member is the weaker half of the rule.** No *Y. lipolytica* data exist for it (E10). Its
  Aco2-clade lysine is shared by the inactive *Aspergillus* AcoB proteins (E9). The curated notes say so
  (`r24_homocitrate_aco_evidence_tier`).

## 2. How it is implemented

| Piece | File | SHA-256 (first 16) |
| --- | --- | --- |
| Curated data | [model/curation/gpr/r24_homocitrate_aconitase.json](../../curation/gpr/r24_homocitrate_aconitase.json), written by [scripts/write_curation.py](scripts/write_curation.py) from [before.json](before.json) | `f2f1b9d03372c6ac` |
| Guarded step | [platform/scripts/gem_annotate/r24_homocitrate_aconitase.py](../../../platform/scripts/gem_annotate/r24_homocitrate_aconitase.py) | `15e4625c19d4b14a` |
| Test | [platform/tests/test_r24_homocitrate_aconitase.py](../../../platform/tests/test_r24_homocitrate_aconitase.py) | `7765013775b5ccfd` (after the audit; the full runs used `60aaedc6f5defca3`, §5) |
| Wiring | `--r24-homocitrate-aconitase` in `cli.py`; `main.py` runs the step right after the aconitase split | — |

**What the step checks before editing** (all checks run first; any mismatch raises and nothing is edited):
- R24's name, bounds, stoichiometry, species (formula, charge, compartment) and annotation.
- R24's rule is either the before rule or the after rule.
- No conflicting notes, and no notes outside this curation.
- Notes and rule match, so a half-applied state is rejected.
- For each of the three genes: its annotation, and its reaction list before or after.

**Why the step must run after `--aconitase-aco2-split`.** The gene guard pins `YALI1E17935g` with no reactions and
`YALI1D11984g` on R1556, R216, R238, R239. That is the state after the aconitase split, so without that flag the
step fails closed (test case `aco2_still_linked`). It also works the other way: the aconitase split's own guard
expects `YALI1E17935g` with no reactions. So `apply_aconitase_aco2_split` fails closed on a model where this step
has already run. That is correct for re-applying it to `f505571d`, and it means the order in `main.py` must stay as
it is.

## 3. Build and whole-model verification

**Builds** ([scripts/run_build.zsh](scripts/run_build.zsh) with the `e1e7c868` argv in
[scripts/args_e1e7c868.txt](scripts/args_e1e7c868.txt); logs in [logs/](logs/); hashes in
[logs/builds_sha256.txt](logs/builds_sha256.txt)):

| Build | SHA-256 |
| --- | --- |
| Control: `e1e7c868` argv, seed 0, after the `cli.py` / `main.py` edits, flag off | `e1e7c868…` (byte-identical to the baseline) |
| Candidate: + `--r24-homocitrate-aconitase`, seed 0 | `f505571d3da6d76920f242e334aa74521d0f754c452b65eb12791e89c079e3c6` |
| Candidate, seed 1 | identical |

- Command, from `platform/`: `PYTHONHASHSEED=0 /opt/miniconda3/bin/python -m scripts.gem_annotate <args_e1e7c868.txt>
  --r24-homocitrate-aconitase --output-model ../model/candidates/<name>.xml`.
- The build record ([.build.json](../../candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_r24_aco_20261030.build.json))
  shows HEAD `199dec2` with this round's files uncommitted.
- The lines of `builds_sha256.txt` are, in order: control seed 0, candidate seed 0, candidate seed 1 (the output
  paths name them). Its `reaction_selection.r24_homocitrate_aconitase` reads
  `applied`.
- Every build log, the control included, carries the libSBML error on R304's notes. It is pre-existing: the
  aconitase round's logs have it too.

**Whole-model diff against `e1e7c868`** ([model_diff_vs_e1e7c868.json](model_diff_vs_e1e7c868.json),
[logs/model_diff.log](logs/model_diff.log)):
- Run with the aconitase round's `model_diff.py` (`8a464ff8…`), on SHA-checked short copies.
- **5 differences, all intended:**
  - R24 `gpr`;
  - R24 `notes`: 6 keys added, none removed or changed;
  - the reaction lists of `YALI1D11984g`, `YALI1E03387g` and `YALI1E17935g`.
- No stoichiometry, bound, formula, charge, annotation or name changed.
- So R24's mass and charge balance is the same as in the baseline. Its 3 H+ artefact is untouched (review option H,
  not in scope).

## 4. Targeted solves

**Setup.** [scripts/r24_candidate_kos.py](scripts/r24_candidate_kos.py), run by
[scripts/run_kos.zsh](scripts/run_kos.zsh) on a SHA-checked short copy of `f505571d`:
- Task A medium + PO1f profile; Gurobi at 1e-9, reset before each solve, targeted check (`checked_solve.py`
  unchanged).
- Input hashes in [logs/ko_inputs_sha256.txt](logs/ko_inputs_sha256.txt). Output in
  [candidate_kos.json](candidate_kos.json) and [logs/candidate_kos.log](logs/candidate_kos.log).
- 5 LPs, all optimal, all passed the check.

| Knockout | Reactions closed | KO/WT | Class 1 / 5 / 10 / 15 % |
| --- | --- | --- | --- |
| none (WT) | — | growth 0.3116608216 | — |
| `YALI1E03387g` (listed) | R395 | 0.0694 | FN / FN / TP / TP |
| `YALI1D11984g` (listed) | R1556, R216, R238, R239 | 0.0244 | FN / TP / TP / TP |
| `YALI1E17935g` (not listed) | none | 1.0 | TN / TN / TN / TN |
| `YALI1D11984g` + `YALI1E17935g` | R1556, R216, R238, R239, R24 | 0.0244 | — |

- Classes: essential when KO/WT < cutoff. Negatives are the user-defined class: model genes not on the
  experimental list. The list was used in development, so these are development scores, not validation.
- WT and all four classes equal the `e1e7c868` values and the review's variant C
  ([logs/kos_vs_review_variant_c.txt](logs/kos_vs_review_variant_c.txt)).
- **No other gene can change class.** Only R24's rule changed, and it names only these three genes. Verified by the
  whole-model diff, so no full sweep was run.
- So the score at every cutoff is that of `e1e7c868`. For example, at 15 %: TP / FN / FP / TN = 125 / 190 / 42 / 692
  (from the aconitase round; not recomputed).

**Mechanism.**
- The `YALI1E03387g` knockout still closes R395. Growth is set by the medium's lysine bound: uptake R1220 at
  0.00595.
- In the `YALI1D11984g` knockout, R24 now stays open through `YALI1E17935g` but carries no flux at this optimum
  (no FVA was run), and growth does not change.
- R24 closes only when both aconitases are out. That double knockout grows exactly like the ACO1-type single knockout.

**Erratum to the review, found here.**
- In the review's in-memory variant C, the single `YALI1D11984g` knockout also closed R24. That is wrong for an OR
  rule. Here, on the built candidate, R24 stays open.
- **Cause (verified, [scripts/repro_functional.py](scripts/repro_functional.py),
  [logs/repro_functional.log](logs/repro_functional.log), cobra 0.30.0):**
  - After an SBML read, a gene with no reactions (`YALI1E17935g` in `e1e7c868`) is not linked to the model.
  - So `knock_out()` inside `with model:` is not undone on exit.
  - The review's variant A knocked it out, and it then stayed non-functional through variants B–E.
- **Effect:** only that one row was affected, and its growth equals the correct value (difference 9e-19). No number,
  class or score in the review changes.
- The review REPORT §2 now carries the erratum; its `variants.json` is kept as run.

## 5. pytest

Python 3.13.5, cobra 0.30.0, pytest 9.0.2 (the venv of the workflow-automation-review worktree;
[logs/pytest_env.txt](logs/pytest_env.txt), recorded after the runs from the same interpreter).

| Run | Where | Result |
| --- | --- | --- |
| Before ([logs/pytest_before.txt](logs/pytest_before.txt)) | `git archive 199dec2` under the scratchpad | 386 passed, 19 failed, 43 errors, 31 skipped |
| After ([logs/pytest_after.txt](logs/pytest_after.txt)) | this worktree | 421 passed, 19 failed, 14 errors, 31 skipped |
| Focused ([logs/pytest_focused.txt](logs/pytest_focused.txt)) | this test plus the aconitase test | 12 passed, 55 subtests |
| Focused after the audit ([logs/pytest_focused2.txt](logs/pytest_focused2.txt)) | the same, with 3 added subcases | 12 passed, 58 subtests |

**Test by test** ([logs/pytest_compare.txt](logs/pytest_compare.txt),
[logs/pytest_compare_vs_aconitase_after2.txt](logs/pytest_compare_vs_aconitase_after2.txt)):
- **6 new tests, all passing.**
- **29 tests error in the before run and pass after.** All 29 fail on cobra's 260-character path limit, or on a
  sibling `iyali26_gem/model.xml` that does not exist under the scratchpad
  ([logs/pytest_path_length_note.txt](logs/pytest_path_length_note.txt)).
  - 23 of them also pass in the aconitase round's full run, from a worktree path.
  - The other 6 (`test_dolichol_route`) error there too, for the same reason: the path is 263 characters in that
    worktree and 254 in this one.
- **No other test changed outcome.** The 19 failures and 14 errors are the same tests in both runs.

**TASK criterion 5 is not met as written.** Besides the 6 new tests, 29 tests differ, each traced to the OSError path
cause above. This is not a clean before run at the same path; the no-regression claim rests on two runs made under
different path lengths, plus the cause of each difference.

- The test file gained 3 subcases after the full runs (audit issue 4: a missing gene, and the gene guard on an
  applied model). Only the focused run covers that version.
- `pytest_compare.txt` counts only PASSED / FAILED / ERROR lines (449 before, 455 after). The 31 SKIPPED lines have
  another format (`SKIPPED [1] file:line`) and are not compared; the skip count is 31 in both runs.

## 6. Independent audit

A fresh read-only reviewer audited TASK.md, the diffs, logs and deliverables: [AUDIT.md](AUDIT.md). No bug
changes the model, the knockout table or the erratum. The fixes it asked for are applied (AUDIT.md, last column).

## 7. Limits and what was not tested

- **One medium and one profile** (Task A + PO1f). No condition where the ACO2-type gene alone would matter was
  tested. On this medium, option C and option B give the same knockouts, so this round adds no testable difference.
- **`YALI1E17935g` activity on homocitrate is unproven** in *Y. lipolytica*. Its place in the OR rule rests on
  *S. cerevisiae* genetics and homology. If later evidence shows it is inactive, the fix is option B, a one-member
  rule.
- **Not changed, by scope:**
  - R24's name "2-methylcitrate dehydratase", its wrong cross-references, `ec-code` 4.2.1.36 and `PROTEIN_CLASS`
    4.2.1.114 (review option N);
  - m76's cross-references (option N);
  - the 3 H+ on R24 and R395 (option H).
- **Wording conflict.** The aconitase round's curation file and step docstring say `YALI1E17935g` is "kept in the
  model with no reaction". That is true of `e1e7c868`, but in `f505571d` the gene is on R24. Those files were left
  as they are because they describe their own round; this report records the change.
- **The curated note** `r24_homocitrate_aco_evidence` says homocitrate "is converted by aconitase plus homoaconitase"
  without the "yeast data not shown" caveat. The JSON evidence row carries it. Changing the note would need a
  rebuild, so it was left as built.
- **The stale-knockout behaviour of cobra** (§4) affects any script that changes rules between knockouts on a model
  that has reactionless genes. Fixed-model sweeps are unaffected, because no rule names such a gene.
- **The model's own UniProt IDs** for the two aconitases (A0A1D8NDY4, A0A1H6Q7Y5) are deleted in UniProt (aconitase
  REPORT). The notes cite Q6C9P6 and Q6C5V2 instead.
- **The rule guard is order-sensitive.** A rule written `YALI1E17935g or YALI1D11984g` fails closed.
- **Not pushed.** The aconitase round reached this branch only by a local merge; it is not on `origin/main`.
