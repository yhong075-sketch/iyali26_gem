# TASK: make R153 gene-less (option A of the R153 ASS review), 2026-10-19

## Authorization (quoted)

- 2026-10-19: "Approve option A for R153: a gene-less GPR through a new guarded post-merge step and opt-in
  flag, with YALI1D17462g kept unlinked. Show me the implementation plan before editing code."
- 2026-10-19: "Approve the R153 implementation plan as written, including the fast-forward of
  claude/eloquent-borg-1bd1b0 to dc471d7. Builds and pytest as budgeted; no LPs; no commit."
- Note: the detailed plan text was not displayed when this approval was given. It was shown afterwards, and
  code edits wait for the user's OK on it.

Evidence for the decision: [../r153_ass_gpr_20261019/REPORT.md](../r153_ass_gpr_20261019/REPORT.md).

## Baseline

| Item | Value |
| --- | --- |
| Branch | `claude/eloquent-borg-1bd1b0`, fast-forwarded 285eaee → `dc471d75480bd688280df0bed5af2db8bbd7cebb` (approved) |
| Baseline model | `model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml`, SHA-256 `51c6fabd486f4e1a9d6365bf7740bde3eca20fa05f1d5454a73586d3276db3fa` |
| Baseline build argv | `--offline --no-solve --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh --restore-selection-cleared-gprs --restore-selection-overwritten-gprs --restore-selection-cpa-ura2` |
| Baseline target state | [before.json](before.json) |
| Research root | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` (exported as `IYALI26_RESEARCH_ROOT` for builds and tests) |
| Interpreter | main checkout `.venv` (as in the selection-fix run) |

## Targets

- **Reaction R153** (argininosuccinate synthase, C_cy, [-1000, 1000]): GPR `YALI1D17462g` → empty.
  - Bounds, stoichiometry, species, annotation and groups are unchanged.
  - Notes: 6 keys are overwritten (old values kept in the curation file) and 8 keys are added.
- **Gene `YALI1D17462g`:** stays in the model with no reactions. Its annotation is unchanged.
- **Nothing else.** No metabolite, other reaction, other gene, bound or objective changes.

## Expected effect

- R153 can no longer be disabled by any gene KO.
- The KO of `YALI1D17462g` disables no reaction.
- The LP is unchanged: a GPR is not part of the WT problem, and bounds are untouched. This is provable;
  no LP is run (not authorized).

## Acceptance criteria

1. **Control.** The flag-off build at `dc471d7` reproduces `51c6fabd…` byte for byte.
2. **Diff.** The flag-on build writes a new file. The whole-model diff against `51c6fabd…` is confined to
   R153's `gpr` and `notes`, and to `YALI1D17462g.reactions` going from `[R153]` to `[]`.
   - Gene count stays 1,073.
   - R153 stays element- and charge-balanced (it is balanced in the baseline).
3. **Focused test.** The new in-memory test passes:
   - each guard fails closed with the model unchanged;
   - the step applies, then a rerun returns `already_correct`;
   - the gene is kept;
   - the flag defaults off.
4. **pytest.** Per test, the full run is identical to the baseline run, except the new tests, which pass.
   The baseline already fails some tests: 8 failures and 14 errors in the selection-fix environment; 11
   failures and 14 errors as measured here, with `IYALI26_RESEARCH_ROOT` set.
5. **Audit.** An independent read-only audit returns no unresolved "unsupported" verdict.

## Baseline runs (2026-10-19, before any code edit)

- **Full pytest** with `IYALI26_RESEARCH_ROOT` set: 181 passed, 11 failed, 14 errors, 31 skipped
  ([logs/pytest_before.txt](logs/pytest_before.txt)).
  - The selection-fix run recorded 177 passed, 8 failed, 14 errors, 38 skipped without that variable.
  - 7 vacuole tests un-skip when the variable is set: 4 pass, 3 fail. These failures predate this work; the
    before/after comparison here uses the same environment.
- **Control (no new flag): stop condition hit.**
  - SHA-256 `8f083627…` ≠ `51c6fabd…`.
  - `model_diff.py` finds **0 differences** in every field.
  - At byte level, only R159's reactant `speciesReference` order differs.
  - Cause, from reading the code: `selection_gpr_restore.py:154` iterates a hash-randomized `set` of species
    IDs. See [logs/control_check.txt](logs/control_check.txt).
  - So criterion 1 as written (byte identity) cannot be met reliably. Work stopped for the user's decision.
- **User decision (2026-10-19, quoted; the message was cut at 300 characters):** "Proceed with the R153 code
  plan from your previous message. Amend criterion 1: the control passes if model_diff.py shows 0 differences
  vs 51c6fabd and bytes differ only in R159 speciesReference order. Rebuild the control and build the
  candidate both with PYTHONHASHSEED=0 (one extra build), so their"
  - Read with the offered button text: "…so their byte diff isolates R153. No LPs, no commit."

**Amended criterion 1.**
- The control passes when `model_diff.py` shows 0 differences against `51c6fabd…`, and its bytes differ only
  in R159's `speciesReference` order.
- The control and the candidate are both built with `PYTHONHASHSEED=0`. The byte diff between them must be
  confined to R153 (its GPR association and notes).
- Budget: one extra build.
- The fixed-seed nondeterminism fix itself is in a separate session
  (`local_6e538fd6…`, "Make selection-restore R159 species order deterministic").

## Compute budget and stop condition

- Budget:
  - 2 offline/no-solve builds (control no-flag; candidate), about 70 s each, plus at most one rebuild of
    each after audit fixes;
  - 2–3 full pytest runs, about 4 min each;
  - 0 LPs; no cluster jobs.
- Stop and report on any of:
  - a guard failure on the real build;
  - any diff outside the targets;
  - a control that does not reproduce `51c6fabd…`;
  - any per-test pytest change other than the new tests.
- No commit or push (not authorized).
