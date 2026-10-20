# Combined candidate of the 2026-10-19 opt-in changes (2026-10-20)

## Authorization (user, this session)

> push the newest model

Asked which model, the user chose "Combine first, then push". The user then approved the plan:

> Approved: run the combine plan as written, including committing the great-wu sort fix and the
> gene-search report on their branches, and push claude/combined-candidate-20261019 (not main).

After the first post-merge test run failed one test (§ "Deviation"), the user approved the fix:

> Approved: extend ORDER_PROBE to apply the TKL change set with the tkl test fixture, rerun the
> tests, then continue the combine plan (seed-0 build, union diff, solves, audit, write-up, push).

Not authorized: pushing to `main`, opening a PR, changing any curation data or guard, or any
LP beyond the 20 below.

## Inputs

- **Baseline model:** `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`
  (`a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`). It is the newest
  candidate, not a designated reference.
- **Base commit:** `origin/main` `285eaee`.
- **Merged branches**, each user-approved on 2026-10-19 as an opt-in change built on a56f3649:

  | Branch (commit) | Change | Flag | Part candidate |
  | --- | --- | --- | --- |
  | `claude/great-wu-64e590` (`e6a86c4`, includes `dc471d7`) | 4-GPR restore (R612, R2274, R1302, R1303); R671 restore; EG-GPR-001 on R159/R190; builder fix for the R765/R766 removal; hash-seed fix | `--restore-selection-cleared-gprs`, `--restore-selection-overwritten-gprs`, `--restore-selection-cpa-ura2` | `f0f24b26` |
  | `claude/cool-nightingale-a2f471` (`7ae6cdd`) | YALI1D02625g fragment removed from R765/R766 | `--remove-selection-tkl-fragment-gpr` | `52d3106b` |
  | `claude/eloquent-borg-1bd1b0` (`87a43ab`) | R153 made gene-less | `--remove-r153-peptidase-gpr` | `d941d2cd` |
  | `claude/bold-chaplygin-1c16df` (`ae54dd8`) | R2041 forward-only | `--r2041-forward` | `da0c50ba` |
  | `claude/clever-austin-f6a04e` (`4e44571`) | R2190 product D-ribulose | `--r2190-ribulose-product` | `36f04193` |
  | `claude/reverent-wilbur-4dc5b6` (`56419bd`) | R1964 provisional GPR | `--r1964-gpr` | `db429743` |
  | `claude/orphan-reaction-newest-model-b8ddc8` (`b9eb2b6`) | Orphan inventory; 4-GPR restore (code also in `dc471d7`) | — | `c98961e9` |
  | `claude/unruffled-shannon-2baf44` (`9c7f691`) | Gene-search report only | — | — |

- **Commits made for the plan:**
  - `e6a86c4` on great-wu, holding that session's uncommitted hash-seed fix;
  - `9c7f691` on the gene-search branch;
  - `b9eb2b6` on this session's branch.

  Each holds that branch's work unchanged.
- **Condition:** SD-Leu + PO1f, the same files as the orphan report (`sd_leu.csv` `ed176d26…`,
  `po1f_sd_leu.json` `35307853…`). GLPK.
- **Interpreter:** the main checkout's `.venv` (Python 3.13.5, cobra 0.30.0, optlang 1.8.3).
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.

## Acceptance criteria

1. A control build with a56f3649's flags on the merged code reproduces `a56f3649…`.
2. The combined build with all 8 flags gives the same bytes under `PYTHONHASHSEED=0` and `1`.
3. Its whole-model diff from a56f3649 is exactly the union of the seven part candidates' diffs:
   nothing extra, nothing missing, and no field changed two ways.
4. Every touched reaction is balanced, or carries only a residual already in a56f3649.
5. pytest per test is the same as `origin/main` plus the merged branches' new tests passing.
6. At most 20 LPs, each recorded with its status.

**Stop:**
- a non-textual merge conflict;
- a control mismatch;
- seed-dependent bytes;
- a guard rejection;
- a non-union diff;
- a test regression.

## Deviation

The first post-merge test run failed
`test_selection_gpr_restore_change_sets.py::…::test_restored_species_order_does_not_depend_on_the_hash_seed`
with `KeyError: 'R765'`.
- That test (great-wu) checks every change set in `CHANGE_SETS`.
- After merging cool-nightingale, `CHANGE_SETS` also holds the R765/R766 fragment change set, but
  the test's probe applied only the first three change sets. Neither branch alone had both.
- Work stopped, and the user approved extending the probe to the TKL change set.
- The R765/R766 change is GPR-only (no species added or removed), so it cannot change species order
  or the model bytes.

## Corrections after the independent audit (2026-10-20)

The text above is kept as written. These corrections apply:

- **Part bases.** Not every part was built directly on a56f3649. `52d3106b` (R765/R766 fragment)
  and `d941d2cd` (R153) are built on `51c6fabd`, which is the selection fixes on a56f3649. The union
  check handles stacked parts.
- **Deviation wording.** "Cannot change species order or the model bytes" is wrong. The R765/R766
  change does change bytes (GPRs and notes). Correct statement: it adds no species, so it cannot
  make species order depend on the hash seed.
- **origin/main has moved.** It was `285eaee` when the merges were made. It is now `1540a7d`
  (PR #7, gem_review literature layer), which is not in this branch, so the branch is 5 commits
  behind. Both sides add entries at the top of `model/STATE.md`, so a later merge into main will
  conflict there.

## Exact commands (from `platform/`; `$V` = the main checkout's `.venv/bin/python`)

`$FLAGS` = `--offline --no-solve --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh`

```
PYTHONHASHSEED=0 $V -m scripts.gem_annotate --research-root "$RESEARCH" $FLAGS \
  --output-model ../.scratch/combined_20261020/control_no_flag.xml
PYTHONHASHSEED=1 $V -m scripts.gem_annotate --research-root "$RESEARCH" $FLAGS \
  --restore-selection-cleared-gprs --restore-selection-overwritten-gprs --restore-selection-cpa-ura2 \
  --remove-selection-tkl-fragment-gpr --remove-r153-peptidase-gpr --r2041-forward \
  --r2190-ribulose-product --r1964-gpr --output-model ../.scratch/combined_20261020/combined_seed1.xml
PYTHONHASHSEED=0 $V -m scripts.gem_annotate --research-root "$RESEARCH" $FLAGS <same 8 flags> \
  --output-model ../model/candidates/E5_coq9_resp_pkg_r132_combined_20261020.xml
IYALI26_RESEARCH_ROOT="$RESEARCH" $V -m pytest -q -rA -p no:cacheprovider
```

`$RESEARCH` = `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`. The build records
do not store `PYTHONHASHSEED`; the values above are the ones set on the command line.

From the worktree root, with `$D` = `model/reports/reference_pipeline_restore_20260909/research/state`:

```
$V model/reports/combined_candidate_20261020/scripts/union_diff.py \
  model/candidates/E5_coq9_resp_pkg_r132_20261012.xml .scratch/combined_20261020/combined_seed1.xml \
  union_diff.json <the 7 part candidates listed in logs/sha256.txt>
$V model/reports/combined_candidate_20261020/scripts/solves_combined.py \
  $D/media/sd_leu.csv $D/strain_profiles/po1f_sd_leu.json \
  model/candidates/E5_coq9_resp_pkg_r132_20261012.xml .scratch/combined_20261020/combined_seed1.xml solves.json
```

The union diff and the solves were run on the seed-1 file, which is byte-identical to the
candidate.
