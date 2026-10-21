# Task: deterministic species order in the selection restore (2026-10-19)

## Authorization

The user approved this in this session on 2026-10-19:

> Approved: fast-forward claude/great-wu-64e590 to dc471d7, add the cross-seed test, apply the
> sorted() fix at line 154, run pytest before/after, and do the two seeded offline builds (new
> candidate + .scratch). No commit or push.

## Scope

- Code: `platform/scripts/gem_annotate/selection_gpr_restore.py`, line 154 only. The restore must
  iterate in sorted ID order, so that species it adds go into the reaction in an order that does not
  depend on `PYTHONHASHSEED`.
- Test: one new test in `platform/tests/test_selection_gpr_restore_change_sets.py`.
- Computation:
  - full pytest before and after the fix;
  - two offline, no-solve builds with the dc471d7 argv, one with `PYTHONHASHSEED=0` (new candidate
    file) and one with `PYTHONHASHSEED=1` (`.scratch/`).
- Out of scope: model science, curation data, other builder steps, commit, push.

## Starting state

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/great-wu-64e590`, fast-forwarded from `285eaee` to `dc471d75480bd688280df0bed5af2db8bbd7cebb` |
| Code under change | `selection_gpr_restore.py` SHA-256 `9bf0f5a63f246f7c122e0c3a9110edf335dccb6c6464809c9b1e96adbc276be7` |
| Baseline model (byte reference) | `model/candidates/E5_coq9_resp_pkg_r132_selection_gpr_fix_v2_20261019.xml`, SHA-256 `51c6fabd486f4e1a9d6365bf7740bde3eca20fa05f1d5454a73586d3276db3fa` |
| Start model | `model/source/iyali26.xml` (as recorded by the build) |
| Research root | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` |
| Interpreter | main checkout `.venv`: Python 3.13.5, cobra 0.30.0 |

## Stop conditions

Stop and report if any of these happens:

- pytest changes for any existing test;
- the two seeded builds differ;
- `model_diff.py` against 51c6fabd reports any difference;
- the bytes differ anywhere outside the `R_R159` speciesReference order.
