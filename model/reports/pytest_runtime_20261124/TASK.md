**Status: changes 1 and 2 done (section 2). The three review follow-ups are done (section 4): five of their six
criteria are met, and B2 was not met as written. The user settled that by choosing the variable's type in each
term (section 5), which is done with all six of its criteria met. No decision is open. Change 3 is not authorized
and not started. Committed as `655b9b68` on branch `claude/pytest-performance-analysis-3e0794`; pull request
to `main`. Results: [REPORT.md](REPORT.md).**

# Shorter pytest runs: parallel test files and a faster state snapshot, 2026-11-24

A tooling round on branch `claude/pytest-performance-analysis-3e0794`. No model, curation file, builder step or test
assertion changes. The only LPs solved are the ones the test suite already solves.

## 0. Authorization

The user, 2026-11-24, in this session:

1. **The question.** "The pytest usually takes long time, can you explain by code, and how can we make it shorter?"
   That covered reading, measuring and full test runs on copies outside the repository.
2. **The choice.** After the findings and three proposed changes, asked "Which of the three proposed changes should I
   implement?", the user answered: "1 and 2".
3. **The follow-ups.** After the independent review, asked "Which of the review's three follow-ups should I apply?",
   the user answered: "All three". Their scope and criteria are in section 4.
4. **The type.** Asked "Keep the bounds in the snapshot terms as they are now?", the user answered: "Add the type
   too". Scope and criteria are in section 5.
5. **The commit.** "Commit the working tree changes with a sensible message." Commit `655b9b68`.
6. **The pull request.** The user chose "Create PR" in the app: push the branch and open a pull request to `main`
   of `yhong075-sketch/iyali26_gem`.

In scope:

1. **Run test files in parallel.**
   - A hook in `platform/tests/conftest.py` that lets pytest-xdist send sub-test reports whose labels hold an object
     it cannot serialise (today: a `Path`).
   - `pytest-xdist` installed into the root `.venv` (it is already in the dev group of `platform/pyproject.toml`).
   - The parallel command in `platform/AGENTS.md`.
2. **Faster state snapshot.** In `platform/tests/test_coq9_curation.py`, `semantics()` reads each constraint's
   coefficients straight from the solver's matrix under Gurobi and GLPK and keeps the symbolic path for any other
   solver. Term keys become variable names instead of the printed variable (name with bounds); the bounds are still
   compared in the `variables` entry.

Not in scope:

- change 3 (GLPK-backed models in the test files that solve no LP);
- any test's assertions, and the 44 test files that import the snapshot helpers;
- the four tests that fail only in parallel runs (two of the three files test lipid-line tools);
- the same snapshot pattern in builder code (`platform/scripts/gem_annotate/energy_candidates.py:124`);
- the temporary-folder leak in `conftest.py`, the rounding of Gurobi-backed `model.copy()` in tools, and other
  sessions' run copies in `/private/tmp/claude-501`: reported only.

No commit or push without a separate instruction.

## 1. Starting state (verified 2026-11-24)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/pytest-performance-analysis-3e0794` @ `2b428e29750ea598294f566f3b80c1dfe1f69a90` (PR #57 merge; equal to `origin/main`). Worktree clean |
| Python environment | root `.venv`: Python 3.13.5, cobra 0.30.0, optlang 1.8.3, gurobipy 13.0.1, swiglpk 5.0.13, scipy 1.18.0, pytest 9.0.2. `pytest-xdist` and `execnet` absent. Default solver Gurobi |
| Files to change | `platform/tests/conftest.py` `bd96427d…`, `platform/tests/test_coq9_curation.py` `167105c8…`, `platform/AGENTS.md` `ac219bb0…` |
| "Before" run, serial | Full suite on a short-path `git archive` copy of HEAD, root `.venv`, with a timing plugin that lives outside the repository: 501 passed, 3 failed, 24 skipped, 18 errors, 1391 sub-tests passed, 2137.84 s. `junit_full_serial.xml` `dad51f6c…` |
| "Before" run, parallel | Same copy, `-n 6 --dist loadfile`, the `pytest-performance-1daa1e` worktree's `.venv` (it has `pytest-xdist` 3.8.0): 759.03 s. Four tests fail only there. `junit_parallel_n6.xml` `b6bcdeef…` |

The measurements behind the proposal (where the time goes, by operation and by test kind) are in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/pytest_runtime_20261124/`, with the scripts that produced them.

## 2. Acceptance criteria

| # | Criterion |
| --- | --- |
| A1 | Full suite, serial, after against before: every test has the same outcome and the sub-test count is the same |
| A2 | Full suite with `-n 6 --dist loadfile` in the root `.venv`: every test has the same outcome as the serial run, including the four tests that failed only in parallel |
| A3 | On every model file under `model/`, under Gurobi and under GLPK, the rows read from the solver matrix equal the rows from the symbolic path: same constraints, same variable names, same coefficients |
| A4 | For a fixed list of edits to a model (cobra-level and solver-only), the new snapshot reports a difference exactly when the old one does, apart from the intended one: a bound that changes only between `0` and `0.0` |
| A5 | Measured times for the serial and the parallel run are reported next to the "before" times |
| A6 | A fresh read-only reviewer has checked the diff; findings are fixed or the disagreement is reported |

## 3. Stop conditions

- Any test outcome that differs between before and after, or between serial and parallel, and is not explained by the
  change: stop and report, no further edit.
- A3 or A4 fails on any model or edit: stop; the snapshot change is withdrawn.
- `pytest-xdist` cannot be installed at the version pinned in `platform/uv.lock` (3.8.0, with `execnet` 2.1.2): stop
  and report.

### Note added during the round (2026-11-24)

The second condition was not applied as written. A3 failed on the first version of the reader: on the five
lipid-unlump candidates it wrote a constraint without terms as an empty list, where the old snapshot writes the
constant 0 (13 rows per model). The difference was a defect of the new code, so the reader was corrected to write
such rows as before and A3 was rerun unchanged on all 61 models: 0 differences. The condition was meant for a
difference that cannot be removed. Details in [REPORT.md](REPORT.md), section 5.

## 4. Review follow-ups (authorized 2026-11-24, item 3 of section 0)

### Scope

1. **A test for the readers**, in `platform/tests/test_coq9_curation.py`. On the fixture, under Gurobi and under
   GLPK, with a one-sided and a two-sided extra constraint and an unused metabolite: the rows read from the solver's
   matrix equal the symbolic rows, and a solver-only coefficient edit changes the snapshot.
2. **Numeric bounds in each term.** Each term of the `objective` and `constraints` entries carries its variable's
   lower and upper bound as numbers. An assertion that compares one entry alone then sees a bound change of that
   entry's variables again, and `0` against `0.0` stays equal. The new test also covers this.
3. **The parallel command in the other documents:** `.agents/skills/curation-change/SKILL.md`, `README.md` and
   `platform/README.md`.

Not in scope, as before: change 3, any assertion of the 44 importing test files, their 30 outdated comments, and
the builder's own snapshot (`platform/scripts/gem_annotate/energy_candidates.py`).

### Starting state (verified 2026-11-24 19:36)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/pytest-performance-analysis-3e0794` @ `2b428e29`, with the uncommitted changes 1 and 2 |
| Files to change | `platform/tests/test_coq9_curation.py` `943c61d5…`, `.agents/skills/curation-change/SKILL.md` `ee3530bf…`, `README.md` `3837e44e…`, `platform/README.md` `7d6a9709…` |
| Runs to compare with | Serial, before: `junit_full_serial.xml` `dad51f6c…`. Six workers, after changes 1 and 2: `round2/junit_after_parallel.xml` `37681336…`. Both hold the same outcome for each of the 546 tests |

### Acceptance criteria

| # | Criterion |
| --- | --- |
| B1 | On every model file under `model/` and the fixture, under Gurobi and under GLPK: each term of the new snapshot is the old snapshot's term, with the same variable, the same coefficient, and bounds that print exactly as the old key printed them |
| B2 | For a fixed list of edits, the entries the new snapshot reports as changed (the objective, each constraint row, each other entry) are the entries the old snapshot reports. The one intended exception: an edit that only reprints an unchanged bound (`0` against `0.0`) is reported by the old snapshot alone |
| B3 | The new test passes, and it fails when either reader is broken on purpose |
| B4 | Full suite with `-n 6 --dist loadfile` in the root `.venv`: each of the 546 earlier tests has its earlier outcome, and the new tests pass |
| B5 | A fresh read-only reviewer has checked the follow-up diff; findings are fixed or the disagreement is reported |
| B6 | The time of one `semantics()` call is reported next to the earlier figure |

### Stop conditions

- B1 or B2 shows a difference that remains after a defect of the new code is corrected: stop and report; the
  bounds are withdrawn from the terms.
- Any earlier test with a different outcome that the change does not explain: stop and report, no further edit.

### Note added during the follow-ups (2026-11-24)

B2 is not met as written in 6 of the 120 cases (30 edits, two models, two solvers). The first stop condition was
not applied as written: the helper was not edited further, and the bounds were not withdrawn, because neither
difference comes from the bounds. The helper without them gives the same six results (the `round2` column of the
check). Whether the bounds stay is the user's decision.

- **A variable's type set to "integer"** (4 cases: both models, both solvers). optlang rounds the bounds of an
  integer variable to ints, so `0 <= R305 <= 1000.0` becomes `0 <= R305 <= 1000`. The old snapshot reports every
  row that holds the variable because of that reprint; the new one reports the change in `variables` only. This is
  the intended exception (a bound that changes only in print) inside an edit that also changes something real. The
  criterion stated the exception per edit; it applies per entry.
- **The objective reaction removed, under Gurobi, on the fixture** (2 cases). The new snapshot reports `objective`.
  The old one does not on its first call, because it reads the objective before optlang applies the queued
  removal. REPORT.md names this difference under change 2; it does not come from the follow-ups.

In the other 114 cases the reported entries are as the criterion says. Details in [REPORT.md](REPORT.md),
section 10.

## 5. The variable's type in each term (authorized 2026-11-24, after section 4)

The report offered three ways to settle B2: keep the snapshot as it was, also carry the variable's type in each
term, or withdraw the bounds. Asked "Keep the bounds in the snapshot terms as they are now?", the user answered:
"Add the type too".

### Scope

In `platform/tests/test_coq9_curation.py` only. Each term of the `objective` and `constraints` entries also carries
its variable's type, so a term holds what the `variables` entry holds for that variable: lower bound, upper bound
and type. The two tests of section 4 are extended to cover it. Nothing else changes.

### Starting state (verified 2026-11-24 21:26)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/pytest-performance-analysis-3e0794` @ `2b428e29`, with the uncommitted work of sections 0 to 4 |
| File to change | `platform/tests/test_coq9_curation.py` `03ceff2a…` |
| Run to compare with | Six workers, after the follow-ups: `round3/junit_parallel.xml` `714ab723…` (548 tests) |

### Acceptance criteria

| # | Criterion |
| --- | --- |
| C1 | On every model file under `model/` and the fixture, under Gurobi and under GLPK: each term is the old snapshot's term, as in B1, and its bounds and type are those of the `variables` entry |
| C2 | For the 30 edits of B2, the new snapshot reports the entries the old one reports, with two exceptions known beforehand: an edit that only reprints an unchanged bound is reported by the old snapshot alone; and under Gurobi the old snapshot's first call after the objective reaction is removed does not report `objective`. For one added edit, a type change of a variable whose bounds are already ints, the new snapshot reports `variables` and every entry that holds the variable, where the old one reports `variables` only |
| C3 | The two tests pass, and fail when the type is left out of the terms or a reader is broken |
| C4 | Full suite with `-n 6 --dist loadfile` in the root `.venv`: each of the 548 earlier tests has its earlier outcome |
| C5 | A fresh read-only reviewer has checked the diff; findings are fixed or the disagreement is reported |
| C6 | The time of one `semantics()` call is reported next to the earlier figures |

### Stop condition

Any difference in C1 or C2 other than those named, or any earlier test with a different outcome: stop and report,
with no further edit.
