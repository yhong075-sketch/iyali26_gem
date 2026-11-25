**Status: changes 1 and 2, the three review follow-ups and the variable's type in each term are in the working
tree. The criteria for changes 1 and 2 (A1 to A6) and for the type (C1 to C6) are met. Of the six for the
follow-ups, B2 was not met as written (section 10); the user settled it by choosing the type (section 11). No
decision is open. Committed on branch `claude/pytest-performance-analysis-3e0794`; not pushed.**

# Shorter pytest runs: results, 2026-11-24

Scope and authorization: [TASK.md](TASK.md). Branch `claude/pytest-performance-analysis-3e0794`: the changes below are
one commit on top of `2b428e29`. No model, curation file or builder step changed, and no existing test
assertion. Two tests were added (sections 10 and 11).

## 1. Result

Full suite of `main` (`2b428e29`), root `.venv`, Gurobi as default solver, on short-path copies:

| Run | Time | Outcomes |
| --- | --- | --- |
| Serial, before | 35.6 min (2,138 s) | 501 passed, 3 failed, 24 skipped, 18 errors, 1,391 sub-tests |
| Serial, after changes 1 and 2 | 31.6 min (1,899 s) | the same for each of the 546 tests |
| Six workers, after changes 1 and 2 | 12.1 min (727 s) | the same for each of the 546 tests |
| Six workers, after the follow-ups | 6.6 min (398 s) | the same for each of the 546 tests; the 2 new tests pass (503 passed, 1,395 sub-tests) |
| Six workers, after the type was added | 8.1 min (489 s) | the same for each of the 548 tests |

- The three six-worker runs differ because of the machine, not the code. The first shared it with three jobs of
  other sessions (load average 5.2 at its start; its test times add up to 4,116 s). The second had it alone (load
  average 1.5 at its start; 2,236 s). The third shared it with one serial pytest of another session (2,706 s).
- Against the serial run before the change, six workers take 34 % of the time on the shared machine and 19 % on
  the idle one. No serial run was made on an idle machine after the change.
- The four tests that failed only in parallel runs now pass there.
- The 3 failures and 18 errors are the standing ones (missing optional packages, files outside the repository,
  recorded hashes). None involves the changed code.
- Timing conditions differ between the runs, so the times are single observations, not a controlled comparison:
  the "before" run and the last run were on mains power; the three runs in between were on battery. During the
  two runs after changes 1 and 2, other sessions ran a serial pytest, an FVA and a knockout screen on the same
  machine.

## 2. Why the suite is slow (measured before the change)

A timing plugin outside the repository wrapped the main operations during the serial "before" run.

| Operation | Time | Share | Calls | Per call |
| --- | --- | --- | --- | --- |
| `model.copy()` | 983 s | 46.0 % | 937 | 1.05 s |
| `read_sbml_model` | 613 s | 28.7 % | 251 | 2.44 s |
| State snapshots (`semantics`, `annotations`) | 391 s | 18.3 % | 1,297 | 0.30 s |
| `write_sbml_model` | 83 s | 3.9 % | 59 | 1.41 s |
| Code under test (curation steps, build chain, LP solves) | 12 s | 0.6 % | | |
| Everything else | 53 s | 2.5 % | | |

- **Copies.** The guard tests copy the whole model once per sub-case. cobra deep-copies the solver, and optlang's
  Gurobi interface does that by writing the problem as LP text and rebuilding every variable and constraint
  (`optlang/gurobi_interface.py`, `Model.__getstate__` / `__setstate__`). 924 of the 937 copies are Gurobi-backed.
- **Reads.** Tests read the 14 MB SBML again instead of copying, because a Gurobi-backed copy is not exact: the LP
  text keeps about 6 significant digits (72 of 18,084 coefficients change on candidate `2f352037`).
- **Snapshots.** `semantics()` built one symbolic expression for each of about 1,900 constraints on every call.
- **Growth.** Each curation step adds a test file with the same template. By kind of test, before the change:
  stale-input guards 784 s (40 tests), round trips 594 s (53), scope guards 295 s (35), all others 463 s (418).

## 3. Change 1: test files in parallel

**What changed.**
- [conftest.py](../../../platform/tests/conftest.py): a `pytest_report_to_serializable` hook and its helper
  `_sendable`.
- Root `.venv`: `pytest-xdist` 3.8.0 and `execnet` 2.1.2, the versions pinned in `platform/uv.lock`, installed from
  the local uv cache (`uv pip install --offline`). Nothing else in the environment changed.
- [platform/AGENTS.md](../../../platform/AGENTS.md): the hand-back run is now
  `python -m pytest -n 6 --dist loadfile`.

**Logic.** pytest-xdist sends each test report from a worker to the controller with execnet, which accepts plain
values only. pytest puts a sub-test's labels into its report unchanged, so `subTest(conflict=some_path)` made the
whole test fail with `DumpError: can't serialize PosixPath`. The hook runs before pytest's own serialising and
replaces a label that execnet cannot send by its `repr`. Labels execnet can send are left as they are. In a serial
run nothing calls the hook: among the installed packages only `xdist/remote.py` does.

**Why a hook and not edits to the four tests.** Two of the three files test lipid-line tools, which this line may
only review. The hook also covers any later test that passes an object as a label.

**`--dist loadfile`** keeps each test file on one worker, because a file's tests share one model read in
`setUpClass` and one module-scoped solver fixture.

## 4. Change 2: the snapshot reads solver rows directly

**What changed.** In [test_coq9_curation.py](../../../platform/tests/test_coq9_curation.py):
- new `constraint_terms(model)`: one pass over the solver's matrix, through `getRow` for Gurobi and
  `glp_get_mat_row` for GLPK; any other interface keeps the symbolic path;
- `semantics()` takes its constraint rows from it;
- `expression_terms()` keys each term by the variable's name. Later rounds added the variable's two bounds
  (section 10) and its type (section 11) to each term.

The 44 test files that import `semantics` are untouched.

**Effect on time** (same 1,297 calls before and after):

| | Before | After |
| --- | --- | --- |
| `semantics()`, total | 308.8 s | 63.9 s |
| `semantics()`, per call | 0.238 s | 0.049 s |
| Test time of the serial run | 2,134 s | 1,897 s (−11.1 %) |

`annotations()` also fell, from 81.6 s to 32.8 s, with unchanged code and call count. The cause was not
investigated.

**What the snapshot compares now** (after sections 10 and 11). It carries the content of the old snapshot with
two differences:

1. **A term holds the variable's own values where it held its print.** A term was keyed by the printed variable
   (`0 <= R305 <= 1000.0`). It is now `("R305", coefficient, 0, 1000.0, "continuous")`: the variable's name, the
   coefficient, and what the `variables` entry holds for the variable (lower bound, upper bound, type).
   - Bounds that are equal as numbers but print differently are no longer a difference: `0`, `0.0` and `-0.0`;
     `1000` and `1000.0`.
   - So under GLPK a `model.copy()` has the same snapshot as its original. Under Gurobi it still differs, because
     of the rounded coefficients.
   - An entry compared alone sees a change of its variables' bounds, as it did before this round. Between change
     2 and the follow-ups a term held the name alone, and three assertions had lost that (section 6, finding 2).
   - An entry compared alone also sees a change of a variable's type. The old snapshot saw one only when optlang's
     rounding of the bounds reprinted them (`1000.0` as `1000`), and not when the bounds were ints already
     (sections 10 and 11).
2. **Queued changes are applied first.** optlang queues added and removed variables and constraints until its next
   `update()`. The old snapshot applied the queue part-way through, after reading the objective. The observed
   effect: under Gurobi, right after the objective reaction is removed, the old snapshot's first call still lists
   the removed variables in `objective`; the new one does not.

Unchanged: a constraint without terms is written as the symbolic path writes it, as the constant 0. Also
unchanged, and a limit of both snapshots: under GLPK, once a snapshot has been taken, later snapshots keep listing
a removed objective reaction in `objective`. optlang's GLPK interface rebuilds its stored objective expression
only after `set_linear_coefficients`, not after a variable is removed (`optlang/glpk_interface.py`, lines 308-320
and 732-748 in optlang 1.8.3).

## 5. Verification of changes 1 and 2

| # | Criterion | Result |
| --- | --- | --- |
| A1 | Serial, after against before | 0 of 546 tests differ; 1,391 sub-tests in both |
| A2 | Six workers against serial | 0 of 546 tests differ; the four parallel-only failures pass |
| A3 | Same rows as the symbolic path on every model | 0 of 61 models differ under Gurobi and under GLPK |
| A4 | Same change detection as before | 108 of 108 cases as expected |
| A5 | Times reported | Section 1 |
| A6 | Independent review | No bug; five findings (section 6) |

**A3 in detail.** 61 models: the 60 SBML files under `model/` and the test fixture (the other 21 XML files under
`model/` are not SBML). Per solver, 118,359 constraints and 1,105,653 terms were compared. The key of each old
term was translated to the variable's name and the terms then had to be equal as they are. The symbolic fallback
was checked on the fixture under the `glpk_exact` and `scipy` interfaces.

**A4 in detail.** 27 kinds of edit, on the fixture and on candidate `2f352037`, under both solvers:

| Kind of edit | Cases | Old snapshot reports it | New snapshot reports it |
| --- | --- | --- | --- |
| Changes the mathematics (coefficients, bounds, constraints, objective, species, gene rule) | 84 | yes | yes |
| Changes nothing the snapshot covers | 16 | no | no |
| Reprints an unchanged bound (`0` ↔ `0.0`) | 8 | yes | no |

After every edit the two snapshots still carried the same content.

**Two defects in my first version, both found before any full run and fixed.**
- It read the solver before optlang's queued additions were applied, so a just-added constraint raised `KeyError`.
  The helper's own test `test_metadata_preserves_full_mathematics_and_extra_constraint` caught it.
- It wrote a constraint without terms as an empty list, where the old snapshot writes the constant 0. A3 caught it
  on the five lipid-unlump candidates (13 unused metabolites each). TASK.md says to withdraw the change if A3
  fails. I fixed the reader and reran A3 unchanged instead, because the difference was a defect of the new code
  and not a property of the approach; the note is added to TASK.md.

## 6. Independent review of changes 1 and 2

A fresh read-only reviewer (the repository's `code-reviewer` instructions) checked the diff and ran its own
in-memory checks. It found no correctness bug and no rule violation. Its five findings (the review of the
follow-ups is in section 10):

| # | Finding | Handling |
| --- | --- | --- |
| 1 | The two new readers have no test in the repository. A reader that returns no terms at all still passes the helper's own tests. With Gurobi as default, the suite never runs the GLPK reader or the fallback | Agreed. Tests added in the follow-ups (section 10) |
| 2 | Three assertions no longer see variable bounds: `test_reaction_selection.py:47` (`objective`), `:49` (`constraints["extra_cap"]`) and `test_gpi_chain_links.py:104` (`objective`) | Agreed; "not weaker" held per test, not per assertion. Restored in the follow-ups: each term carries its bounds (section 10) |
| 3 | Comments in importing test files give the print format as the reason for fresh reads | Agreed; by my count 30 comment lines in 28 files. Out of scope, not edited |
| 4 | `.agents/skills/curation-change/SKILL.md:32` still names the serial command | Agreed; `README.md:25` and `platform/README.md:14,60` did too. All three updated in the follow-ups (section 10) |
| 5 | The hand-off said "61 SBML files under `model/`"; it is 60 plus the fixture. Evidence was only in the scratchpad | Corrected here; evidence moved (section 9) |

For finding 2, as things stood after change 2 (since the follow-ups the three assertions see the change
themselves), the later lines that still caught the change, checked in the code:
- a reaction's bounds: `test_reaction_selection.py:52` checks every reaction of that test's model against the
  curated after-state, and `test_gpi_chain_links.py:108` compares the whole `reactions` entry;
- a bound changed only in the solver: the SBML round trip at `test_reaction_selection.py:65` and
  `test_gpi_chain_links.py:128` compares the whole snapshot, including `variables`. The reviewer replayed this; I
  did not rerun it.

## 7. Decisions and what is not done

**The three review follow-ups are done** (section 10): tests for the readers, bounds in each term, and the
parallel command in the `curation-change` skill and the two READMEs.

**No decision is open.** Criterion B2 of the follow-ups was not met as written in 6 of 120 cases (section 10), and
TASK.md said to withdraw the bounds in that case. The report offered three options:

1. keep the snapshot as it was (my recommendation then): an entry compared alone would not see a change of a
   variable's type alone;
2. also carry the variable's type in each term;
3. withdraw the bounds, which removes neither of the two differences.

The user chose the second: "Add the type too". It is done (section 11).

**Not done, by scope:**
- Change 3 (GLPK-backed models in the test files that solve no LP), estimated at −38 % of the serial time. Under
  GLPK the `objective` entry can stay stale after the objective reaction is removed (section 4); that needs a look
  before change 3.
- The 30 outdated comments.
- The same snapshot pattern in builder code (`platform/scripts/gem_annotate/energy_candidates.py:124`).
- `conftest.py` creates a temporary folder per run and never removes it.
- Gurobi-backed `model.copy()` outside tests: a solve on such a copy uses rounded coefficients. A text scan of
  `platform/scripts` and `platform/tools` finds 16 call sites in 11 files (a name containing `model`, then
  `.copy()`). None was examined for whether a solve follows.
- `platform/README.md`, "Test-order coupling": it says `test_vacuole_supply::test_xml_only_four_bounds…` fails
  after `test_lipid_moiety_ledger`. That test passed in every full run of this round, serial and parallel. The
  note was not edited.

**Limits of the evidence:**
- One run per condition; machine shared; power source differs (section 1).
- Only six workers were tried.
- The follow-ups and the type were run in full on six workers only, not serially.
- A variable named `1` in a model would make the snapshot's sort of terms raise (section 11, third review). No
  model file here has one.
- The other active line (`gene-essentiality-reactions`, 59 to 74 min serial) was not run with these changes.
- The hook has no test of its own. The four tests with `Path` labels fail in a parallel run without it.

## 8. Could it be faster, and how to test it

After the change, the serial run splits into copies 54 %, reads 33 %, snapshots 5 %, writes 5 %.

| Idea | Evidence so far | Test |
| --- | --- | --- |
| GLPK-backed models where no LP is solved (change 3) | A GLPK copy takes 0.12 s against 0.9 s and is exact | Trial run with GLPK as default and solving blocked; compare test by test |
| symengine for optlang | optlang here runs on sympy (`USE_SYMENGINE` is False, symengine not installed); 55 % of a copy's own time is inside sympy. How much symengine saves is untested | Scratch environment with symengine; repeat the operation benchmark and the full suite |
| One copy shared by a scope test's sub-cases | Those cases never reach a model edit | Per file; estimated 245 s in total today |
| GLPK copies instead of fresh reads | Exact under GLPK | Depends on change 3 |
| More workers | Not tried; the machine has 4 performance and 6 efficiency cores and 16 GB | Runs with 8 and 10 workers, watching memory |

## 9. Reproduce and evidence

From `platform/` of a short-path copy (`git archive HEAD` plus the changed files), with the root `.venv`:

```bash
python -m pytest -n 6 --dist loadfile --continue-on-collection-errors -p no:cacheprovider -q --junitxml=after_parallel.xml
```

```bash
python -m pytest --continue-on-collection-errors -p no:cacheprovider -q --junitxml=after_serial.xml
```

Evidence, in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/pytest_runtime_20261124/`:
- the "before" measurements and their scripts (junit files, timing data, benchmarks);
- `round2/`: the junit files, logs and timing data of both runs after changes 1 and 2, the equivalence check
  (`check_snapshot.py`, `check_gurobi.json`, `check_glpk.json`), the queue check, the launcher, and
  `SHA256SUMS.txt`. `debug_differences.txt` shows the empty-row defect of the first version;
- `round3/`: everything for the follow-ups (section 10).

Code identity:
- runs after changes 1 and 2: `conftest.py` `3a24b8f0…`, `test_coq9_curation.py` `943c61d5…`;
- run after the follow-ups: `conftest.py` `3a24b8f0…`, `test_coq9_curation.py` `03ceff2a…`;
- run after the type was added: `conftest.py` `3a24b8f0…`, `test_coq9_curation.py` `d42ce045…`, equal to the
  working tree.

## 10. Review follow-ups (2026-11-24, same day)

Authorized by the user's "All three" ([TASK.md](TASK.md), section 4). This section describes the state after the
follow-ups. Section 11 then added the variable's type to each term, which settles the open point of B2 below.

### What changed

- **Bounds in each term** ([test_coq9_curation.py](../../../platform/tests/test_coq9_curation.py)).
  `expression_terms()` writes a term as (variable name, coefficient, lower bound, upper bound). The two readers in
  `constraint_terms()` attach the same bounds, taken from the model's variables in one pass. A constant term has
  no bounds.
- **Two tests for the snapshot**, in the same file, class `SnapshotTests`. Both use the 8-reaction fixture under
  Gurobi and under GLPK, with a one-sided and a two-sided extra constraint and an unused metabolite.
  - `test_matrix_rows_equal_the_symbolic_rows`: the rows read from the matrix equal the symbolic rows, also while
    a row is still queued in optlang; and a snapshot builds no symbolic row (the test fails if it does).
  - `test_solver_only_edits_show_in_every_entry_that_holds_them`: a coefficient changed in the solver alone shows
    in its row and nowhere else; a bound changed in the solver alone shows in `variables`, in `objective` and in
    every row that holds the variable; the same bounds written in the other numeric type are no change.
- **The parallel command** in `.agents/skills/curation-change/SKILL.md`, `README.md` and `platform/README.md`.
  `platform/README.md` also says what the two options do and that they need `pytest-xdist`.

The 44 test files that import `semantics` are untouched. None of them reads inside a term (scan of all 44).

### Verification

| # | Criterion | Result |
| --- | --- | --- |
| B1 | Each new term is the old term, on every model | Met: 0 of 61 models differ, under Gurobi and under GLPK |
| B2 | The same entries reported as changed, edit by edit | **Not met as written in 6 of 120 cases** (below) |
| B3 | The new tests pass, and fail on a broken reader | Met: 13 of 13 deliberate breakages fail |
| B4 | Full suite on six workers | Met: 0 of 546 earlier tests differ; the 2 new tests pass |
| B5 | Independent review | Done: no bug; two test findings (below) |
| B6 | Time of one snapshot | 25.7 ms, against 22.3 ms before the follow-ups and 139.5 ms at HEAD (Gurobi) |

**B1 in detail.** The same 61 models as in A3. Per solver: 118,359 constraints, 1,105,653 terms and 122 objective
terms. Each new term was printed the way optlang prints a variable (`lb <= name <= ub`) and had to give the old
key exactly, with an equal coefficient. The 65 terms without bounds are the constant terms of the 65 rows without
variables. The symbolic fallback was checked on the fixture under `glpk_exact` and `scipy`: no difference.

**B2 in detail.** 30 kinds of edit: the 27 of A4, and three that put the edited variable into the objective and
into an extra constraint. On the fixture and on candidate `2f352037`, under both solvers, that is 120 cases. For
each case the entries that differ before and after the edit were listed for the old and for the new snapshot:
`objective`, each constraint row, and each other entry.

| Kind of edit | Cases | Result |
| --- | --- | --- |
| Changes the mathematics | 92 | 86: the same entries in both. 6: see below |
| Changes nothing the snapshot covers | 16 | none reported by either |
| Only reprints an unchanged bound | 12 | reported by the old snapshot alone, as intended |

The six cases:

- **A variable's type set to "integer"** (4 cases: both models, both solvers). Both snapshots report `variables`.
  The old one also reports the six rows that hold the variable (and `objective` on the fixture); the new one does
  not.
  - Cause: optlang's type setter rounds the bounds of an integer variable to ints (`optlang/interface.py`,
    `Variable.type`), so `0 <= R305 <= 1000.0` becomes `0 <= R305 <= 1000`. In those rows nothing changed but
    that print. It is the intended exception, inside an edit that also changes something real.
  - Consequence at this stage: an assertion that compares `objective` or one row alone does not see a change of
    type alone. A whole-snapshot comparison does. Section 11 changes this.
  - No code under `platform/scripts`, `platform/tools` or `platform/tests` sets a variable's type (text scan, 0
    matches).
- **The objective reaction removed, under Gurobi, on the fixture** (2 cases: "remove a reaction" and "remove
  every reaction of a metabolite"). The new snapshot reports `objective`; the old one does not on its first call.
  This is difference 2 of section 4. It dates from change 2.

Neither comes from the bounds. The snapshot without them, as it was after change 2, gives the same result in all
six cases. Against the old snapshot it reports fewer entries in 20 of the 92 cases (the bound edits and the type
edit); the new one does so in 4 (the type edit).

TASK.md said to withdraw the bounds if B2 failed. I did not, for the reason just given, and left the decision to
the user, who chose to add the type as well (section 11).

B1 and B2 were run twice: before the second review's edits to the tests (file `961461c9…`) and on the final file
(`03ceff2a…`). Apart from timings the two results are identical.

**The three assertions, replayed.** On the models of their own tests, with one bound changed in the solver alone,
the compared entry differs with the old snapshot, did not differ after change 2, and differs now. This holds for
`test_reaction_selection.py:47` and `:49` and for `test_gpi_chain_links.py:104`, under both solvers.

**B3 in detail.** 13 deliberate breakages of the helper, in a scratch copy:
- a reader that returns no terms, keeps the helper column of a two-sided constraint, doubles the coefficients or
  drops the last term of each row;
- the queue not applied, an empty row written differently, the bounds swapped or left out;
- three that bypass a reader while the result stays right.

Each fails at least one of the two new tests. 12 of the 13 pass all nine earlier tests of the file, which is the
gap the first review named.

**B4 in detail.** `python -m pytest -n 6 --dist loadfile` on a short-path copy with the root `.venv`: 398 s; 503
passed, 3 failed, 24 skipped, 18 errors, 1,395 sub-tests. Compared test by test with the serial run before any
change and with the six-worker run after changes 1 and 2: no outcome differs among the 546 shared tests.

**B6 in detail.** On candidate `2f352037` (1,883 rows, 4,692 columns), best of nine calls:

| Solver | At HEAD | Before the follow-ups | Now |
| --- | --- | --- | --- |
| Gurobi | 139.5 ms | 22.3 ms | 25.7 ms |
| GLPK | 134.6 ms | 12.9 ms | 16.8 ms |

The bounds cost about 4 ms per call.

### Second independent review

A fresh read-only reviewer checked the follow-up diff and ran its own in-memory checks: 27 model states under
both solvers, 4 full models and 9 breakages. It found no bug, no rule violation and nothing out of scope. Two
findings on the tests:

| # | Finding | Handling |
| --- | --- | --- |
| 1 | B2 is not met as written (the six cases above). They need a decision, not a silent pass | Agreed. Reported here and in TASK.md; decided by the user (section 11) |
| 2 | Neither new test could tell whether the matrix readers are used at all: both passed with the rows built symbolically | Fixed. The first test now fails if a constraint's symbolic expression is built during a snapshot; three breakages of that kind were added to B3 |

Its notes, none a defect:

- Under Gurobi, a two-sided constraint renamed after it was added makes the snapshot raise `KeyError`. optlang's
  own `constraint.expression` raises the same error, so the old snapshot did too; after change 2 the helper column
  came back as a term instead. No test renames such a constraint. The reviewer ran this; I did not.
- The test's value for an empty row was the sympy form. It is now taken from the symbolic library itself.
- "About a third" of the serial time, as `platform/AGENTS.md` and `platform/README.md` say, was 38 % on the shared
  machine (727 s against 1,899 s). The later run on an idle machine took 398 s. The wording was kept.
- `platform/README.md` still describes a test-order failure that no run of this round shows (section 7).

### Evidence

`$IYALI26_RESEARCH_ROOT/artifacts/tasks/pytest_runtime_20261124/round3/`, with `SHA256SUMS.txt`:

| File | Content |
| --- | --- |
| `check_snapshot.py`, `check_gurobi*.json/.txt`, `check_glpk*.json/.txt` | B1 and B2 |
| `explain_b2_cases.py/.txt`, `explain_objective_removal.py/.txt` | The six cases of B2, shown on the fixture |
| `replay_single_entry.py/.txt` | The three assertions, replayed |
| `mutation_check.py/.json/.txt`, `mutation_check_before_review.*` | B3, after and before the second review |
| `junit_parallel.xml`, `pytest_parallel.txt`, `parallel_started.txt`, `run_full.zsh` | B4 |
| `bench_snapshot.py/.txt` | B6 |
| `check_queued.py/.txt` | The queue check of round 2, repeated: the same lines |
| `head_test_coq9_curation.py`, `round2_test_coq9_curation.py` | The helper at HEAD and after change 2 |
| `reviewer/` | The second reviewer's own scripts |

## 11. The variable's type in each term (2026-11-24, same day)

Authorized by the user's "Add the type too" ([TASK.md](TASK.md), section 5).

### What changed

In [test_coq9_curation.py](../../../platform/tests/test_coq9_curation.py) only:

- A term is (variable name, coefficient, lower bound, upper bound, type), in `expression_terms()` and in both
  readers. A constant term has none of the three.
- The two tests cover it:
  - the pinned terms include the type;
  - a type changed in the solver alone must show in `variables`, in `objective` and in every row that holds the
    variable;
  - with an integer variable in the model, the objective's term is pinned and the rows read from the matrix must
    still equal the symbolic rows.

### Verification

| # | Criterion | Result |
| --- | --- | --- |
| C1 | Each new term is the old term, with the bounds and type of `variables` | Met: 0 of 61 models differ, under Gurobi and under GLPK |
| C2 | The same entries reported as changed, with the exceptions named beforehand | Met: 124 of 124 cases |
| C3 | The tests fail when the type is left out or a reader is broken | Met: 20 of 20 deliberate breakages fail |
| C4 | Full suite on six workers | Met: 0 of 548 earlier tests differ |
| C5 | Independent review | Done: no bug that can occur here; three low findings (below) |
| C6 | Time of one snapshot | 27.8 ms, against 27.0 ms with the bounds alone (Gurobi) |

**C1 in detail.** The same 61 models. Per solver: 118,359 constraints, 1,105,653 terms and 122 objective terms.
As in B1, each term printed the way optlang prints a variable gives the old key. In addition, each term's bounds
and type equal the snapshot's own `variables` entry. Every variable in these models is continuous: 1,105,588 terms
say "continuous", and the other 65 are the constant terms of the rows without variables. So on the model files C1
can only confirm the value "continuous". Other types are covered by C2 and C3.

**C2 in detail.** The 30 edits of B2 and one more, a type change of a variable whose bounds are ints already. On
the two models under both solvers that is 124 cases.

| Kind of edit | Cases | Result |
| --- | --- | --- |
| Changes the mathematics | 92 | The same entries in both. In 2 of them the new snapshot also reports `objective` (the Gurobi exception named beforehand) |
| Changes nothing the snapshot covers | 16 | None reported by either |
| Only reprints an unchanged bound | 12 | Reported by the old snapshot alone (named beforehand) |
| Type change that reprints nothing | 4 | The old snapshot reports `variables` only; the new one also every entry that holds the variable (named beforehand) |

The four type-change cases that B2 did not meet now give the same entries in both snapshots.

**The three assertions, replayed again.** A variable of the compared entry was changed in the solver alone, on the
models of the assertions' own tests. The result is the same for all three assertions and both solvers:

| Change | Old snapshot | After change 2 | With bounds | With bounds and type |
| --- | --- | --- | --- | --- |
| Upper bound | seen | not seen | seen | seen |
| Type, with the bounds reprinted by optlang | seen | not seen | not seen | seen |
| Type, with bounds that are ints already | not seen | not seen | not seen | seen |

**C3 in detail.** The 13 breakages of B3; four that leave the type out or write "continuous" for every variable;
and three that write an integer variable as binary. Each fails at least one of the two tests. 19 of the 20 pass
all nine earlier tests of the file.

**C4 in detail.** 489 s on six workers: 503 passed, 3 failed, 24 skipped, 18 errors, 1,395 sub-tests. Compared
test by test with the run after the follow-ups (548 tests) and with the serial run before any change (546 tests):
no outcome differs. Another session's serial pytest ran on the machine throughout.

**C6 in detail.** Candidate `2f352037`, best of nine calls, all four versions measured in one run while another
session's pytest was running:

| Solver | At HEAD | After change 2 | With bounds | With bounds and type |
| --- | --- | --- | --- | --- |
| Gurobi | 148.4 ms | 22.8 ms | 27.0 ms | 27.8 ms |
| GLPK | 144.9 ms | 14.0 ms | 16.9 ms | 17.8 ms |

The type costs under 1 ms per call.

C1 and C2 were run on the file `54a3bce7…`. The final file `d42ce045…` differs from it only by two assertions
added after the review; the helper's code, everything above the first test class, is byte-identical.

### Third independent review

A fresh read-only reviewer checked this change: all 61 models under both solvers, 78 runs on constructed model
states, and its own breakages. It found no bug that can occur in this repository. Three findings, all low:

| # | Finding | Handling |
| --- | --- | --- |
| 1 | The tests pinned a term's content for continuous variables only. After the type change they checked which entries changed, not what the terms say. The empty row is compared with the helper's own output | Partly fixed. After the type change the tests now pin the objective's term and compare the matrix rows with the symbolic rows again; three breakages of that kind were added to C3. Not adopted: a literal for the empty row. The second review had it replaced, because the literal is sympy's form of an empty expression, and the snapshot only needs both paths to agree |
| 2 | REPORT.md, TASK.md and the STATE.md entry still described the state before this change | Updated |
| 3 | A variable named `1` without bounds, with the same coefficient as a constant in the same expression, makes the sort of the terms raise `TypeError` | Not changed. The reviewer checked that it cannot occur here: no model file has a reaction id `1`, cobra gives every reaction variable a lower bound, and no test adds a variable. With a bounded variable `1` the version of section 10 already raised |

Other observations of the reviewer, which I did not rerun:

- Under GLPK, an objective that still lists a removed variable (section 4) carries that variable's own type.
- The symbolic path for other interfaces runs in no test of the repository. `glpk_exact` and `scipy` refuse
  integer and binary variables, so the type is always "continuous" there.
- Apart from the new test, none of the repository's 701 Python files sets a variable's type or creates an integer
  or binary variable.

### Evidence

`$IYALI26_RESEARCH_ROOT/artifacts/tasks/pytest_runtime_20261124/round4/`, with `SHA256SUMS.txt`:

| File | Content |
| --- | --- |
| `check_snapshot.py`, `check_gurobi.json/.txt`, `check_glpk.json/.txt` | C1 and C2 |
| `replay_single_entry.py/.txt` | The three assertions, replayed |
| `mutation_check.py/.json/.txt`, `mutation_check_before_review.*` | C3, after and before the third review |
| `junit_parallel.xml`, `pytest_parallel.txt`, `parallel_started.txt`, `run_full.zsh` | C4 |
| `bench_snapshot.py/.txt`, `bench_snapshot_loaded.txt` | C6, and an earlier run of it under load |
| `check_queued.py/.txt` | The queue check of round 2, repeated: the same lines |
| `head_`, `round2_`, `round3_`, `checked_test_coq9_curation.py` | The helper at HEAD, after change 2, after the follow-ups, and as checked by C1 and C2 |
| `reviewer/` | The third reviewer's own scripts |
