# Independent review of the screen script (2026-10-30)

Reviewer: general-purpose agent following `.claude/agents/code-reviewer.md`. Static and read-only: no solves, no
tests. It read the two scripts, TASK.md, `ko_compare_two.py`, `checked_solve.py`, `essentiality_simulation_context.py`
and the cobra 0.30.0 / optlang sources in /opt/miniconda3. Result: no high-severity finding, 4 low, 1 note.

## Verified correct (reviewer)

- **Rule A matches the old scorer.**
  - The scoring code differs from `ko_compare_two.py` lines 84–121 only by the added `rule` and `excluded`
    parameters and the rule-B branch.
  - The sweep loop, the fresh-re-solve gate and the mtDNA set are the same. In the earlier candidate run all 20
    non-optimal statuses were `infeasible`, and there were no fresh re-solves. A2 is therefore a fair test.
- **Context handling.**
  - Nested `with model:` blocks unwind in order. The objective reset restores both the expression and the direction.
  - The `xMAINTENANCE` lower bound is read before it is set to 0.
  - `slim_optimize` then returns the maximum `xMAINTENANCE` flux.
  - The ceiling check uses the knockout's own closed reactions, and `forced()` is evaluated inside the knockout
    context.
- **Rule B is limited to certified knockouts.** It returns 0.0 only for an infeasible row whose certificate passed,
  so A4 holds by construction.
- **Solver failures cannot count as lethal.** optlang maps only Gurobi `INFEASIBLE` to `infeasible`.
  `INF_OR_UNBD`, `TIME_LIMIT`, `NUMERIC` and `SUBOPTIMAL` stay unresolved, and a non-optimal or NaN ceiling makes
  condition 3 false.
  - Limit: the ceiling maximum relies on Gurobi's optimality status. `check()` tests the primal point only.

## Findings and what was done

| # | Finding | Action |
| --- | --- | --- |
| 1 | `run_screen.zsh` used `$short` as given, so a relative path would break after `cd platform/` | Fixed: `short=${1:A}` |
| 2 | WT and WT-ceiling asserts ran only after the full sweep, so a failure discarded all outputs | Fixed: WT is checked before the loop (stop as declared); the WT ceiling is recorded, not asserted |
| 3 | Research root hard-coded; input hashes recorded but not gated | Fixed: `IYALI26_RESEARCH_ROOT` required (fail closed); medium, profile, list and `checked_solve.py` gated on the TASK §1 SHAs |
| 4 | Under `set -e`, `echo "sweep exit $?"` printed only 0 | Fixed: `set +e` around the Python step, real exit code echoed and returned |
| note | If the certificate re-solve came back optimal and passed, that growth is not used; the knockout stays unresolved under both rules | Kept: consistent with the rule text |
