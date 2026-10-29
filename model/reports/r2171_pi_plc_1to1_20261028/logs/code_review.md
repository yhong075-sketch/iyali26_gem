# Code review of the R2171 1:1 step (2026-10-28)

One round by a fresh general-purpose agent following `.claude/agents/code-reviewer.md`. It was given the task, the
files and the tests already run, not the reasoning, and was read-only (it ran the focused test and in-memory probes).

**Verdict:** no bugs in the shipped path. On `9d655c2b` with the shipped JSON the step changes only R2171's
stoichiometry and notes, fails closed before any edit on every model-side mismatch tried, and reruns (also after an
SBML save and reload) as `already_correct`. Two low findings, both reachable only by editing the curated JSON:

1. **`before_stoichiometry` not pinned.** With `combine=False`, a species added to the JSON's before-equation (with a
   matching species row and model) would survive the rewrite while the step reported `applied`; a missing
   `before_stoichiometry` raised `KeyError` instead of the scope error. **Fixed:** a `BEFORE` constant is pinned in
   code beside `AFTER`, the scope check requires both, and the species set must equal `AFTER`'s keys.
2. **`notes` not checked.** Empty notes would void the "after-equation without these notes" guard, and a note with
   `<`, `>` or `&` would be dropped by libSBML. **Fixed:** the scope check requires a non-empty dict of strings without
   those characters.

**Untested guards it listed, now tested:** a conflicting note on the after-state; the H2O formula; the DAG
`metanetx.chemical` identity; the deeper flag chain (all flags up to `--r1726-cpt-open` without `--dolichol-route`);
the two scope fixes (a changed before-equation, empty notes, a note containing `->`).

**Also reported by the reviewer (inferred, not a model check):** with the common DAG formula C5H6O5R2 the 1:1
equation with 2 H+ balances C, H, O, P and charge; the model's DAG `m1994` has no formula, so this stays inferred.

**Focused tests after the fixes:** `tests/test_r2171_pi_plc.py` + `tests/test_dolichol_route.py`, 10 passed, 82
subtests ([pytest_focused_after_review.log](pytest_focused_after_review.log)). The fixes touch only the scope check,
so the model bytes cannot change; the round-2 builds check that.
