# Code reviews of the R1726 step (2026-10-27)

Both rounds were done by fresh general-purpose agents following `.claude/agents/code-reviewer.md`. Each was
given the task, the files and the tests run, and not the reasoning. Both were read-only.

## Round 1 (first run, on 00e8e397)

**Verdict:** no bugs, four low findings, all fixed:
1. A stale "runs last" comment.
2. The R1725 guard did not pin stoichiometry and species. It now does, from the curation file's
   `kept_closed`.
3. Untested guards: the scope fields, a missing R1726, an R1725 name change, and guard messages. All are
   now tested by message.
4. No build-context requirement. A flag requirement was added (first `--cofactor-biomass-terms`; in the
   re-run `--r2135-acyl-coa-pool`).

## Round 2 (re-run, on e9d35857)

**Verdict:** no bugs. The reviewer checked all 8 combinations of lipoylation / cofactor / R2135 flags. The
first error always names the deepest missing flag, and every check runs before file I/O. It found three low
items, all fixed:
1. A closed R1726 that still carries this curation's notes was silently reopened. A new guard now raises
   ("closed but carries this curation's notes"), matching the R2135 step.
2. Test gaps: no subtest changed an R1725-only species, and none covered closed-with-notes. Both were added
   (an `m1633` formula change; a notes-only revert).
3. A stale comment ("after the cofactor step"), now "requires and runs after the R2135 step".

**Focused tests after round 2:** 9 passed, 54 subtests ([pytest_focused_after_review.log](pytest_focused_after_review.log)).
The final build is byte-identical to round 1 (`ce9e1112…`).
