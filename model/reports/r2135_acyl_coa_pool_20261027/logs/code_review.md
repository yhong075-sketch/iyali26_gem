# Code review: R2135 step (2026-10-27)

**Reviewer:** a fresh general-purpose agent following `.claude/agents/code-reviewer.md`. It was given the
task, the files and the tests run, and not the reasoning. It ran read-only in-memory probes and the focused
test.

**Verdict:** no bugs. One medium-low fail-closed finding and three low findings, all fixed before the final
build. The final build is byte-identical to round 1.

| # | Finding (reviewer) | Fix |
| --- | --- | --- |
| 1 | The scope check did not pin the after-equations or check species coverage. An unknown id in R2145's after-state raised `KeyError` after R2135 had already been rewritten (half-applied in memory). A swapped after-state `m68 <=> m567` was accepted. | `AFTER` pinned in code; `spec["species"]` must equal the union of all before / after ids; both checked before any edit. Subtests added: swapped after-state, unknown id, missing species row. |
| 2 | This curation's notes on a reaction still in the before-state were accepted. | New guard "carries this curation's notes without its equation"; subtest added. |
| 3 | Missing tests: a charge change; the m1984 MetaNetX id; conflicting R2145 notes; the reverse half-applied pair; partial notes; a missing R2135. The CLI `getattr` forwarding is not exercised. | All added except the forwarding: read correct, not tested, reported. |
| 4 | `m1984`'s reaction context (only R2134 / R2135 / R2168) was not pinned. | `reactions` added to the m1984 species row and checked; subtest with a new consumer. |
| 5 | Out of scope: `--r1726-cpt-open` does not require `--r2135-acyl-coa-pool`. | Not changed. It belongs to the paused R1726 change, which needs the user's approval. |

**Found while fixing:** `Model.copy()` shares each metabolite's annotation dict with the original. An
identity subtest that mutated it corrupted the class fixture. The test now replaces the dict.

**Focused tests after the fixes:** `test_r2135_acyl_coa_pool.py` + `test_r1726_cpt_open.py`, 9 passed, 52
subtests ([pytest_focused_after_review.log](pytest_focused_after_review.log)).
