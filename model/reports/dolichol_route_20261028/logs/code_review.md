# Code review of the dolichol-route step (2026-10-28)

One round, by a fresh general-purpose agent following `.claude/agents/code-reviewer.md`. It was given the
task, the changed files and the tests run, and not the reasoning. It was read-only (it ran an in-memory
probe and the flag enumeration, and edited nothing).

**Verdict:** no bug that can produce a wrong model file. Three items:

1. **Low, fixed.** `_scope_ok` did not check that every species of the new reaction's stoichiometry, plus
   `m503[C_er]`, has a row in `species`. A metabolite missing from the model would raise `KeyError` only
   after R2178 was closed and R_DOLSYN added. In the builder that aborts before any write, but an in-memory
   caller could keep a half-edited model. A metabolite in the model but without a row would enter R_DOLSYN
   unchecked. The reviewer reproduced both in memory. **Fix:** the scope check now requires
   `set(new.stoichiometry) | {m503[C_er]} <= set(species)`; a scope test adds `m1[C_cy]` to the stoichiometry.
2. **Rule, disclosed (REPORT.md).** R_DOLSYN adds 2 H+ while the charge does not close (−2), and its mass
   balance cannot be checked (m2019 has no formula, m499 has a C15 formula). platform/AGENTS.md says no
   H+ is added unless H, O and charge all close. Here the charge cannot close in the model's own NAD / NADP
   convention (NADP and NADPH both charge 0); the 1 H+ form closes neither H nor charge. See TASK.md
   section 10.
3. **Minor, fixed.** Untested guards: `{mid} missing`, `biomass_C missing`, and the reaction-context guard
   for m2019 and m503. Four subtests were added. The `Group ... missing` message is unreachable through
   the model (R287's group guard fires first); the duplicate-name case tests the uniqueness half.

**Checked and fine (reviewer):** every model-side check runs before the first edit; the three edits are
all present, all absent, or rejected; idempotent after a deterministic SBML write and reload; all 15
combinations of missing predecessor flags raise before any build, naming the missing predecessor of the
most downstream flag, as the R1726 step does; no `<`, `>` or `&` in notes; no PYTHONHASHSEED dependence
(no genes added, metabolites added sorted); the step runs after every other biomass_C edit.

The builds were made after these fixes (first and only round).
