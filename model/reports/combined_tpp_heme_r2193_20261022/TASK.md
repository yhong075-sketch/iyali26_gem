# Combined candidate: TPP line + heme/iron + R2193 line (2026-10-22)

## Authorization

User, 2026-10-22, after reading the plan in this session:

> Approved: combine the TPP and heme/iron + R2193 lines as planned (heme guard accepts the reviewed
> TPP-line biomass as a second before-state, combined build with seed and control checks, three-model
> knockout sweep on Gurobi at 1e-9 with the knocked-flux check, audit, report). No commit or push.

## Starting state

- **Branch:** `claude/science-test-gap-fix-5eaeae` at `8e4e879`. Its tree is identical to `origin/main`
  `2913ac0` (PR #17 merged).
- **Models** (none is a reference model):

| Role | File | SHA-256 |
| --- | --- | --- |
| Base | `E5_coq9_resp_pkg_r132_combined_gpi_coproducts_20261021.xml` | `cf0825b86c3533d8a12f697a994eba3d740671a68a367143ee32cbfdc8686b7c` |
| TPP parent | `E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_v2_20261022.xml` (PR #16) | `2fe51cc2b3834bb7303d6702f4f68ac7bec2dd5353cafcd635803fd6e74c174e` |
| Heme/iron + R2193 parent | `E5_coq9_resp_pkg_r132_combined_gpi_coproducts_heme_iron_r2193_gpr_v3_20261022.xml` (PR #17) | `e1586027cd5e3727a75b237097d3578d2eba976b9162294f0d6ae2b8be1a127b` |

- **Before-state:** [before.json](before.json), from [scripts/before_state.py](scripts/before_state.py).
  - The TPP parent's `biomass_C` is the base's plus `m975[C_cy]` (TPP) at −1.2e-6 and 7
    `tpp_requirement_*` notes, equal to `tpp_requirement.json`'s recorded after-state.
  - A read-only comparison in the planning turn showed the TPP line changes only R1342, `biomass_C`
    and 7 thiamine-phosphate species. It touches nothing else that the heme/iron or R2193 guard checks.
  - The heme/iron guard is the only clash: it requires the base's exact `biomass_C` composition and
    note keys.

## Target

1. **Guard.** `heme_iron_biomass.json` and `apply_heme_iron_biomass` accept exactly two reviewed
   starting states for `biomass_C`:
   - the base (as before);
   - the base plus the TPP term and its 7 note keys (the TPP line's reviewed after-state).

   Any other composition still stops the build. No model content changes, and existing builds
   reproduce.
2. **Combined candidate.** The base argv plus
   `--tpp-requirement --thiamine-phosphate-chemistry --thiamine-proton-lock --heme-iron-biomass --r2193-iron-gpr`,
   written to
   `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml`.

## Expected effect (inferred)

- `biomass_C` = base + TPP + heme a + iron, with all three curations' notes.
- The other changes are the union of both lines: R1342 closed, the thiamine-phosphate formulas, the
  proton-lock notes, the R2193 GPR, and four added genes.
- WT growth slightly below both parents: TPP and heme each cost about 4e-9 to 5e-7.
- Knockout classes = the union of both lines' changes against `cf0825b8`: the 13 heme/iron genes, the
  THI80-like YALI1E25443g, and the R2193 gene swap.

## Acceptance criteria

1. **Tests.**
   - Focused tests of both lines pass, including a new case: the heme/iron step applies on top of an
     in-memory TPP-requirement model, and still rejects any other composition.
   - The full suite is identical test by test to the last full run on the pre-merge code
     ([pytest_after_v3.log](../r2193_iron_gpr_change_20261022/logs/pytest_after_v3.log)), plus `main`'s
     TPP/thiamine tests and the new cases.
2. **Builds:**
   - TPP-only argv: reproduces `2fe51cc2`;
   - heme-only argv: reproduces `e1586027`;
   - combined argv: identical bytes under `PYTHONHASHSEED=0` and `1`.
3. **Whole-model diffs:**
   - combined against `2fe51cc2`: only the heme/iron and R2193 targets;
   - combined against `e1586027`: only the TPP-line targets.
4. **WT** on all three models (Task A medium + PO1f): optimal, with combined ≤ each parent.
5. **Knockouts** on all three models (Task A medium + PO1f):
   - Gurobi, FeasibilityTol and OptimalityTol 1e-9, 60 s per LP, model reset before each knockout;
   - targeted check: knocked-out reactions |v| ≤ 1e-9, other bound or |S v| violations ≤ 1e-6, all
     magnitudes recorded;
   - flagged points get a fresh re-solve (at most 20);
   - class changes at 1/5/10/15 % reported pairwise; scores per benchmark.md.
6. Independent read-only audit, REPORT.md, a STATE.md entry and a candidates README row.

## Compute budget and stop condition

- **Budget:**
  - 4 builds;
  - about 10 diagnostic LPs;
  - 3 × (about 1,076 + 1) knockout solves;
  - at most 20 fresh re-solves.
- **Stop and report on:**
  - a guard rejection in the combined build;
  - a control build that differs from its parent;
  - a diff outside the targets;
  - seed-dependent bytes;
  - more than 20 flagged knockouts;
  - a test regression.

## Environment

- `.venv` (Python 3.13.5, cobra 0.30.0, optlang 1.8.3, pytest 9.0.2) and Gurobi (academic license).
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- Medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) and profile `po1f_sd_leu.json`
  (`35307853…`).

## Definitions fixed before the sweep (2026-10-22)

- **Flagged point** (gets a fresh re-solve; at most 20 in total):
  - an optimal point that fails the targeted check; or
  - any status other than optimal or infeasible, for example `time_limit`.
- **Infeasible knockouts** are the solver's verdict, not a stale point. They are reported as unresolved
  without a re-solve. This is the same rule as the R2193 sweep, which its audit noted was not
  separately approved; it is stated here up front.
- **The solves above** used 12 LPs: 3 per model, including the base `cf0825b8` for reference. The
  estimate was about 10.

## After the independent audit (2026-10-22)

See [AUDIT.md](AUDIT.md) and [REPORT.md §6](REPORT.md#6-independent-audit-and-responses).

- **Builder rule added** (`main.py`): `--heme-iron-biomass` with `--tpp-requirement` requires
  `--thiamine-phosphate-chemistry` and `--thiamine-proton-lock`. The audit showed that the TPP
  requirement alone plus heme passed every guard, a combination nobody reviewed.
  - Focused tests rerun ([log](logs/pytest_focused_after_audit.log)).
  - The combined rebuild still gives `2c9cebd0`.
  - The unreviewed mix is refused by the real CLI.
- **Extra work beyond the plan:** 2 builds (the rebuild and the refused mix, which wrote nothing), and a
  second full test-suite run.
