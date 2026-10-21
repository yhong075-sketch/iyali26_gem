# R2190-RIBULOSE: write D-ribulose as R2190's product (2026-10-19)

## Authorization

The user, in this session on 2026-10-19, after the read-only review
[r2190_identity_review_20261019](../r2190_identity_review_20261019/REPORT.md):

> Approve Option A: opt-in guarded step R2190-RIBULOSE that writes a new D-ribulose[C_cy] as R2190's
> product, built on 285eaee in this worktree, plus the 8 LPs in REPORT §7.3. No commit or push.

The approved design is the review's §7.2 (opt-in flag, wiring choice (i)). The LPs are its §7.3.
Out of scope: a commit or push, any other reaction, the R2041 direction, the R1964 GPR (a separate
session), and the optional single-KO screen.

## Starting state

- Worktree `claude/clever-austin-f6a04e` at `285eaee`. At start the only untracked folder was
  `model/reports/r2190_identity_review_20261019/` (this session's review).
- Baseline model actually used: `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`, SHA-256
  `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`. It is not a designated
  reference model (`model/candidates/README.md`); no reference model is designated.
- Baseline build command: `model/candidates/E5_coq9_resp_pkg_r132_20261012.build.json` (argv, head
  `7dc9511`, dirty status recorded there).
- Research root: `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, passed as
  `--research-root` (the shell's `IYALI26_RESEARCH_ROOT` is unset).
- Baseline target state: [before.json](before.json), written by
  [scripts/record_before.py](scripts/record_before.py).

## Targets

| Entity | Before (a56f3649) | After |
| --- | --- | --- |
| R2190 `D-arabinitol:NAD+ 2-oxidoreductase` | m122 + m1823 ⇌ m123 + **m454** (D-arabinose); [−1000, 1000]; GPR `YALI1F03452g`; EC 1.1.1.250; notes empty | m122 + m1823 ⇌ m123 + **m2061** (D-ribulose); bounds, GPR, name and annotations unchanged; curation notes added |
| `m2061[C_cy]` | absent (no `m2061` in any compartment; highest numeric ID m2060) | new species `D-ribulose_C5H10O5`, C5H10O5, charge 0, `C_cy`; KEGG C00309, CHEBI:17173, MNXM1371690, SBO:0000247 |

Nothing else may change. `m454` keeps its other reactions (R1063, R2041, R254, R255).

## Expected effect

- R2190 has a single-reaction product, so steady state forces it to 0 (provable).
- Closed-boundary ranges of R1964, R2041 and R2042 become [0, 0] (provable from connectivity; review §6).
- WT growth equal to a56f3649's (inferred; all 109 recorded R2190 fluxes are 0).
- Element and charge balance of R2190 unchanged (C5H10O5/0 replaces C5H10O5/0).

## Acceptance criteria

1. The focused test passes, and the full `pytest` matches the baseline run test by test except the
   new tests.
2. A build without the flag reproduces `a56f3649…` byte for byte.
3. The whole-model diff of the new candidate against a56f3649 shows only R2190's product, its notes
   and the one new species.
4. R2190 balances in elements and charge.
5. LPs: WT unchanged within 1e-9 relative; the six closed ranges are 0 within 1e-7.
6. An independent read-only audit finds no undeclared change.

## Compute budget and stop conditions

- **Builds:** two offline `--no-solve` builds (with the flag, and the control without it).
- **pytest:** the full suite twice (baseline and after the change), as the approved plan's step 8.
  Some repository tests run solver calls of their own. They are counted and reported separately
  from the 8 LPs.
- **LPs (8):** context of `model/reports/coq9_respiratory_package_20261008/run_controls.py` (the
  `coq_candidate_checks_20261005/manifest.json` config: PO1f SD-Leu, Gurobi with 1 thread,
  Presolve 0, tolerances 1e-7, 60 s per solve). M0 = a56f3649, M1 = the new candidate.
  - 1–2: WT growth on M0 and M1.
  - 3–8: on M1, all boundary reactions at [0, 0] and forced fluxes relaxed to include 0 (the orphan
    line's `loop_check.py` protocol), then minimize and maximize each of R1964, R2041 and R2042.
    This is done as six explicit LPs rather than cobra FVA, which would add a seventh solve.
- **Stop conditions:**
  - any guard failure or an unexpected diff: stop and report;
  - WT M0 ≠ M1 beyond 1e-9 relative: stop and report R2190's role;
  - a closed range ≠ 0 beyond 1e-7: stop and report the witness flux;
  - a non-optimal status: record the raw value; it is not biological evidence.
- No other LP, screen, commit or push.
