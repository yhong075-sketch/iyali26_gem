# TASK: fresh benchmark screen of the newest candidate 2c9cebd0 (2026-10-23)

## Authorization (user, this conversation)

"can you do a screen test for the newest model?" After the plan the user approved: "fresh benchmark
screen of 2c9cebd0 on Task A medium + PO1f, Gurobi at 1e-9 with the knocked-flux check, about 1,083
LPs; compare with the build's sweep and b3b9fe22; figure vs iYali21, report and STATE entry, no
commit." Computation and a figure only; no model change.

## Starting state

- Branch `claude/confusion-matrices-iyali21-7e63ca`, HEAD `ce3f341` (PR #12), clean. `origin/main`
  at `4c4248d` (PR #18).
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.

## Model

`E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml`, SHA-256
`2c9cebd0a5c53ec18bf1a583327099dc9e53fcf6b2fcbedd66cc7b9c54ea5a7b`, the newest candidate on
`origin/main` (commit `47cdd5b`, PR #18). Read from a byte-identical `git show origin/main:` export in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_2c9cebd0_20261022/input/`. Not a designated reference.

## Method, budget, stop

- [run_screen.py](run_screen.py): the b3b9fe22 no-overlay runner with the model changed, Gurobi
  FeasibilityTol/OptimalityTol 1e-9, and a knocked-flux check (closed reactions |v| ≤ 1e-9, else
  not scorable).
- Medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`, Task A medium, benchmark default since
  2026-10-21), PO1f profile `35307853…`, consensus list `1e887f5a…`, mtDNA rule.
- Budget WT + 1,078 KOs = 1,079 LPs; 1 thread, 60 s per LP, 900 s wall, no retries. Stop on hash
  change, invalid WT or budget.

## Scoring and comparisons

Benchmark rule primary; labelled lethal view; four cutoffs; mtDNA out of scope; user-defined
negatives; development reference. No-solve comparisons with the build's sweep
(`reports/combined_tpp_heme_r2193_20261022/ko_compare.tsv` on main) and with the b3b9fe22 screen.
