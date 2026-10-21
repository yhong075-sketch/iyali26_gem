# GPI chain links: step 1, three ER links and an open exit (2026-10-21)

## Authorization

User, 2026-10-21, after reading the plan in this session:

> Approved: implement the GPI chain plan with the defaults as written (baseline 164f235f, open
> demand, co-product fix as a later change, KO comparison included). No commit or push.

The plan answered the user's request to "fix the gap as the figure showed". The figure was a
Claude Science in-memory test: the GPI chain R371…R376 is blocked, and it carries flux once five
links are added.

## Starting state

- **Branch:** `claude/gpi-chain-gaps-398472` at `a78a49d` (main after PR #9). Clean worktree.
- **Baseline model:** `model/candidates/E5_coq9_resp_pkg_r132_combined_20261020.xml`,
  SHA-256 `164f235f4fe823181df0cc83082d6058ee9cf59675777adfe93557e9a77583a6`. It has 2,306
  reactions, 1,876 metabolites and 1,074 genes. It is a candidate, not a designated reference.
- **Before-state of every target:** [before.json](before.json), written by
  [scripts/before_state.py](scripts/before_state.py).
- **Claude Science's numbers:** growth 0.31240 and R1839 = 0.00154 do not reproduce on this file.
  Here the growth is 1.355248 on the model's own bounds, 1.415816 on SD-Leu and 1.430806 on
  SD-Leu + PO1f, and R1839 reaches about 0.007. The model or medium it ran is unknown. Its
  structural finding does reproduce (see below).

## What was reproduced before any change (2026-10-20, in memory, GLPK)

- With fraction of optimum 0, the FVA maximum of R371, R113 and R372–R376 is 0 in all three
  conditions.
- With the five links added and the exit forced to ≥ 1e-5, the solve is optimal and every chain
  reaction carries 1e-5. Removing step 1, the PE link, the DAG link or the acetate link makes it
  infeasible.
- With the exit open (0–1000), growth is unchanged and the chain carries 0 at the optimum.

## Targets

Five new reactions. No existing reaction, metabolite or gene is edited.

| ID | Equation | Bounds | Group |
| --- | --- | --- | --- |
| `R_GPIGNT` | m403[C_cy] + m1648[C_em] + 2 m10[C_cy] → m11[C_cy] + m596[C_er] | 0, 1000 | GPI-anchor biosynthesis (group37) |
| `R_PEt_em_er` | m1633[C_em] ⇌ m599[C_er] | −1000, 1000 | Transport, ER to ER membrane (group87) |
| `R_DAGt_er_em` | m604[C_er] ⇌ m1649[C_em] | −1000, 1000 | Transport, ER to ER membrane (group87) |
| `R_ACt_er_cy` | m597[C_er] ⇌ m202[C_cy] | −1000, 1000 | Transport, ER (group79) |
| `R_DM_GPI` | m603[C_er] → | 0, 1000 | GPI-anchor biosynthesis (group37) |

- All five are gene-less.
- The data file is `model/curation/gapfill/gpi_chain_links.json`. It is applied by
  `apply_gpi_chain_links()` in `platform/scripts/gem_annotate/gpi_chain_links.py` under the
  opt-in flag `--gpi-chain-links`.

## Expected effect

- The seven GPI reactions and the five new ones can carry flux (FVA maximum > 0 at fraction 0).
- WT growth is unchanged, because the exit is open, not forced.
- Single-gene knockouts are expected to be unchanged. New outputs (DAG, Dol-P, CoA, UDP, acetate)
  could in principle rescue a knockout. That is untested, so all genes are compared.

## Acceptance criteria

1. pytest per test equals the baseline run plus the new test file passing.
2. The build with the combined argv plus `--gpi-chain-links` gives the same bytes under
   `PYTHONHASHSEED=0` and `1`.
3. Against `164f235f`, the whole-model diff holds only the five new reactions, their membership in
   groups 37, 87 and 79, and the touched metabolites' reaction links. Anything else stops the work.
4. All five reactions are mass- and charge-balanced, except the boundary reaction `R_DM_GPI`.
   Mass balance is under the model's stored formulas and charges.
5. On SD-Leu + PO1f with GLPK:
   - WT growth equals the baseline (1.430806);
   - the chain FVA maximum is > 0;
   - forcing `R_DM_GPI` to ≥ 1e-5 is optimal;
   - closing each of the four links in turn makes that infeasible.
6. Single-gene knockouts on both models give the same class at 1%, 5%, 10% and 15% KO/WT. Any
   gene whose KO growth changes by more than 1e-6 is listed with its cause.

## Compute budget and stop condition

- **Budget:** about 30 diagnostic LPs plus FVA on 12 reactions. The KO comparison is 2 × (1,074
  genes + WT) single solves. Everything uses the same medium, strain profile, objective and solver
  on both models.
- **Stop and report on:**
  - a guard rejection;
  - a diff outside the targets;
  - seed-dependent bytes;
  - a non-optimal WT;
  - a test regression;
  - any knockout class change (it is reported, never "fixed").

## Environment

- **Interpreter:** the main checkout's `.venv` (Python 3.13.5, cobra 0.30.0, optlang 1.8.3,
  pytest 9.0.2). `pydantic_ai` is not installed there, so three gem_review lit test files fail to
  collect. pytest runs with `--continue-on-collection-errors`, both before and after.
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- **Condition:** `state/media/sd_leu.csv` (`ed176d26…`) and
  `state/strain_profiles/po1f_sd_leu.json` (`35307853…`).

## Not in scope

- The R373 DAG and R375 Dol-P co-products: a later change.
- The R371 extra H⁺.
- The `P*2` formula notation of PI and GlcNAc-PI.
- R113's formula-less product.
- The GPI chain's directions.
- A step-1 GPR.
- The GPI7 and fourth-mannose steps.
- A GPI term in biomass.
