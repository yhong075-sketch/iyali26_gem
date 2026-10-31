# Certified-infeasible knockouts scored as growth 0; COX9 out of scope (2026-10-30)

Plan and authorization: [TASK.md](TASK.md) (approved 2026-10-30, D1–D3 yes). Script review: [REVIEW.md](REVIEW.md).
No model change, no commit or push.

## 1. Question

On every recent candidate, 16 scored knockouts (all complex III R305 or complex IV R304 members) have an infeasible
LP. Under the earlier rule they were "unresolved" and left out of every score, although the model gives a definite
answer for them: the cell cannot cover its maintenance ATP.

## 2. Code logic

[scripts/ko_screen_certified.py](scripts/ko_screen_certified.py) is a one-model copy of the scorer used since the
dolichol round (`ko_compare_two.py`, `bb3118a6…`) with the same `checked_solve.py` (`69c42690…`): Gurobi
FeasibilityTol/OptimalityTol 1e-9, a basis reset per solve, the targeted check, at most 20 fresh re-solves and a 60 s
timeout. Added: each knockout with sweep status `infeasible` gets a certificate in a new context. The certificate
checks three things:

1. A re-solve of the knockout is infeasible again.
2. The only forced flux is `xMAINTENANCE`.
3. With the `xMAINTENANCE` lower bound set to 0, the maximum of `xMAINTENANCE` is optimal and passes the targeted
   check, and that maximum stays below the original bound by more than 1e-6.

WT must be optimal (checked before the sweep). Rule A scores every non-optimal knockout as unresolved (earlier
rule). Rule B scores a certified knockout as growth 0. Both rules are scored with the mtDNA scope and with the
mtDNA + untargeted scope. [scripts/run_screen.zsh](scripts/run_screen.zsh) gates the model copy and every input on
its SHA-256 (TASK §1) and requires `IYALI26_RESEARCH_ROOT`.

## 3. Results (verified: observed in this run)

Model `e1e7c868…` (not a designated reference); Task A medium `6eca70fd…`; PO1f `35307853…`; list `1e887f5a…`;
1,123 LPs (1,082 sweep, 1 WT ceiling, 0 fresh, 40 certificate), within the 1,200 budget.
Outputs: [ko_screen.tsv](ko_screen.tsv), [ko_screen_summary.json](ko_screen_summary.json), [logs/](logs/).

**Acceptance:**

| Check | Result |
| --- | --- |
| A1: WT growth | 0.3116608216, identical to the aconitase run (difference 0.0) |
| A2: rule A reproduces the earlier scores | identical at all four cutoffs; the same 20 non-optimal knockouts |
| A3: certificates | all 20 infeasible knockouts certified (16 scored + 4 mtDNA). Each: re-solve `infeasible`; forced flux only `xMAINTENANCE`; lower bound 7.8625; maximum maintenance 4.828; margin 3.0345; ceiling check passed (knocked flux 0, bound violation 0, residual ≤ 2.9e-15) |
| A4: only certified knockouts change class | the 20 genes whose class changes are exactly the 20 certified |
| WT ceiling | maximum maintenance 63.333 (optimal, check passed) |

**Scores** (negatives are the user-defined class; mtDNA genes out of scope; development reference, not
independent validation):

| Cutoff | Rule A, mtDNA scope (earlier) | Rule B, mtDNA scope | **Rule B, mtDNA + untargeted (new primary)** |
| --- | --- | --- | --- |
| 1 % | 97 / 218 / 35 / 699, 16 unres. | 108 / 218 / 40 / 699 | **108 / 218 / 39 / 699** |
| 5 % | 104 / 211 / 37 / 697, 16 unres. | 115 / 211 / 42 / 697 | **115 / 211 / 41 / 697** |
| 10 % | 108 / 207 / 42 / 692, 16 unres. | 119 / 207 / 47 / 692 | **119 / 207 / 46 / 692** |
| 15 % | 125 / 190 / 42 / 692, 16 unres. | 136 / 190 / 47 / 692 | **136 / 190 / 46 / 692** |

At 15 %, recall goes 0.3968 → 0.4172, precision 0.7485 → 0.7473 and MCC 0.4255 → 0.4344 (new primary). These match
the prediction written in TASK.md §2 before the run. Under rule A with the new scope, the counts are unchanged and
15 knockouts are unresolved. No knockout remains unresolved under rule B.

## 4. What it means

- The infeasibility is not numerical. All 20 knockouts fall short by the same certified margin: they can make at
  most 4.828 maintenance ATP against a required 7.8625. That is a model prediction ("cannot sustain itself"), now
  scored as growth 0.
- The 16 become 11 TP and 5 FP, or 4 FP once COX9 is out of scope. The remaining FP are CYT1 (A17570g), QCR10
  (C17100g), QCR6 (F32010g) and COX4 (E23635g).
  - CYT1 and COX4 are core subunits with borderline screen calls. Whether the screens missed them is inferred.
  - QCR6 and QCR10 are accessory subunits. Their yeast deletions keep respiration at 30 °C (abstracts read
    2026-10-30).
- The maximum 4.828 = 2 × the glucose cap 2.414 (inferred: fermentative ATP only). The result therefore depends on
  NGAM 7.8625 (inherited from `iyli21.xml`; origin not traced) and on the glucose uptake, which was measured on a
  different medium.

## 5. Could the method be better?

- The certificate covers maintenance shortfalls only. A future infeasibility from another bound stays unresolved.
  A general alternative is a Farkas certificate (Gurobi `InfUnbdInfo`) recorded per knockout. Test: compare it with
  this certificate on the same 20 knockouts.
- The ceiling maximum relies on Gurobi's optimality status. `check()` tests the primal point only. A dual-bound
  check, or an exact re-solve at a looser tolerance, would make condition 3 independent of one solver call.

## 6. Verification level

- **Provable:** a certified knockout has no flux state with non-negative growth, because maintenance alone is
  unreachable.
- **Statically verified:** rule B can change only certified infeasible knockouts (independent review).
- **Observed in this run:** A1–A4 and the numbers above.
- **Untested:** the rule on other media or models; whether complex III/IV loss is lethal in *Y. lipolytica*
  (screen data only); the complex III per-subunit GPR review (not done; complex IV was reviewed on 2026-10-08).

## 7. Files changed outside this folder

- [`expected/benchmark.md`](../../expected/benchmark.md): the Classification bullet replaced by the approved rule
  text (D1), and the untargeted-gene scope paragraph added under Labels (D2).
- [`expected/untargeted_genes_out_of_scope.json`](../../expected/untargeted_genes_out_of_scope.json): byte-identical
  copy of [scripts/untargeted_genes_out_of_scope.json](scripts/untargeted_genes_out_of_scope.json) (`1719dd53…`).
- D3: later rounds copy `scripts/ko_screen_certified.py` instead of `ko_compare_two.py`. A two-model comparison
  needs this script run once per model, or a two-model adaptation that adds the certificate.
