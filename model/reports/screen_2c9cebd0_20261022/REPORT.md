# Fresh benchmark screen of the newest candidate 2c9cebd0, Task A medium + PO1f (2026-10-23)

Scope and authorization: [TASK.md](TASK.md). Computation and a figure only; no model change; no commit.
Figure (iYali21 vs 2c9cebd0, no overlay, benchmark rule):
[confusion_matrices_15pct_2c9cebd0.png](../figure_iyali21_vs_a56f3649_20261014/confusion_matrices_15pct_2c9cebd0.png).
Scores and comparisons: [summary.json](summary.json).

## Run

- Model `E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml`,
  `2c9cebd0…`, newest candidate on `origin/main` (`47cdd5b`, PR #18); byte-identical export, SHA checked.
  1,078 native genes. Not a designated reference.
- Task A medium (`po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd…`) + PO1f profile; consensus list
  `1e887f5a…`; mtDNA rule.
- Gurobi 13.0.1, FeasibilityTol / OptimalityTol 1e-9, presolve off, 1 thread; knocked-flux check
  (|v| ≤ 1e-9 on every reaction a knockout closes).
- 1,079 LPs (the approved budget), 12 s, no retries, no stop condition: 1,058 optimal, 20 infeasible.
- Numerics over the 1,058 optimal knockouts: largest mass residual 4.1e-10, largest bound violation
  4.2e-10, every closed reaction at exactly 0.

## Results (mtDNA out of scope; verified this session)

| Cutoff | TP | FN | FP | TN | Unresolved | Recall | Precision | Specificity | MCC |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 % | 75 | 238 | 29 | 704 | 16 | 24.0 % | 72.1 % | 96.0 % | 0.306 |
| 5 % | 81 | 232 | 31 | 702 | 16 | 25.9 % | 72.3 % | 95.8 % | 0.321 |
| 10 % | 85 | 228 | 36 | 697 | 16 | 27.2 % | 70.2 % | 95.1 % | 0.318 |
| **15 %** | **93** | **220** | **36** | **697** | **16** | **29.7 %** | **72.1 %** | **95.1 %** | **0.345** |

- WT 0.3124262915 h⁻¹.
- Labelled lethal view at 15 % (the 16 no-solution knockouts counted as lethal): 104 / 220 / 41 / 697,
  MCC 0.356.
- The 20 no-solution knockouts (16 in scope) are the complex III (R305) and complex IV (R304) genes,
  as in the earlier candidates on this medium.
- Positives: 324 of the 1,612 listed genes are model genes; 1,288 are out of scope (not FN).

## Reproduction of the build's own sweep (no new solves)

Against `reports/combined_tpp_heme_r2193_20261022/ko_compare.tsv` on main (`combined_*` columns; export
`4aee38ac…`): no solver status changes, largest growth difference 5.2e-14 over the 1,058 solved
knockouts, **no class change at any cutoff**. The build's 93 / 220 / 36 / 697 / 16 is reproduced by an
independent runner.

## Comparison with b3b9fe22 (no new solves)

**2c9cebd0 is not built on b3b9fe22.** Its build argv has the GPI, TPP/thiamine and heme/iron + R2193
flags but **not** `--lipoylation-module` or `--acyl-acp-chemistry`; those exist only in `9f2356b`
on `claude/combined-candidate-20261019`, which is not on main or any remote branch. The two models
are parallel lines from `164f235f`. The near-identical 15 % counts (93 / 220 / 36 / 697 vs
93 / 220 / 36 / 696) are a coincidence of different genes:

- **Gene sets:** 2c9cebd0 adds the four R2193 iron-import genes (YALI1D08564g, YALI1D08684g,
  YALI1D09277g, YALI1D09308g; all non-essential, the two D09 genes are listed → FN) and lacks the three
  lipoylation genes of b3b9fe22 (YALI1C04991g, YALI1E13472g, YALI1F34729g).
- **23 shared genes change call, the same 23 at every cutoff:**
  - 10 essential in b3b9fe22 → non-essential here (the acyl-ACP / lipoyl genes: YALI1A20089g,
    C07638g, C15991g, C26939g, D18037g, D26360g, E00574g, E22262g, F37498g, F38317g). Their
    b3b9fe22 calls came from the lipoylation and acyl-ACP steps this line does not have.
  - 13 non-essential in b3b9fe22 → essential here (zero growth): 12 heme/iron-route genes on R383,
    R384, R302, R790, R61, R674, R240, R792, R667, R409 (YALI1B03314g, B19490g, C02391g, C07816g,
    D30196g, E35410g, F13039g, F25924g, F27686g, F31135g, F32002g, F34051g) and THI80-like
    YALI1E25443g (R752/R753, TPP line).

  These match the class changes the build's STATE entry reports for the TPP and heme lines. Gene
  names are not verified here; reaction assignments are model GPRs, not experimental validation.

## Comparison with iYali21 (figure)

Same medium, no overlay, benchmark rule: iYali21 58 / 255 / 25 / 694, 10 unresolved, MCC 0.254
(screen of 2026-10-14, Gurobi at 1e-7) vs 2c9cebd0 93 / 220 / 36 / 697, 16 unresolved, MCC 0.345.

## Verification and limits

- Verified this session: all counts, WT, numerics, the build-sweep reproduction and the b3b9fe22
  comparison (raw tables hash-checked).
- The heme open problems of the build (STATE 2026-10-22, items 1–7) still apply; this screen only
  confirms that Gurobi at 1e-9 with the knocked-flux check shows no solver artifact.
- The iYali21 panel was screened at 1e-7 (no near-tolerance fluxes in that model were checked).
- Negatives are a user-defined class; development reference, not independent validation.
- No independent audit.

## Commands (repository root, `IYALI26_RESEARCH_ROOT` set)

```
git show origin/main:model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml > "$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_2c9cebd0_20261022/input/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_20261022.xml"
git show origin/main:model/reports/combined_tpp_heme_r2193_20261022/ko_compare.tsv > "$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_2c9cebd0_20261022/input/ko_compare_build_sweep.tsv"
"$PY" -B model/reports/screen_2c9cebd0_20261022/run_screen.py 2c9cebd0
/opt/miniconda3/bin/python -B model/reports/screen_2c9cebd0_20261022/summarize.py
/opt/miniconda3/bin/python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure_2c9cebd0.py
```

Raw per-knockout table, solver journal and WT fluxes:
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_2c9cebd0_20261022/screen_2c9cebd0/`.
