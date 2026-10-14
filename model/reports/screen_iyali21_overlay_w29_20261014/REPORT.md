# iYali21 under the adapted diagnostic overlay, PO1f + CSM-Leu, glucose 2.414 (2026-10-14)

Scope and authorization: [TASK.md](TASK.md). Computation only; overlay in memory only; no commit.
Figure (both models, overlay, lethal view):
[confusion_matrices_15pct_overlay_lethal.png](../figure_iyali21_vs_a56f3649_20261014/confusion_matrices_15pct_overlay_lethal.png).

## Run

- `model/source/iyli21.xml` (`6974b758…`), `po1f_csm_leu_w29uptake_20261012.csv`, PO1f profile;
  overlay [holistic_overlay_iyali21.py](holistic_overlay_iyali21.py) (`11e28764…`) with the two approved
  adaptations (R2069 from iYali21's form to KEGG R01752; R_MDHy / R_ALCDH_nadp_hi absent). All
  overlay assertions passed; only forced bound afterwards is `xMAINTENANCE` ≥ 7.8625.
- 1,144 of 1,213 LPs (1,063 screen + 81 diagnostic), no retries, no stop condition.

## Results at 15 % (mtDNA out of scope; verified this session)

| View | TP | FN | FP | TN | Unresolved | Recall | Precision | Specificity | MCC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Benchmark rule (primary) | 58 | 246 | 25 | 692 | 21 | 19.1 % | 69.9 % | 96.5 % | 0.261 |
| Lethal view (labelled) | 74 | 246 | 30 | 692 | 0 | 23.1 % | 71.2 % | 95.8 % | 0.292 |

Other cutoffs, lethal view: 56/264/23/699 (1 %), 63/257/25/697 (5 %), 67/253/30/692 (10 %). All views
are in [summary.json](summary.json). The maintenance view equals the lethal view.

- **Controls:** WT 0.2982130 (identical to no overlay); complex I KO/WT 1.000; complex IV 0.936;
  complex III and ATP synthase no solution.
- **No-solution knockouts:** 25 genes (21 in scope: 10 complex III R305, 11 ATP synthase R171; 16
  listed) plus 2 reaction controls. All 27 "cannot meet maintenance": maximum maintenance ATP 1.938
  (R305) and 5.179 (R171), below 7.8625.
- **Versus no overlay** (same model and medium; no solves): the only call changes at every cutoff are
  the 14 R171 genes (non-essential at KO/WT 1.0 → no solution), as in a56f3649. Three ratios move
  without crossing a cutoff: YALI1C08262g (R344) and YALI1D15424g (R642) 0.84299 → 0.84295,
  YALI1F01981g (R940) 0.95871 → 0.95639.

## Comparison with a56f3649 under the same overlay (lethal view, 15 %)

iYali21 74 / 246 / 30 / 692 (MCC 0.292) against a56f3649 102 / 220 / 40 / 695 (MCC 0.354). The
overlay adds the same 14 R171 genes in both models, so the gap (0.266 → 0.331 without the overlay,
lethal view) carries over: 0.292 → 0.354.

Differences that remain under the "same" overlay: iYali21 survives a complex IV knockout (0.936)
and complex I knockout (1.000); its R171 moves 3 H⁺ per ATP (4 in a56f3649); and the overlay steps
for R_MDHy / R_ALCDH_nadp_hi have no iYali21 counterpart. The caps were set for iYali26.

## Verification and limits

- Verified this session: all counts, controls, diagnostic labels and the no-overlay comparison (raw
  table hash checked).
- The overlay is diagnostic; the caps are assumptions; Task A audit F3 (R1864/R1866 +1 H/+1 charge;
  m1949 no formula) applies. The two adaptations are this task's, approved by the user.
- Negatives are a user-defined class; development reference, not independent validation.
- No independent audit.

## Commands (repository root, `IYALI26_RESEARCH_ROOT` set)

```
"$PY" -B model/reports/screen_iyali21_overlay_w29_20261014/run_screen.py --self-check
"$PY" -B model/reports/screen_iyali21_overlay_w29_20261014/run_screen.py run
/opt/miniconda3/bin/python -B model/reports/screen_iyali21_overlay_w29_20261014/summarize.py
/opt/miniconda3/bin/python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure_overlay_lethal.py
```

Raw tables: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_iyali21_overlay_w29_20261014/`.
