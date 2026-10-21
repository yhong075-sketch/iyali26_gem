# Newest candidate b3b9fe22: screens without and with the diagnostic overlay, PO1f + CSM-Leu, glucose 2.414 (2026-10-21)

Scope and authorization: [TASK.md](TASK.md). Computation and a figure only; no model change; no commit.
Figure (iYali21 vs b3b9fe22, both overlay + lethal view):
[confusion_matrices_15pct_b3b9fe22_overlay_lethal.png](../figure_iyali21_vs_a56f3649_20261014/confusion_matrices_15pct_b3b9fe22_overlay_lethal.png).
Scores and comparisons: [../screen_b3b9fe22_overlay_w29_20261021/summary.json](../screen_b3b9fe22_overlay_w29_20261021/summary.json).

## Run

- Model `E5_coq9_resp_pkg_r132_combined_lipoylation_acylacp_20261020.xml`, `b3b9fe22…` (from `9f2356b`
  on the unmerged branch `claude/combined-candidate-20261019`; byte-identical export, SHA checked).
  1,077 native genes. Not a designated reference.
- Medium `po1f_csm_leu_w29uptake_20261012.csv`, PO1f profile, consensus list `1e887f5a…`; overlay
  `holistic_overlay.py` (`bf6c7ce5…`) applied unchanged (assertions pass; only forced bound NGAM).
- No overlay: 1,078 LPs. Overlay: 1,082 + 111 diagnostic = 1,193 LPs. Total 2,271 of the approved
  ~2,310. All solves optimal or infeasible; no retries; no stop condition.

## Results at 15 % (mtDNA out of scope; verified this session)

| Condition | View | TP | FN | FP | TN | Unresolved | MCC |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| No overlay | Benchmark rule (primary) | 93 | 220 | 36 | 696 | 16 | 0.345 |
| No overlay | Lethal view | 104 | 220 | 41 | 696 | 0 | 0.356 |
| Overlay | Benchmark rule | 93 | 211 | 36 | 694 | 27 | 0.354 |
| **Overlay** | **Lethal view (figure)** | **113** | **211** | **43** | **694** | **0** | **0.378** |

All four cutoffs are in summary.json. Overlay lethal view: recall 34.9 %, precision 72.4 %,
specificity 94.2 %.

- **WT:** 0.3124213 (no overlay; a56f3649 0.3124267); 0.3117468 (overlay; a56f3649 0.3117523).
- **Overlay controls:** complex I KO/WT 0.759; complexes III, IV and ATP synthase no solution. All
  37 no-solution knockouts (34 genes, 3 reactions) "cannot meet maintenance"; the maintenance view
  equals the lethal view.
- The no-solution genes are the same 27 in-scope complex III / IV / ATP synthase genes as in
  a56f3649 (16 without the overlay).

## What changed versus a56f3649 (no new solves; same at both conditions)

**4 genes new to the model**, from the lipoylation module and the R1964 GPR:
YALI1C04991g (R2288, listed, essential → TP), YALI1E13472g (R2287, listed, essential → TP),
YALI1F34729g (R2286, unlisted, essential → FP), YALI1E06616g (R1964, unlisted, non-essential → TN).

**11 shared genes change call, the same 11 at 1, 5, 10 and 15 %:** each goes from non-essential to
an optimal solve at **exactly zero growth** (KO/WT 0.0). 9 are listed (→ TP), 2 unlisted (→ FP).

| Gene | Listed | a56f3649 KO/WT | Model reactions (b3b9fe22) |
| --- | --- | ---: | --- |
| YALI1A20089g | yes | 1.000 | R78, R1400–R1402 |
| YALI1C07638g | yes | 1.000 | R1434, R1435, R1869 |
| YALI1C15991g | yes | 1.000 | R87, R88, R1393 |
| YALI1C26939g | no | 1.000 | R78, R1403–R1405 |
| YALI1D18037g | yes | 0.758 | R2062, R1394–R1396, R2218–R2233, … |
| YALI1D26360g | yes | 0.974 | R357, R360, R361, R364, R365, R613, R1968–R1976, R2281, … |
| YALI1E00574g | yes | 0.983 | R2288, R357, R360, R361, R365 |
| YALI1E07744g | yes | 1.000 | R765, R766 |
| YALI1E22262g | yes | 1.000 | R78, R539 |
| YALI1F37498g | no | 1.000 | R78, R1397–R1399 |
| YALI1F38317g | yes | 1.000 | R78, R1394–R1396, R2218–R2233 |

No gene changes in the other direction. Arithmetic check (no overlay, benchmark rule): a56f3649
82 / 229 / 33 / 697 → +9 TP −9 FN, +2 FP −2 TN from the 11, + 2 TP + 1 FP + 1 TN from the new genes →
93 / 220 / 36 / 696.

**Likely causes (inferred from which reactions the genes carry; mechanisms not traced):**
- R357/R360/R361/R365 and R2288 are the lipoyl-dependent dehydrogenase and lipoylation reactions;
  the lipoylation module adds a lipoyl growth-dilution requirement.
- R78, R1393–R1405 and R2218–R2233 are acyl-ACP / ACP reactions; the acyl-ACP chemistry step changed
  their formulas and H⁺ and made 7 hydrolases forward-only (STATE 2026-10-21).
- R765/R766 (transketolase) lost the fragment YALI1D02625g from its OR rule in the combined
  candidate, leaving YALI1E07744g alone.
Gene names and functions are not verified here; reaction assignments are model GPRs, not
experimental validation. These calls are model consequences of changes already accepted on that
branch, not evidence that the changes are right.

## Comparison to the figure's left panel

iYali21 under the adapted overlay, lethal view: 74 / 246 / 30 / 692 (MCC 0.292)
([../screen_iyali21_overlay_w29_20261014/REPORT.md](../screen_iyali21_overlay_w29_20261014/REPORT.md)).
b3b9fe22: 113 / 211 / 43 / 694 (MCC 0.378). Same caveats as before: the overlay is diagnostic, its
caps were set for iYali26, and iYali21 needed two adaptations.

## Verification and limits

- Verified this session: all counts, WT, controls, diagnostic labels and the gene comparisons (raw
  tables hash-checked against their manifests).
- b3b9fe22 is on an unmerged, unpushed branch; if that branch changes, these results refer to
  `9f2356b`.
- Negatives are a user-defined class; development reference, not independent validation.
- No independent audit.

## Commands (repository root, `IYALI26_RESEARCH_ROOT` set)

```
git show 9f2356b:model/candidates/E5_coq9_resp_pkg_r132_combined_lipoylation_acylacp_20261020.xml > "$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_b3b9fe22_w29_20261021/input/E5_coq9_resp_pkg_r132_combined_lipoylation_acylacp_20261020.xml"
"$PY" -B model/reports/screen_b3b9fe22_w29_20261021/run_screen.py b3b9fe22
"$PY" -B model/reports/screen_b3b9fe22_overlay_w29_20261021/run_screen.py run
/opt/miniconda3/bin/python -B model/reports/screen_b3b9fe22_overlay_w29_20261021/summarize.py
/opt/miniconda3/bin/python -B model/reports/figure_iyali21_vs_a56f3649_20261014/make_figure_b3b9fe22_overlay_lethal.py
```
