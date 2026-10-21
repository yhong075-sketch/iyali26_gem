# TASK: iYali21 screen under the diagnostic overlay, PO1f + CSM-Leu, glucose 2.414 (2026-10-14)

## Authorization (user, 2026-10-14, this conversation)

The user asked for the paired figure with both models "under iYali26 condition" (diagnostic overlay,
no-solution knockouts counted as lethal). After a read-only applicability check and plan review they
approved: "screen iYali21 on PO1f + CSM-Leu glucose 2.414 under the overlay with the two adaptations
(R2069 from iYali21's form to R01752; skip absent R_MDHy/R_ALCDH_nadp_hi), 4 controls, NGAM
diagnostic, at most 1,213 LPs; then draw both panels as overlay lethal view."

Computation only; the overlay exists in memory only. No model file changes; no commit or push.

## Starting state

- Branch `claude/confusion-matrices-iyali21-7e63ca`, HEAD `f45eb87`, clean.
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.

## Inputs

| Input | Path | SHA-256 |
| --- | --- | --- |
| Model | `model/source/iyli21.xml` | `6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631` |
| Medium | `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` | `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| Strain profile | `state/strain_profiles/po1f_sd_leu.json` | `35307853…` |
| Experimental positives | `state/essentiality/repository/consensus_essential_genes.csv` | `1e887f5a…` |
| Overlay (adapted) | [holistic_overlay_iyali21.py](holistic_overlay_iyali21.py) | `11e28764cadf30f0d916cf18c67f64bd54c085ca18121c6be03eb02fd70aa023` |

## Applicability check (read-only, this session)

Compared with a56f3649 after the medium and PO1f profile:
- R1864, R1866, the six FAO reactions, R1113, R171 GPR, R1384 and R860 are identical → overlay
  steps unchanged.
- R2069 starts from a different equation in iYali21 (extra H⁺ left, H₂O right) → adaptation 1:
  rewritten to the same KEGG R01752 end equation.
- R_MDHy and R_ALCDH_nadp_hi do not exist in iYali21 → adaptation 2: nothing to close.
- Not touched by the overlay but relevant: R171 moves 3 matrix H⁺ per ATP in iYali21 vs 4 in
  a56f3649.
- Only forced bound: `xMAINTENANCE` ≥ 7.8625 in both.

## Method, budget, scoring

[run_screen.py](run_screen.py) is the a56f3649 overlay runner with the model, overlay file and gene
count changed and the iYali26-only assertions removed; [summarize.py](summarize.py) scores and
compares. Budget: WT + 1,058 KOs + 4 controls = 1,063 LPs, NGAM diagnostic ≤ 150, total ≤ 1,213;
1 thread, 900 s wall, no retries; stop on hash change, assertion failure, invalid WT or budget.
Scoring as in the a56f3649 overlay screen (benchmark rule primary; lethal and maintenance views
labelled; mtDNA out of scope; user-defined negatives; development reference).
