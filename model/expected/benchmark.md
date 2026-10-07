# Essentiality benchmark

How single-gene knockout predictions are scored against experiment. The historical contract
this replaces is kept unchanged in
[benchmark_contract_20260905.md](benchmark_contract_20260905.md), with its input identities in
[baseline_manifest_20260905.json](baseline_manifest_20260905.json).

## Reference

| Item | Value |
| --- | --- |
| Experimental positives | `consensus_essential_genes.csv`, 1,612 genes (consensus essential in at least 2 of 3 screens, per the file's `source`/`confidence` columns) |
| Location | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv` (SHA-256 `1e887f5a…`) |
| Medium | SD-Leu, `state/media/sd_leu.csv` (`ed176d26…`) |
| Strain | PO1f runtime profile, `state/strain_profiles/po1f_sd_leu.json` (`35307853…`); applied in memory only |

Copies of the medium and profile used by recorded runs are tracked under
`reports/reference_pipeline_restore_20260909/research/state/`; they are byte-identical to the
research-workspace files above.

## Labels (decided 2026-10-07)

- **Positive (essential in experiment):** a model gene on the consensus list.
- **Negative (non-essential in experiment):** a model gene **not** on the consensus list. This is a
  user-defined negative class; state it whenever FP or TN are reported.
- **Out of scope:** a listed gene with no model counterpart. Report the count and classify the
  reason; it is never counted as FN.

This replaces the earlier positive-only rule ("unlabelled genes are not negatives; no FP/TN").

## Classification

- Essential when the unrounded ratio `KO growth / WT growth` is **strictly below** the cutoff;
  equal to the cutoff is non-essential.
- Cutoffs 1%, 5%, 10% and 15%; **15% is primary** (decided 2026-10-07; used in STATE.md since
  2026-09-15).
- WT and every knockout use the same model file, medium, strain profile and objective.
- Keep raw solver status and growth for every knockout. A failed, missing or non-finite solve is
  not biological evidence; report such cases separately instead of folding them into a class.

## What to report

At each cutoff:
- the counts TP, FN, FP and TN;
- recall `TP/(TP+FN)`;
- specificity `TN/(TN+FP)`;
- precision `TP/(TP+FP)`;
- the number of out-of-scope positives.

Every report also states three things:
- the model SHA-256, medium and profile used;
- that negatives are the user-defined class above;
- that this list has been used during model development, so the result is a development
  reference, not independent validation.
