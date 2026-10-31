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
| Medium | Task A medium, [`model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`](../conditions/media/po1f_csm_leu_w29uptake_20261012.csv) (`6eca70fd…`); default since 2026-10-21 (see below) |
| Strain | PO1f runtime profile, `state/strain_profiles/po1f_sd_leu.json` (`35307853…`); applied in memory only |

Copies of the SD-Leu medium (the default before 2026-10-21) and of the profile are tracked under
`reports/reference_pipeline_restore_20260909/research/state/`; they are byte-identical to the
research-workspace files.

## Medium (decided 2026-10-21)

- The default is the Task A medium with the PO1f profile above. It lists the same 35 exchanges as
  SD-Leu. Glucose uptake is 2.414 instead of 10, and adenine, the 11 amino acids and iron are
  scaled by 2.414/10. Uracil stays 0.01607; the profile checks it and raises it to 1000 for
  static FBA.
- The 2.414 is a W29 uptake measured on a different defined medium (Worland et al. 2024; the
  link to the data repository is inferred, pending user confirmation). The finite amino-acid
  bounds are concentration-scaled assumptions, not uptake rates.
- Runs before 2026-10-21 used SD-Leu, `state/media/sd_leu.csv` (`ed176d26…`), unless they say
  otherwise. Scores on the two media are not directly comparable; when comparing, state which
  medium each result used.

## Labels (decided 2026-10-07)

- **Positive (essential in experiment):** a model gene on the consensus list.
- **Negative (non-essential in experiment):** a model gene **not** on the consensus list. This is a
  user-defined negative class; state it whenever FP or TN are reported.
- **Out of scope:** a listed gene with no model counterpart. Report the count and classify the
  reason; it is never counted as FN.

This replaces the earlier positive-only rule ("unlabelled genes are not negatives; no FP/TN").

**Mitochondrially encoded genes (decided 2026-10-08):** the 16 model genes encoded in the
mitochondrial genome ([mtdna_genes_out_of_scope.json](mtdna_genes_out_of_scope.json)) are out of
scope, because the nuclear CRISPR screens cannot target them. They are excluded from TP, FN, FP
and TN and reported separately; they stay in the model and in GPRs. Reports made before this
date counted them as negatives; when comparing with those, also give the scores without this
exclusion.

**Genes no screen guide targets (decided 2026-10-30):** model genes that no guide of the Cas9 or
Cas12a library targets ([untargeted_genes_out_of_scope.json](untargeted_genes_out_of_scope.json),
currently `YALI0F04114g`, COX9) are out of scope in the same way: excluded from TP, FN, FP and TN,
reported separately, kept in the model and in GPRs. Evidence and decision:
[infeasible_ko_rule_20261030](../reports/infeasible_ko_rule_20261030/REPORT.md).

## Classification

- Essential when the unrounded ratio `KO growth / WT growth` is **strictly below** the cutoff;
  equal to the cutoff is non-essential.
- Cutoffs 1%, 5%, 10% and 15%; **15% is primary** (decided 2026-10-07; used in STATE.md since
  2026-09-15).
- WT and every knockout use the same model file, medium, strain profile and objective.
- Keep raw solver status and growth for every knockout. A knockout whose LP is infeasible is scored
  as growth 0 (KO/WT = 0, essential at every cutoff) only when the infeasibility is **certified as
  a maintenance shortfall**:
  1. the sweep solve and a fresh re-solve (new context) both return `infeasible`;
  2. `xMAINTENANCE` is the only reaction with a lower bound above 0 or an upper bound below 0 in
     that context;
  3. in the same knockout context with the `xMAINTENANCE` lower bound set to 0, maximizing
     `xMAINTENANCE` is optimal, passes the targeted check, and its maximum is below the original
     lower bound by more than 1e-6;
  4. WT is optimal in the same context.

  Record the raw status, the maximum maintenance flux and the margin for each such knockout. Any
  other failed, missing or non-finite solve is not biological evidence and stays unresolved,
  reported separately. Reports give the counts under the earlier rule (every non-optimal knockout
  unresolved) alongside (decided 2026-10-30; reference implementation
  [ko_screen_certified.py](../reports/infeasible_ko_rule_20261030/scripts/ko_screen_certified.py)).

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
