# R673 prolyl-tRNA synthetase pair: presequence prediction before any rule change (2026-10-29)

## 0. Authorization

**Brief** (item 2): the localisation evidence for the prolyl pair `R673` (`YALI1D19110g or YALI1E05929g`) is the
weakest of the three synthetase pairs (synteny plus an N-terminal heuristic). Run a presequence predictor (TargetP 2.0 /
DeepLoc 2.1, as in `reports/glypro_localization_predictors_20260924/`) only after the user approves the web
submission; if it agrees that `YALI1E05929g` is cytosolic, propose R673 ← `YALI1E05929g` with `YALI1D19110g` unlinked.

**User, 2026-10-29**, answering "Submit the two prolyl-tRNA synthetase sequences to TargetP 2.0 and DeepLoc 2.1?":

> "Yes, run predictors" (reply button whose text was: "Approved: for round 2, submit the full sequences of
> YALI1E05929g and YALI1D19110g to TargetP 2.0 (Non-Plant) and DeepLoc 2.1 through the web forms, as in
> reports/glypro_localization_predictors_20260924. Draft TASK.md first and keep the raw outputs. No model change until
> I approve the plan.")

Scope of this stage: the two submissions and their raw outputs. No model, GPR, data, code or test change; no commit.

## 1. Starting state

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-r171-pro-7e2477` @ `ff7fc62` with the uncommitted R171 round |
| Model consulted (read only) | `ca43ae0c…` (R171 round); `R673` = `YALI1D19110g or YALI1E05929g`, identical to the review's `965ab947` rule |
| Prior evidence | OR-rule review §7.2: `YALI1E05929g` "similar to YHR020W" (cytoplasmic ProRS), `YALI1D19110g` "similar to YER087W" (AIM10, mitochondrial ProRS), from the annotation strings; N-terminal heuristic. No *Y. lipolytica* experimental localisation found. |

## 2. Sequences (verified this session; [sequence_identity.json](sequence_identity.json), [target.fasta](target.fasta))

The model genes are W29 (`YALI1`) loci. Each sequence is taken from the W29 RefSeq proteome (GCF_001761485.1, cached
2026-10-19, SHA-256 `1863842c…`) and is identical to the W29 UniProt proteome UP000182444 entry named in the model's
gene annotation (`a0bacab5…`). Script: [scripts/pin_sequences.py](scripts/pin_sequences.py).

| Gene | W29 protein | Length | vs the CLIB122 entry used in the review |
| --- | --- | ---: | --- |
| `YALI1E05929g` | XP_503566.2 = A0A1H6Q518 | 528 | **differs at the N-terminus**: CLIB122 Q6C6Z6 (517 aa) starts `MPPKAAK…`; W29 starts `MTRKTNDATNKVPPKAAK…`. Residues 13–528 of W29 equal 2–517 of Q6C6Z6, so the two gene models differ only in the start site (11 extra residues, and V12 vs M1). The review's N-terminal heuristic used the 517-aa model. |
| `YALI1D19110g` | XP_502869.1 = A0A1D8NEQ1 | 445 | identical to Q6C8Z3 |

The submitted sequences are the W29 ones (the model's genes). The start-site difference matters for a presequence
predictor; which start is used in vivo is unknown.

## 3. Plan of this stage

1. Submit [target.fasta](target.fasta) (both sequences in one job) to
   [TargetP 2.0](https://services.healthtech.dtu.dk/services/TargetP-2.0/) (organism group Non-Plant, long output) and
   [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) (High-quality / ProtT5, long output), through
   the web forms in the built-in browser. No e-mail address, no account, one submission each.
2. Save the raw outputs (result pages and downloadable files) to [raw/](raw/), with job IDs, dates and SHA-256.
3. Record the scores in a small table and stop. Then present the result and, if it supports the synteny assignment,
   a separate curation-change plan for R673 for approval.

**Decision rule, stated before the run:** the R673 edit is proposed only if `YALI1E05929g` gets no mitochondrial
targeting call (TargetP not mTP; DeepLoc not Mitochondrion above its threshold) **and** `YALI1D19110g` gets a
mitochondrial call from at least one predictor. If the predictions disagree with the synteny assignment or are
inconclusive, stop and report. Predictions are labelled as predictions; they do not establish localisation.

**Not submitted (would need approval):** the 517-aa CLIB122 start-site variant of `YALI1E05929g` as a sensitivity
check; SignalP / DeepTMHMM.

**Outcome (2026-10-29, [REPORT.md](REPORT.md)):** the rule is met. TargetP: `YALI1E05929g` Other 0.99995,
`YALI1D19110g` mTP 0.625. DeepLoc: `YALI1E05929g` Cytoplasm 0.752 (Mitochondrion 0.222 < 0.622),
`YALI1D19110g` Mitochondrion 0.844.

## 4. Proposed curation-change plan for R673 (for approval; nothing written yet)

**Change:** `R673` (prolyl-tRNA synthetase, cytosol) `YALI1D19110g or YALI1E05929g` → `YALI1E05929g`.
`YALI1D19110g` (AIM10-like, predicted mitochondrial) is kept in the model with no reaction. The model has no
mitochondrial prolyl-tRNA reaction (checked: R673 is the only prolyl-tRNA synthetase reaction in `ca43ae0c`), so
`YALI1D19110g` stays FN (category d). Notes record the curation id, authorization, previous rule, the unlinked gene
and the evidence (synteny annotation, TargetP / DeepLoc predictions, start-site caveat).

**Baseline:** `ca43ae0c…` (R171 round, uncommitted), to keep one line; `d8d1bac7` is the alternative.

**Expected effect** (review variant `aaRS_split`, verified in memory on `965ab947`; R673 closed → growth 0):
`YALI1E05929g` FN → TP at every cutoff, nothing else. 15 %: 123 / 192 / 42 / 692 → 124 / 191 / 42 / 692.

**Steps** (same pattern as `aars_thr_ser_split.py` and the R171 round):
1. `model/reports/r673_prolyl_split_20261029/` with TASK.md and before.json.
2. Data `model/curation/gpr/r673_prolyl_split.json`; step `platform/scripts/gem_annotate/r673_prolyl_split.py`
   (`apply_r673_prolyl_split`, pins R673 and both genes, fails closed); flag `--r673-prolyl-split` after the R171 step;
   test `platform/tests/test_r673_prolyl_split.py`; full pytest before and after.
3. Builds with `/opt/miniconda3/bin/python`: control (the `ca43ae0c` argv) byte-identical; candidate + flag; seed-1
   check.
4. Whole-model diff (only R673 rule and notes, `YALI1D19110g` link); checked knockout sweep (2,164 LPs); four-cutoff
   scores.
5. Independent audit, REPORT, STATE, candidates README row, Escher before/after figure of R673 (knockout of
   `YALI1E05929g`). No commit or push.

**Budget:** about 2,164 LPs for the sweep, 40 for the audit, 20 for the figure; 30 minutes. Same stop conditions as
the R171 round.
