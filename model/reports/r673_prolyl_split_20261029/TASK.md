# R673 prolyl-tRNA synthetase rule: cytosolic reaction assigned to YALI1E05929g (2026-10-29)

## 0. Authorization

**User, 2026-10-29**, answering "Approve the R673 plan (R673 ← YALI1E05929g, built on ca43ae0c)?":

> "Approve as written" (reply button whose text was: "Approved: run the R673 plan as written in
> model/reports/r673_prolyl_localization_20261029/TASK.md section 4 on baseline ca43ae0c. No commit or push.")

The plan approved is [../r673_prolyl_localization_20261029/TASK.md](../r673_prolyl_localization_20261029/TASK.md) §4;
the evidence is that stage's [REPORT.md](../r673_prolyl_localization_20261029/REPORT.md) and the OR-rule review §7.2.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-r171-pro-7e2477` @ `ff7fc62` plus the uncommitted R171 round |
| Baseline | `candidates/E5_…_aars_thr_ser_split_r171_atp_synthase_gpr_20261029.xml`, SHA-256 `ca43ae0c12e01c13c7b2ca4792b745d00d4a28ed950f8e401aece37b2a0843d0` (R171 round, uncommitted). Not a designated reference. |
| Its build | `/opt/miniconda3/bin/python`, `PYTHONHASHSEED=0`, the `d8d1bac7` argv + `--r171-atp-synthase-gpr` |
| Before-state | [before.json](before.json): R673 `prolyl-tRNA synthetase`, cytosol, bounds (0, 1000), rule `YALI1D19110g or YALI1E05929g`, notes `PROTEIN_CLASS: 6.1.1.15`; each gene linked to R673 only. R673 is the only prolyl-tRNA synthetase reaction in the model (checked). |
| Code since the R171 pytest after-run | unchanged (`r171_atp_synthase_gpr.py` `315646ab`, `main.py` `70693d3b`, `cli.py` `c023d389`, data `63a7c6ba`; last edit 14:10:56, after-run started about 14:12:40), so that run is this round's "before" ([../r171_atp_synthase_gpr_20261029/logs/pytest_after.txt](../r171_atp_synthase_gpr_20261029/logs/pytest_after.txt): 330 passed, 19 failed, 20 errors) |

## 2. The change

`R673` `YALI1D19110g or YALI1E05929g` → `YALI1E05929g`. `YALI1D19110g` (AIM10-like, predicted mitochondrial) is kept in
the model with no reaction. Notes: curation id, status, authorization, previous rule, unlinked gene, evidence, tier.

Evidence (verified in the prediction stage unless marked):
- `YALI1E05929g` (W29 XP_503566.2, 528 aa): TargetP 2.0 Other 0.99995; DeepLoc 2.1 Cytoplasm 0.752, Mitochondrion
  0.222 (threshold 0.622). Genome annotation "similar to YHR020W" (UniProt: "Putative" ProRS, no location; that it
  is the cytoplasmic enzyme is recalled, uncited), ancestral locus Anc_1.352. *Correction after the audit:*
  `YALI1D19110g`'s annotation also reads "similar to *C. albicans* ProRS cytoplasmic" (REPORT §1).
- `YALI1D19110g` (XP_502869.1, 445 aa): TargetP mTP 0.625 (cleavage 13–14); DeepLoc Mitochondrion 0.844. Annotation
  "similar to YER087W" (AIM10, mitochondrial ProRS), ancestral locus Anc_7.368 (review §7.2).
- Predictions only; no *Y. lipolytica* localisation experiment; start-site difference between the W29 and CLIB122
  gene models of `YALI1E05929g` (prediction REPORT §3).

## 3. Expected effect

`YALI1E05929g` FN → TP at every cutoff (R673 closed → growth 0, review §5); `YALI1D19110g` stays FN (no mitochondrial
prolyl-tRNA reaction, category d). 15 %: 123 / 192 / 42 / 692 → 124 / 191 / 42 / 692. WT unchanged.

## 4. Acceptance criteria

1. Control build with the `ca43ae0c` argv (flag off) reproduces `ca43ae0c…` byte for byte.
2. Candidate seeds 0 and 1 identical.
3. Whole-model diff vs `ca43ae0c`: only R673's rule and notes and `YALI1D19110g`'s reaction link.
4. Sweep: WT unchanged, exactly the §3 transition at every cutoff, unresolved set unchanged, no failing point.
5. Full pytest: only the new tests differ from the before run.

## 5. Plan

As approved: data `model/curation/gpr/r673_prolyl_split.json`; step `platform/scripts/gem_annotate/r673_prolyl_split.py`
(`apply_r673_prolyl_split`); flag `--r673-prolyl-split` after the R171 step; test `platform/tests/test_r673_prolyl_split.py`;
builds; diff; sweep; audit; REPORT; STATE; README row; Escher figure (knockout of `YALI1E05929g`, before vs after).
Budget about 2,164 + 40 + 20 LPs, 30 minutes; stop on any failed criterion. No commit or push.
