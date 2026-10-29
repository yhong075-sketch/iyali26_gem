# Independent audit: R673 prolyl-tRNA synthetase round, prediction and curation stages (2026-10-29)

A fresh read-only general-purpose subagent with no conversation history audited both stages
(`../r673_prolyl_localization_20261029/` and this folder). It edited nothing (`git status --short` identical before and
after), used 8 of 30 LPs, ran the focused test with `PYTHONDONTWRITEBYTECODE=1`, and fetched three UniProt REST entries.
Its scratch scripts (`seq.py`, `mdiff.py`, `spot.py`) were in the session scratchpad and are not kept. Summary of its
report, with the values it returned; fixes at the end.

## 1. Prediction stage — supported (minor caveats)

Source SHAs match `sequence_identity.json` (`w29.faa` `1863842c…`, `UP000182444.fasta` `a0bacab5…`, review TSV
`303f78b1…`, `target.fasta` `c45c02c3…`). `YALI1E05929g` 528 aa (`676045c5…`) = RefSeq XP_503566.2 ("YALI1_E05929g") =
UniProt A0A1H6Q518; `YALI1D19110g` 445 aa (`61876bbf…`) = XP_502869.1 = A0A1D8NEQ1 = CLIB122 Q6C8Z3. The model's gene
annotations name the same accessions. Start site: W29[12:] equals Q6C6Z6[1:]; W29 residue 12 is V, Q6C6Z6 residue 1
is M. All 14 raw-file SHAs match; every score in the prediction REPORT and TSV matches the raw files (TargetP
`output_protein_type.txt`, DeepLoc CSV, thresholds from `results.json`); attention CSVs spell the full sequences;
TargetP per-residue files equal the first 200 residues; zip `testzip` None; PNG SHAs match. Decision rule: written at
21:35:56Z, the first TargetP job ID appears at 21:36:57Z (session transcript), and the outputs meet it. Caveat: the
DeepLoc encoder (Slow / ProtT5) and the TargetP `non-pl` setting are evidenced only by the form reads in the session
and the page text, not by a file in `raw/`.

## 2. SHAs and builds — supported

Baseline `ca43ae0c…`, candidate `6af03c40…` (and their short copies); scratch control = `ca43ae0c…`, seed 1 =
`6af03c40…`. `.build.json` argv without `--output-model`: control = the `ca43ae0c` argv; candidate and seed 1 = that plus
`--r673-prolyl-split`; control `options.r673_prolyl_split=False`. Source and data SHAs in the candidate `.build.json`
match the files; `input_unchanged` and `source_unchanged_during_build` true. The failed first attempt logged "applied"
then `OSError 63`; no temporary or long-name file exists in `model/candidates/`. Seed not in `.build.json` (disclosed).

## 3. Whole-model comparison — supported

Raw XML diff: 11 changed lines (six R673 note lines; the `fbc:or` becomes one `geneProductRef`). cobra comparison of
every reaction, metabolite, gene, compartment, model notes / annotation and the objective: exactly 3 differences (R673
GPR, R673 notes, `YALI1D19110g` reactions → none). ID sets and gene list (1,081) identical. R673 is the only EC 6.1.1.15
reaction; tRNA(Pro) m905 and Pro-tRNA m906 exist only in the cytosol.

## 4. Guarded step and test — supported

The step pins the spec and the full before-state of R673 and both genes, rejects partial and conflicting states, has a
post-condition and restores notes on context exit; scope R673 and two genes; wiring as described. Test: 13 stale-input
+ 5 scope + 4 knockout subtests; run: 6 passed, 22 subtests (40.8 s). Full pytest after-run: 19 failed / 336 passed /
20 errors / 931 subtests vs 19 / 330 / 20 / 909 before; FAILED / ERROR line sets identical (44 lines each).

## 5. Spot check and scores — supported

Own script, basis reset per solve, 8 LPs: WT 0.3116608216336117 on both (R673 0.06432); `YALI1E05929g` KO 0.31166 before,
0.0 after (R673 closed, knocked flux 0); `YALI1D19110g` KO 0.31166 on both; unrelated `YALI1A15879g` 0.038418 on both;
all equal to `ko_compare.tsv`. Matrices recomputed from the TSV: baseline 96/219/35/699, 102/213/37/697, 106/209/42/692,
123/192/42/692; candidate 97/218/35/699, 103/212/37/697, 107/208/42/692, 124/191/42/692; 0 class mismatches; only
`YALI1E05929g` differs; 16 unresolved (11 positives) + 4 mtDNA; coverage 1,612 / 326 / 315 / 1,286; 2,164 LPs, no fresh
re-solve, no failing point.

## 6. Escher numbers — supported

`flux.json`: all four scenarios pass; KO after growth 0.0, R673 (0, 0); KO before = WT for all 10 reactions (R673,
`TRNA_BIOMASS_R673`, R689 0.06432); after: R1225 0.004615, R689 0.005437. The hard-coded reaction list equals the
membership rule stated in the README. R1225 is cytosol → mitochondrion for positive flux. PNG 317 KB, SVG 78 KB.

## 7. Evidence statements and AGENTS.md — partially supported

Correct: predictions called predictions; shared features, dual targeting and start site disclosed; no verified
*Y. lipolytica* name stated; "not a reason: the TP gain" stated. Gaps (UniProt REST, fetched by the auditor):
`YALI1D19110g`'s annotation begins "similar to uniprot|P78600 *Candida albicans* Prolyl-tRNA synthetase cytoplasmic"
(P78600: location Cytoplasm, ECO:0000305), recorded in the review but omitted from both stages, the data and the notes;
"YHR020W = the cytoplasmic ProRS" stated as fact although P38708 is "Putative proline--tRNA ligase YHR020W" with no
location; AIM10 (P39965) is "Probable proline--tRNA ligase, mitochondrial" (Mitochondrion, ECO:0000269, PMID 16823961),
and "probable" was dropped. STATE entry pending (not counted).

## Defects and actions

| # | Defect | Action |
| --- | --- | --- |
| 1 | *C. albicans* "cytoplasmic" annotation hit omitted | Fixed in this REPORT §1, TASK §2 and the prediction REPORT §3; recorded as a limit for the model note (§9) |
| 2 | YHR020W / AIM10 locations stated as fact | Fixed in REPORT §1 and TASK §2 (labelled recalled / uncited, UniProt wording given); model note and step docstring not changed (§9) |
| 3 | pytest "pending" | Fixed: REPORT §7 and acceptance #5 |
| 4 | AUDIT.md missing, §10 pending | Fixed (this file, REPORT §10) |
| 5 | Provenance gaps: encoder and form settings only in the session; seed not in `.build.json`; stale copied docstring | Recorded in REPORT §10; not changed |
