# R673 prolyl-tRNA synthetase rule: cytosolic reaction assigned to YALI1E05929g, candidate 6af03c40 (2026-10-29)

Scope and authorization: [TASK.md](TASK.md) (user, 2026-10-29: "Approve as written" for the plan in
[../r673_prolyl_localization_20261029/TASK.md](../r673_prolyl_localization_20261029/TASK.md) §4 on `ca43ae0c`).
Predictions: [../r673_prolyl_localization_20261029/REPORT.md](../r673_prolyl_localization_20261029/REPORT.md). Audit:
[AUDIT.md](AUDIT.md). Figure: [escher/](escher/README.md). Nothing committed or pushed.

## 1. What changed and why

R673 (prolyl-tRNA synthetase, cytosol) read `YALI1D19110g or YALI1E05929g`, inherited from iYli21: a cytosolic-type
enzyme OR its mitochondrial-type paralog, so neither single knockout closed the reaction. The cytosolic reaction now has
the cytosolic enzyme only: `R673` ← `YALI1E05929g`. `YALI1D19110g` is kept in the model with no reaction.

| Gene | Identity | Localisation evidence | Level |
| --- | --- | --- | --- |
| `YALI1E05929g` (YALI0E05027g; W29 XP_503566.2, 528 aa) | proline–tRNA ligase; annotation "similar to uniprot\|P38708 *S. cerevisiae* Putative prolyl-tRNA synthetase YHR020W … ancestral locus Anc_1.352" | TargetP 2.0 Other 0.99995; DeepLoc 2.1 Cytoplasm 0.752, Mitochondrion 0.222 (threshold 0.622) | annotation string verified this session (`or_rule_review_genes.csv`); predictions verified |
| `YALI1D19110g` (YALI0D15620g; XP_502869.1, 445 aa) | proline–tRNA ligase; annotation "similar to uniprot\|P78600 *Candida albicans* Prolyl-tRNA synthetase cytoplasmic, similar to *S. cerevisiae* YER087W … Anc_7.368" | TargetP mTP 0.625 (cleavage 13–14); DeepLoc Mitochondrion 0.844 | as above |

**The annotation signal is mixed and is not the ground for the change.** `YALI1E05929g`'s best *S. cerevisiae* hit is
YHR020W, which UniProt ([P38708](https://rest.uniprot.org/uniprotkb/P38708)) names "Putative proline--tRNA ligase
YHR020W" with no subcellular location; that it is the cytoplasmic enzyme is established / uncited here (recalled), not
a database fact. `YALI1D19110g`'s annotation names both a *C. albicans* ProRS annotated cytoplasmic
([P78600](https://rest.uniprot.org/uniprotkb/P78600), location Cytoplasm, ECO:0000305) and YER087W / AIM10
([P39965](https://rest.uniprot.org/uniprotkb/P39965), "Probable proline--tRNA ligase, mitochondrial", location
Mitochondrion, ECO:0000269, PMID 16823961) as its ancestral-locus (synteny) partner. Both UniProt entries were fetched
by the independent audit this round (verified there, not by the main session). The first report summarised only the
YER087W part; that was a selective summary, corrected here. The localisation call therefore rests on the two
predictors, not on homology names. *Clarification (note-revision audit, same day):* the "synteny" assignment used
in this round and the review is the annotation's "ancestral locus Anc_x" tag, which is part of the same "similar to
*S. cerevisiae* …" clause; it is not an independent *Y. lipolytica* gene-order comparison (inferred; that
*Y. lipolytica* is not in YGOB is established / uncited).

Neither name is a verified *Y. lipolytica* name; the identities are database / homology assignments, and the
localisations are predictions. No *Y. lipolytica* localisation experiment was found. The decision rule was fixed before
the predictor run (prediction TASK §3) and was met. Not a reason: the TP gain. The model has no mitochondrial
prolyl-tRNA reaction (R673 is the only prolyl-tRNA synthetase reaction, checked), so `YALI1D19110g`'s own
essentiality cannot be modelled by this rule (category d of the 2026-10-27 study).

## 2. Implementation

| Part | File |
| --- | --- |
| Decision and evidence (data) | [model/curation/gpr/r673_prolyl_split.json](../../curation/gpr/r673_prolyl_split.json), written by [scripts/write_curation.py](scripts/write_curation.py) from [before.json](before.json) |
| Guarded step | `platform/scripts/gem_annotate/r673_prolyl_split.py`, `apply_r673_prolyl_split` |
| Flag | `--r673-prolyl-split` (`cli.py`), called in `main.py` after the R171 step; no flag dependency |
| Test | `platform/tests/test_r673_prolyl_split.py` (6 tests, 22 subtests: 13 stale-input cases, 5 scope cases, 4 knockout cases) |

The step is the single-reaction pattern of the R171 step (itself from `aars_thr_ser_split.py`): spec and the full
before-state of R673 and both genes pinned; stale, foreign or partly applied states raise `ValueError` with no edit;
the curated state returns `already_correct`; notes are restored on cobra context exit. The data generator asserts that
no note needs an XML entity (the R171 round's roundtrip lesson).

## 3. Builds (verified)

From `platform/` at HEAD `ff7fc62` plus the uncommitted R171 and R673 changes, `/opt/miniconda3/bin/python`;
commands in [logs/build_commands.txt](logs/build_commands.txt).

| Build | Seed | Output SHA-256 |
| --- | --- | --- |
| Control: `ca43ae0c` argv (flag off) | 0 | `ca43ae0c12e01c13c7b2ca4792b745d00d4a28ed950f8e401aece37b2a0843d0`, byte-identical to the baseline |
| Candidate: + `--r673-prolyl-split` | 0 | `6af03c405a691a76aa5857e23667dd9d34642b9c7d383f4ec17159966fb6f25a`, [candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_20261029.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_20261029.xml) |
| Candidate | 1 | `6af03c40…`, identical |

**File name.** The first candidate build failed at write time (`OSError 63`, file name too long): appending
`_r673_prolyl_split` to the line's name plus the builder's temporary prefix exceeded macOS's 255-byte limit; nothing was
written ([logs/build_candidate_seed0_attempt1_name_too_long.log](logs/build_candidate_seed0_attempt1_name_too_long.log)).
The candidate was rebuilt with the same argv to a shorter name that drops the middle flag tokens
(`gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt`); the full
flag list is in its `.build.json`. The seed-1 build to a scratch path gives the same bytes, so the output name does not
enter the file. Both candidate logs report "R673 prolyl split: applied"; the pre-existing libSBML `R_R304` notes error
appears in the control, the seed-1 log and the second seed-0 log.

## 4. Whole-model diff (verified)

`model_diff.py` (unchanged copy, `8a464ff8…`) vs `ca43ae0c` ([model_diff_vs_ca43ae0c.json](model_diff_vs_ca43ae0c.json)):
3 differences, all intended: R673 rule, R673 notes, `YALI1D19110g` reaction list (`R673` → none). Raw line diff: 11
changed lines. No stoichiometry, bound, formula or charge change.

## 5. Knockouts (verified)

Task A medium (`6eca70fd…`) + PO1f profile (`35307853…`), list `1e887f5a…`
([logs/sweep_inputs_sha256.txt](logs/sweep_inputs_sha256.txt)); Gurobi 13.0.2 at 1e-9, reset per knockout, targeted
check (`ko_compare_two.py`, unchanged copy, `bb3118a6…`); SHA-checked short copies; 2,164 LPs;
[ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv).

- WT 0.3116608216 on both.
- Growth changes above 1e-6: only `YALI1E05929g`, 0.31166 → 0.
- Class change: `YALI1E05929g` FN → TP at every cutoff; nothing else. `YALI1D19110g` stays FN (KO/WT 1.0).
- No failing point, no fresh re-solve, identical non-optimal sets (16 unresolved + 4 mtDNA), max knocked flux 0.

| Cutoff | ca43ae0c TP / FN / FP / TN | 6af03c40 TP / FN / FP / TN | Recall | Specificity | Precision | MCC |
| --- | --- | --- | --- | --- | --- | --- |
| 1 % | 96 / 219 / 35 / 699 | 97 / 218 / 35 / 699 | 0.305 → 0.308 | 0.952 | 0.733 → 0.735 | 0.356 → 0.360 |
| 5 % | 102 / 213 / 37 / 697 | 103 / 212 / 37 / 697 | 0.324 → 0.327 | 0.950 | 0.734 → 0.736 | 0.370 → 0.373 |
| 10 % | 106 / 209 / 42 / 692 | 107 / 208 / 42 / 692 | 0.337 → 0.340 | 0.943 | 0.716 → 0.718 | 0.368 → 0.371 |
| **15 %** | 123 / 192 / 42 / 692 | **124 / 191 / 42 / 692** | 0.390 → **0.394** | 0.943 | 0.745 → **0.747** | 0.420 → **0.423** |

16 unresolved at each cutoff. Positive coverage: 326 of 1,612 listed genes are model genes, 315 scored, 11 among the
unresolved; 1,286 listed positives are out of scope. Negatives are the user-defined class; mtDNA genes are out of
scope; the list was used during development, so the scores are a development reference.

## 6. Mechanism (verified, figure run)

With `YALI1E05929g` knocked out ([escher/](escher/README.md), 12 LPs): before, R673 stays open through
`YALI1D19110g` and every mapped flux equals WT (R673 0.0643); after, R673 is (0, 0) and growth is 0, because
Pro-tRNA(Pro) for protein synthesis has no other source.

## 7. pytest

Full suite, workflow-automation venv, before = the R171 round's after-run (code unchanged since, TASK §1), after =
[logs/pytest_after.txt](logs/pytest_after.txt); comparison [logs/pytest_compare.txt](logs/pytest_compare.txt): 330 → 336
passed, the same 19 failed and 20 errors, subtests 909 → 931. Test by test, the only differences are the 6 new tests,
all passing; no existing test changed outcome.

## 8. Acceptance (TASK.md §4)

| # | Criterion | Result |
| --- | --- | --- |
| 1 | Control reproduces `ca43ae0c` | met |
| 2 | Seeds identical | met |
| 3 | Diff = R673 rule and notes + `YALI1D19110g` link | met (3 differences) |
| 4 | Sweep = baseline + the one transition; WT, unresolved unchanged; no failing point | met |
| 5 | Only new tests differ | met (§7) |

## 9. Limits, untested, next gates

- Localisation is predicted (TargetP and DeepLoc share sequence features; not independent); TargetP's mTP call for
  `YALI1D19110g` is moderate (0.625). Dual targeting of one product is not detectable by these tools.
- The W29 and CLIB122 gene models of `YALI1E05929g` differ in start site (11 residues); the 517-aa variant was not run.
- The gain rests on the consensus list for `YALI1E05929g`'s essentiality; it is the benchmark list.
- **Model notes and docstring overstate slightly (audit finding, not rebuilt):** the R673 note
  `r673_prolyl_split_evidence` calls YHR020W "the cytoplasmic ProRS" and YER087W / AIM10 "the mitochondrial ProRS"
  (UniProt: "Putative …" with no location; "Probable …, mitochondrial"), and omits the *C. albicans* "cytoplasmic" hit
  in `YALI1D19110g`'s annotation; the step's docstring has the same wording. Fixing the note is a data change and a
  rebuild to a new file (as for the R171 note); left for a combined note revision or the user's call.
  **Fixed later the same day** in the curated data, with the user's approval:
  [../r171_r673_note_revision_20261029/](../r171_r673_note_revision_20261029/REPORT.md), candidate `25dabf62`. The step
  docstring was not overstated on re-reading, so it was not changed.
- `YALI1D19110g` stays FN: the missing mitochondrial Pro charging and the dead-end mitochondrial tRNAs are a separate
  line (category d).
- One medium. Not traced: a small zero-growth pFBA flux through R689 / R1225 after the knockout (escher/README.md).
- Provable: only the two R673 genes can change class under this edit. Statically verified: the diff. Observed: sweep
  and figure. Untested: other media.

## 10. Independent audit

A fresh read-only subagent audited both stages with its own scripts and 8 LPs ([AUDIT.md](AUDIT.md)). Supported: the
prediction stage (sequences, raw outputs, every score, the start-site claim, and that the decision rule was written
about a minute before the first result, from the session transcript), SHAs and builds (including the failed first
attempt leaving no file), the whole-model diff, the guard and test (6 passed, 22 subtests), the spot check and
recomputed matrices, the full pytest comparison, the figure numbers. Partially supported: the evidence statements
(the contrary *C. albicans* annotation hit omitted; YHR020W "cytoplasmic" and AIM10 without "probable" stated as fact).
Fixed: §1 now gives the full annotation strings and labels the ortholog locations; the same correction is in the
prediction REPORT; §9 records the note wording; §7 and acceptance #5 were filled. Not changed: the model note and the
step docstring (need a rebuild / would change the recorded source SHA); the DeepLoc encoder and TargetP settings are
evidenced by the form reads recorded in the prediction REPORT §2 and the session, not by a file in `raw/`; the stale
docstring of the copied `ko_compare_two.py`.
