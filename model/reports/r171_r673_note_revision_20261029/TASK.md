# R171 and R673 note wording revised in the curated data, rebuilt (2026-10-29)

## 0. Authorization

**User, 2026-10-29**, answering "Fix the R171 and R673 note wording in one rebuild now?": "Yes, fix notes" (reply
text: fix the R171 note "catalytically active" → "hydrolytically active" and the R673 note (YHR020W putative, AIM10
probable, add the *C. albicans* cytoplasmic hit) through the curated data, rebuild to a new candidate, verify by
whole-model diff that only those notes differ; plan first, no commit). Then, answering "Approve the note-revision
plan with option A?": **"Approve, option A"** (reply text: "Approved: run the note-revision plan with option A (edit
the two curated data files in place, rebuild to a new candidate, mark ca43ae0c and 6af03c40 superseded). No commit.").

## 1. Why

Both independent audits found note text that overstates the evidence:
- R171 ([../r171_atp_synthase_gpr_20261029/AUDIT.md](../r171_atp_synthase_gpr_20261029/AUDIT.md) §7): the monomer of
  Hahn 2016 was measured hydrolysing ATP, not synthesising it; "catalytically active" overstates.
- R673 ([../r673_prolyl_split_20261029/AUDIT.md](../r673_prolyl_split_20261029/AUDIT.md) §7): YHR020W is "Putative
  proline--tRNA ligase" with no UniProt location; AIM10 is "Probable … mitochondrial"; `YALI1D19110g`'s annotation also
  names *C. albicans* P78600, annotated cytoplasmic.

UniProt records re-fetched this session (verified; [evidence/](evidence/)): P38708 "Putative proline--tRNA ligase
YHR020W", no location; P39965 "Probable proline--tRNA ligase, mitochondrial", Mitochondrion (ECO:0000269, PMID
16823961); P78600 *Candida albicans* "Proline--tRNA ligase, cytoplasmic", Cytoplasm (ECO:0000305).

## 2. Starting state

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-r171-pro-7e2477` @ `ff7fc62`, uncommitted R171 and R673 rounds |
| Baseline for the diff | `6af03c405a691a76aa5857e23667dd9d34642b9c7d383f4ec17159966fb6f25a` (R673 round) |
| Data before | `r171_atp_synthase_gpr.json` `63a7c6ba…`, `r673_prolyl_split.json` `2446217a…`, step `r171_atp_synthase_gpr.py` `315646ab…`; copies in [before/](before/) |

## 3. The edits (data only, plus one docstring)

| Where | Before | After |
| --- | --- | --- |
| R171 note `r171_atp_synthase_gpr_evidence` | "…the Y. lipolytica monomer lacking them is catalytically active (Hahn 2016, PMID 27373333)…" | "…the Y. lipolytica monomer lacking them hydrolyses ATP at about the rate of the dimer (ATP hydrolysis, not synthesis, was measured; Hahn 2016, PMID 27373333)…" |
| R673 note `r673_prolyl_split_evidence` | "annotation similar to S. cerevisiae YHR020W, the cytoplasmic ProRS" / "annotation similar to YER087W AIM10, the mitochondrial ProRS" | "annotation similar to S. cerevisiae YHR020W, a putative prolyl-tRNA synthetase with no UniProt location, ancestral locus Anc_1.352" / "annotation similar to C. albicans P78600, a prolyl-tRNA synthetase annotated cytoplasmic, and to S. cerevisiae YER087W AIM10, a probable mitochondrial prolyl-tRNA synthetase, ancestral locus Anc_7.368"; plus "Homology names are mixed; the assignment rests on the two predictors and synteny." |
| R673 data `evidence` list (not in the model) | P38708 "cytoplasmic prolyl-tRNA synthetase"; P39965 "mitochondrial prolyl-tRNA synthetase" | UniProt wording as in §1; a P78600 entry added |
| `r171_atp_synthase_gpr.py` docstring (not in the model) | "subunits g and k, which cannot synthesise ATP alone" | "subunits g and k, with no evidence that either synthesises ATP alone" |

Both edits go through each round's `scripts/write_curation.py`, re-run to regenerate the JSON. Under option A,
rebuilding the old flags no longer reproduces `ca43ae0c` or `6af03c40`; both files stay on disk and are marked
superseded (note wording) in `candidates/README.md`.

## 4. Acceptance criteria

1. Focused R171 and R673 tests pass; full pytest: no test changes outcome versus the R673 round's after-run.
2. New candidate seeds 0 and 1 identical.
3. Whole-model diff vs `6af03c40`: exactly 2 differences, R171 notes and R673 notes; within them only
   `r171_atp_synthase_gpr_evidence` and `r673_prolyl_split_evidence` change, and the new values equal §3.
4. Sanity solves on the new candidate (Task A + PO1f, Gurobi 1e-9, basis reset): WT 0.3116608216, `YALI1F04433g` KO/WT
   0.1097, `YALI1E05929g` growth 0. No sweep: a notes-only diff leaves the optimisation problem unchanged.

## 5. Plan and budget

Edit, tests, build (`…_r673_prolyl_split_notes_v2_20261029.xml`) and seed check, diff, 3 sanity LPs, short independent
audit, REPORT, STATE, README rows. No Escher figure (no flux change). Budget: 2 builds, 3 LPs (+ at most 10 for the
audit), one pytest run. Stop if the diff shows anything else. No commit.

## 6. Second pass (authorized 2026-10-29)

**User, 2026-10-29**, answering "Apply the two remaining wording fixes in one more rebuild?": "Yes, apply both" (reply
text: YHR020W "no UniProt subcellular-location annotation", and replace "rests on the two predictors and synteny" with
"rests on the two predictors; the ancestral-locus tags come from the same homology annotation" (also in the R673
evidence_tier note); rebuild, diff, no commit; plan first). Then, answering "Approve the second-pass wording plan?":
**"Approve as written"** (reply text: "Approved: run the second-pass wording plan as written (R673 data only, v3
candidate, notes-only diff vs 25dabf62, mark 25dabf62 superseded). No commit.").

Why: this round's audit (AUDIT.md §4) found "no UniProt location" exact only for the Subcellular location section
(P38708 carries GO cytoplasm, IEA), and "synteny" not independent of the homology names (the "ancestral locus Anc_x"
tags are part of the same annotation clause; inferred).

| Field (R673 data `r673_prolyl_split.json`, before `e1d151a6…`, copy in [pass2/before/](pass2/before/)) | Now | New |
| --- | --- | --- |
| note `r673_prolyl_split_evidence` | "…putative prolyl-tRNA synthetase with no UniProt location…" | "…with no UniProt subcellular-location annotation…" |
| same note | "Homology names are mixed; the assignment rests on the two predictors and synteny." | "Homology names are mixed; the assignment rests on the two predictors; the ancestral-locus tags come from the same homology annotation." |
| note `r673_prolyl_split_evidence_tier` and top-level `evidence_tier` | "…on the W29 sequences plus genome-annotation synteny; …" | "…on the W29 sequences plus genome-annotation homology names (mixed; the ancestral-locus tags come from the same annotation); …" |

Acceptance: only those fields change in the data; focused R673 test passes and full pytest is unchanged test by test;
v3 seeds 0 and 1 identical; whole-model diff vs `25dabf62` = R673 notes only, keys `r673_prolyl_split_evidence` and
`r673_prolyl_split_evidence_tier`; sanity solves as before. Output
`candidates/…_r673_prolyl_split_notes_v3_20261029.xml`; `25dabf62` marked superseded. Budget 2 builds, 3 LPs, one
pytest run, a short audit. No commit.
