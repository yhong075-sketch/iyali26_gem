# Restore four reviewed GPRs cleared by the metadata selection (2026-10-19)

## Authorization (user, this session)

> Plan a curation change that re-applies the R612, R2274, R1302 and R1303 GPRs (and balanced
> R1303/R2274 stoichiometry) after the metadata selection step. Show me the plan before editing
> anything.

> Approved: implement the plan for R612, R2274, R1302 and R1303 as described.

The approved plan: an opt-in builder flag applied at the end of the build, a curated data file, a
guarded `apply_*` step, a focused test, a control build and a candidate build, a whole-model diff,
8 LP solves, an independent audit, then the write-up. No commit or push.

## Baseline

- Worktree branch `claude/orphan-reaction-newest-model-b8ddc8` at `4d97b828ac2b0dc95d592060605af75069705c09`.
- Baseline model: `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`,
  SHA-256 `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`.
  This is the newest candidate, not a designated reference.
- Builder code identity: all `platform/scripts/gem_annotate/*.py` hashes equal the
  `source_sha256` recorded in `E5_coq9_resp_pkg_r132_20261012.build.json` (checked 2026-10-19).
- Data added since that build: `model/conditions/assumptions/r825_reversible_20261013.json`.
  The build hashes it but does not use it.
- State of every target: [before.json](before.json).

## Targets and intended change

| Reaction | Now (selection "after") | Target (selection "before") |
| --- | --- | --- |
| R612 | GPR `""` | GPR `YALI1E31685g` |
| R2274 | GPR `""`; triolein → glycerol + 3 oleate (short by H6 O3) | GPR `YALI1A21372g`; add 3 `m214[C_ex]` H₂O as substrate |
| R1302 | GPR `""` | GPR `YALI1A08512g` |
| R1303 | GPR `""`; + `m10[C_cy]` H⁺ product (H +1, charge +1) | GPR `YALI1A08512g`; remove the H⁺ |

- No gene, metabolite, formula, charge or bound change.
- Each target equation balances for mass and charge with the current formulas (checked in memory
  2026-10-19).
- Notes: add the restore record; keep the existing `metadata_previous_*` notes as history.
- The 2026-09-10 selection record (`metadata_reaction_selection.json`) is not modified.

**Expected effect (inferred, to be checked by solves).** No growth change on SD-Leu + PO1f:
- R612 is closed by the PO1f overlay;
- R2274 is blocked because there is no triolein;
- R1302 is blocked;
- R1303 carries 0 at the optimum.

## Acceptance criteria

1. A control build without the new flag (a56f3649's exact flags) reproduces `a56f3649…`.
2. A build with the new flag differs from a56f3649 only in these 4 reactions (GPR, stoichiometry,
   notes) and the gene–reaction links of the 3 genes.
3. All 4 reactions are mass and charge balanced.
4. Guard behaviour, shown by tests:
   - a rerun is a no-op;
   - a stale or unexpected state raises an error without edits;
   - the flag is rejected without the selection having completed.
5. `python -m pytest`: same results as before, plus the new tests passing.
6. Solves are recorded with status, and any growth or KO change is reported.

## Compute budget and stop condition

- **Builds:** 2 (control to `platform/.scratch/`, candidate to `model/candidates/`), offline, no solve.
  - Research root: `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, as recorded in
    a56f3649's build.
- **Solves:** 8 LPs with GLPK on SD-Leu + PO1f.
  - Models: baseline and new.
  - Per model: WT growth and single-gene KO of YALI1E31685g, YALI1A21372g and YALI1A08512g.
- **Stop:**
  - the control build ≠ a56f3649;
  - any non-target diff;
  - a balance failure;
  - a test regression;
  - a source that contradicts the assignment.

## Evidence (checked 2026-10-19)

- **R612 URA3:**
  - Swietalski et al. 2022, *Yeast* 39:230, [doi:10.1002/yea.3673](https://doi.org/10.1002/yea.3673),
    [PMID 34648204](https://pubmed.ncbi.nlm.nih.gov/34648204/). Abstract read via PubMed (the Wiley
    full text returned 403): the PO1f URA3 deletion removes orotidine-5'-phosphate decarboxylase,
    and a URA3 rescue mutant prevents orotic-acid secretion.
  - Gene identity: NCBI Gene 90949730 (YALI1_E31685g) cross-references UniProt
    [Q12724](https://rest.uniprot.org/uniprotkb/Q12724) PYRF_YARLI (reviewed; URA3, EC 4.1.1.23).
- **R2274 LIP2:**
  - Pignède et al. 2000, *J. Bacteriol.* 182:2802, [PMC101989](https://pmc.ncbi.nlm.nih.gov/articles/PMC101989/)
    (full text read): the lip2 KO strain has no extracellular lipase activity, and the KO gives no
    halo on triolein plates. Activity was assayed on an olive-oil emulsion.
  - The KO still grows on triglycerides, so another route for triglyceride use exists.
  - Gene identity: NCBI Gene 2905927 (YALI1_A21372g) cross-references UniProt
    [Q9P8F7](https://rest.uniprot.org/uniprotkb/Q9P8F7) lip2 (unreviewed, EC 3.1.1.3).
- **R1302, R1303 YALI1A08512g:**
  - Annotation level only. NCBI Gene 2905869 cross-references UniProt
    [Q6CHH8](https://rest.uniprot.org/uniprotkb/Q6CHH8) "pyridoxal kinase" (unreviewed, EC 2.7.1.35).
  - The model gene carries `uniprot A0A1D8N468` and EC 2.7.1.35, and is already on R1306.
  - No experimental evidence in *Y. lipolytica*.

## Corrections after the independent audit (2026-10-19)

The text above is kept as written. These corrections apply:

- **Control output path.** The control build was written to the repo-root `.scratch/` (git-ignored),
  not to `platform/.scratch/`.
- **Data added since the baseline build.** The baseline build's `data_sha256` lacks five files that
  the current builds hash. Four are condition files:
  - `assumptions/r825_reversible_20261013.json`
  - `assumptions/w29_13c_upper_caps_20261012.json`
  - `media/po1f_csm_leu_w29glycerol_20261012.csv`
  - `media/po1f_csm_leu_w29oleate_20261012.csv`

  The fifth is this change's curation file. `model/curation/README.md` also differs. The builder
  does not read any of these, and the control build is byte-identical to a56f3649.
- **Expected effect.** The 8 solves check only WT growth and the KO ratios. "R1302 blocked" and
  "R1303 carries 0 at the optimum" come from the orphan report's baseline FVA and pFBA and were not
  re-run on the new model.
- **R612 evidence wording.** The first candidate's note said a URA3 rescue mutant "restores" the
  activity. The abstract supports only that it prevents orotic-acid secretion. The data file was
  corrected and the candidate rebuilt (see REPORT.md).
- **Rerun after a save.** In the first implementation, rerunning the step on a saved and reloaded
  model raised an error instead of doing nothing, because SBML returns `&apos;` in notes. The step
  now compares unescaped note text, and a save-reload-rerun test was added.
