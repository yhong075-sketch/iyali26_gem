# R1964 provisional GPR and cross-reference correction through the pipeline (2026-10-19)

## Authorization

The user, in this session, after reading
[r1964_gpr_review_20261019/REPORT.md](../r1964_gpr_review_20261019/REPORT.md) §7:

> Approve the R1964 plan: opt-in step with provisional GPR YALI1F03452g or YALI1E06616g, plus the
> MNXR198271/Rhea 17921/EC 1.1.1.11 annotation fix, and the 5 LPs.

Not authorized: commit, push, merge, the optional sixth LP (double KO), the optional evidence steps
(supplementary download, Liu 2026), any KO screen, any other model change.

## Starting state

- Worktree `claude/reverent-wilbur-4dc5b6`, HEAD `285eaee`. At start the only untracked items were the
  review folder `model/reports/r1964_gpr_review_20261019/` and this folder.
- Baseline model: a56f3649 `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`, SHA-256
  `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`.
  - It is not a designated reference model.
  - Recorded build options: `--offline --no-solve --coq9-curation qcycle --energy-candidate E5
    --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh`.
  - Recorded research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
    `IYALI26_RESEARCH_ROOT` is unset in this shell, so it is passed with `--research-root`.
- Interpreter: main checkout `.venv` (cobra 0.30.0, Gurobi academic licence).
- Target state in a56f3649: [before.json](before.json).
- Sibling work not used as a base: R2041 forward-only (`claude/bold-chaplygin-1c16df`, now committed
  as `ae54dd8`), GPR fixes (`claude/clever-kare-7d2e1a`), and the R2190 plan
  (`claude/clever-austin-f6a04e`). All three touch or plan to touch `main.py`/`cli.py`.

## Targets and change

- **Reaction R1964** "arabinitol dehydrogenase", D-xylulose + NADH ⇌ NAD + D-arabinitol, C_cy,
  [−1000, 1000].
  - GPR: "" → `YALI1F03452g or YALI1E06616g`, status `user_authorized_provisional_assignment`.
  - Cross-references: replace the L-xylulose record (bigg/vmh ABTD1, MNXR154181) with the D-xylulose
    reaction:
    - metanetx MNXR198271;
    - rhea 17921–17924;
    - kegg R05604;
    - ec-code 1.1.1.11;
    - bigg ABTDG;
    - vmh ABTDG (MNXref's VMH counterpart of the ABTD1 entry being replaced);
    - metacyc D-ARABINITOL-4-DEHYDROGENASE-RXN.

    The model's paired short keys (`biggr`, `vmhr`, `keggr`, `metacycr`, `rh`) follow the same
    values. SBO:0000176 is kept.
  - New notes record the curation, its evidence and its limits. `PROTEIN_CLASS 1.1.1.11` is kept.
  - Bounds, stoichiometry, species and name stay unchanged.
- **New gene YALI1E06616g** (GeneID 2912324, XP_503594.1, identical to YALI0E05643g/Q6C6W8). It
  gets annotations sbo SBO:0000243, ncbigene 2912324, kegg.genes yli:2912324, refseq XP_503594.1 and
  uniprot Q6C6W8. No gene name is set, because none is verified.
- **YALI1F03452g:** unchanged. It already exists, so it only gains R1964 as a second reaction.
- **Mechanism:**
  - data file `model/curation/gpr/r1964_gpr_assignment.json`;
  - guarded `apply_r1964_gpr_assignment` in `platform/scripts/gem_annotate/r1964_gpr.py`;
  - opt-in flag `--r1964-gpr`, applied in `main.py` after the `--r132-mito-adh` step;
  - focused test `platform/tests/test_r1964_gpr.py`.
- **Not touched:** every other reaction, metabolite and gene; `_apply_gpr_assignment`; the metadata
  selection.

## Expected effect

- WT growth is unchanged, because a GPR enters the LP only through knockouts (provable).
- KO of YALI1F03452g: R1964 stays open through YALI1E06616g, so KO/WT should equal M0's (inferred).
- KO of YALI1E06616g: KO/WT should be 1 (inferred).
- The whole-model diff should be R1964 (GPR, annotation, notes) plus one added gene.

## Acceptance criteria

1. The control build without the flag reproduces `a56f3649…` byte for byte.
2. The whole-model diff against a56f3649 contains only R1964's `gpr`, `annotation` and `notes`, the
   added gene YALI1E06616g, and the gene-reaction membership change of YALI1F03452g.
3. R1964's mass and charge residuals stay {} and 0.
4. The focused test passes. The full suite is unchanged test by test apart from the new tests.
5. The 5 LPs are optimal, and their results are within the expectations above.
6. An independent read-only audit finds no undeclared change.

## Compute budget and stop conditions

- **Builds:** 2 (candidate and control). The control is written to `.scratch/`.
- **LPs:** exactly 5, in the context of `coq9_respiratory_package_20261008/run_controls.py` (PO1f
  SD-Leu, the recorded solver configuration):
  1. WT on M0 = a56f3649;
  2. WT on M1 = the new candidate;
  3. KO YALI1F03452g on M0;
  4. KO YALI1F03452g on M1;
  5. KO YALI1E06616g on M1.
- **Stop conditions:**
  - WT M0 ≠ M1 beyond 1e-9 relative: stop, since that is impossible for a GPR-only diff;
  - a KO/WT that changes the call at any of the 1/5/10/15 % cutoffs: stop and report the R1964
    flux;
  - a non-optimal status: record the raw value and stop. It is not biological evidence and is never
    retried.
- **Outputs:**
  - this folder: `before.json`, `model_diff.json`, `lp_results.json`, logs, scripts and REPORT.md;
  - the new candidate in `model/candidates/`;
  - a `model/STATE.md` entry. No commit.
