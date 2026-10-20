# Restore of four reviewed GPRs cleared by the metadata selection

**New candidate:** `model/candidates/E5_coq9_resp_pkg_r132_gpr_restore_v2_20261019.xml`
(SHA-256 `c98961e960c6a0f8283f485e32daf20b04dd70d998fad86fc09dca0d14ccae4c`). It was built from
a56f3649's chain plus the new opt-in flag. It is not a reference model.
Scope and authorization: [TASK.md](TASK.md), including the corrections added after the audit.

**Superseded first build:** `E5_coq9_resp_pkg_r132_gpr_restore_20261019.xml` (`38e257b6…`).
It had the same optimization problem, but its R612 note overstated the URA3 evidence (§5). Its
three files (`.xml`, `.build.json`, `.coq9_genes.tsv`) were deleted on the user's instruction
(2026-10-19). That happened after the v1/v2 comparison and the solves below were done.
Its SHA stays in [logs/sha256.txt](logs/sha256.txt).

## 1. What changed and why

The 2026-09-10 metadata reaction selection (`metadata_reaction_fields_preserve_trna_20260910`)
runs after the builder's reviewed GPR patches. For four reactions it writes back the source
model's fields:
- an empty GPR on R612, R2274, R1302 and R1303;
- an R2274 equation without its 3 H₂O;
- an R1303 equation with an extra H⁺.

The 2026-10-19 orphan inventory ([report](../orphan_reactions_20261019/REPORT.md), group 1) found
that this left four orphan reactions that already had reviewed genes, and two mass-imbalanced
equations.

Changes in the new candidate:

| Reaction | GPR | Equation |
| --- | --- | --- |
| R612 orotidine-5'-P decarboxylase | `YALI1E31685g` (URA3) | unchanged |
| R2274 triolein hydrolysis (ex) | `YALI1A21372g` (LIP2) | 3 H₂O + triolein → glycerol + 3 oleate |
| R1302 pyridoxal kinase | `YALI1A08512g` | unchanged |
| R1303 pyridoxamine kinase | `YALI1A08512g` | ATP + pyridoxamine → ADP + pyridoxamine-5'-P (H⁺ removed) |

Each reaction also gets four notes: the change ID, the evidence, the evidence status and the
sources. The selection's `metadata_previous_*` notes stay as history. No gene, metabolite,
formula, charge or bound changed, and the selection record itself is unchanged.

## 2. Mechanism (code)

**Data:** [`model/curation/gpr/selection_cleared_gpr_restore.json`](../../curation/gpr/selection_cleared_gpr_restore.json).
For each reaction it holds:
- the cleared state (= the selection's "after");
- the restored state (= the selection's "before");
- the gene identity, evidence tier and status, and the sources;
- the formula, charge and compartment of every species involved.

**Step:** `apply_selection_gpr_restore()` in
`platform/scripts/gem_annotate/selection_gpr_restore.py`. Before it edits anything, it checks:
- the scope (change ID, status and exactly the four reactions);
- every species against the record;
- each restored equation balances, computed from the recorded formulas;
- each gene exists with the recorded name and identifiers;
- no conflicting notes are present (compared after XML unescaping, so a model saved and
  reloaded through SBML compares equal);
- each reaction's definition is exactly its cleared or its restored state.

Anything else raises `ValueError` before any edit. A rerun returns `already_correct`, including on
the saved and reloaded candidate.

**Wiring:**
- The flag `--restore-selection-cleared-gprs` (`cli.py`) is off by default.
- `main.py` calls `restore_after_selection()` after the R132 step. It refuses to run unless the
  metadata selection completed.
- The step writes one line to the build log and records its result under
  `reaction_selection.selection_gpr_restore` in the build record.

**Test:** `platform/tests/test_selection_gpr_restore.py` (6 tests, 13 subtests):
- restore, balance, KO propagation and idempotence;
- rerun after an SBML save and reload;
- ten stale-input cases that must raise with no edits;
- wrong scope and an unbalanced target;
- the flag is opt-in and needs a complete selection;
- on the real a56f3649 file, only the four targets change.

## 3. Evidence for the gene assignments (opened 2026-10-19; also opened by the auditor)

| Reaction | Gene | Evidence | Tier |
| --- | --- | --- | --- |
| R612 | YALI1E31685g, URA3 | [Swietalski et al. 2022, *Yeast* 39:230](https://pubmed.ncbi.nlm.nih.gov/34648204/) (abstract only; the Wiley full text returned 403 and it is not in PMC). In PO1f the URA3 deletion leads to orotic-acid accumulation and secretion, which the authors attribute to the missing orotidine-5'-phosphate decarboxylase. A URA3 rescue mutant completely prevents the secretion. NCBI Gene 90949730 → Swiss-Prot [Q12724](https://rest.uniprot.org/uniprotkb/Q12724) PYRF_YARLI (URA3, EC 4.1.1.23, reviewed, protein evidence by homology). The model's UniProt A0A1H6PUU4 is deleted from UniProtKB but has the same sequence (UniParc UPI0000132D4C). | A, verified |
| R2274 | YALI1A21372g, LIP2 | [Pignède et al. 2000, *J. Bacteriol.* 182:2802](https://pmc.ncbi.nlm.nih.gov/articles/PMC101989/) (full text): no extracellular lipase activity detected in the lip2 KO, and no halo on triolein plates. NCBI Gene 2905927 → UniProt [Q9P8F7](https://rest.uniprot.org/uniprotkb/Q9P8F7) lip2 (EC 3.1.1.3, unreviewed). | A, verified |
| R1302, R1303 | YALI1A08512g, name not verified | NCBI Gene 2905869 → UniProt [Q6CHH8](https://rest.uniprot.org/uniprotkb/Q6CHH8) "pyridoxal kinase" (EC 2.7.1.35, unreviewed, PE3). The model's A0A1D8N468 is deleted from UniProtKB but has the same sequence (UniParc UPI00003BB843). The same gene is already on R1306. | C, annotation only |

**Limits on the evidence:**
- **LIP2:** the activity assay used an olive-oil emulsion, and triolein appears only in the plate
  assay. The lip2 KO still grows on triglycerides, so LIP2 is not the only route to use them.
- **R1302/R1303:** substrate specificity for pyridoxal and pyridoxamine is not shown in
  *Y. lipolytica*. It is an inference from the EC class.

## 4. Verification

| Check | Result | Level |
| --- | --- | --- |
| Builder sources = a56f3649 build record | 0 files differ (before editing) | verified |
| Control builds, no new flag (before and after the post-audit fix) | Both SHA-256 `a56f3649…`, byte-identical to the baseline | verified |
| v2 candidate build | `c98961e9…`. Build record: input and sources unchanged during the build, build complete, `Selection GPR restore: applied` in the log | verified |
| Whole-model diff vs a56f3649 ([model_diff.json](model_diff.json)) | 13 differences: GPR and notes of the 4 reactions, the R1303 and R2274 equations, and the reaction membership of the 3 genes. Nothing else differs. The raw XML differs in 30 lines, all inside these reactions. `.coq9_genes.tsv` is identical. Counts are 2306 / 1875 / 1073. | verified |
| v2 vs v1 | Only the R612 and R1302 evidence notes differ. `model_definition` and `solver_definition` are equal. | verified |
| Balance | All 4 reactions return `check_mass_balance() == {}` | verified |
| Rerun on the saved v2 file | `already_correct` for all 4 | verified |
| pytest | Before: 165 passed, 11 failed, 14 errors, 31 skipped. Final: 171 passed, 11 failed, 14 errors, 31 skipped. Per test: identical except the 6 new tests, which pass ([logs](logs/)). The 25 failures and errors are pre-existing. | verified |
| Solves ([solves.csv](solves.csv)), GLPK, SD-Leu + PO1f | All 8 optimal. WT 1.430806 h⁻¹ on both models. KO/WT = 1.0 for YALI1E31685g, YALI1A21372g and YALI1A08512g on both models; the largest difference is 2.6e-14. These were solved on v1, which has the same solver definition as v2. The auditor reproduced the CSV byte for byte. | observed in runs |

**What the solves mean:**
- **URA3 KO:** it now switches R612 off. In PO1f, R612 is already closed by the strain profile, so
  there is no change.
- **LIP2 KO:** it now also switches off R2274. That matters only when triolein is supplied.
- **YALI1A08512g KO:** it now switches off R1302/R1303 as well as R1306. It was dispensable on
  SD-Leu before and still is.

TASK.md expected R1302 to stay blocked and R1303 to carry no flux. Those statements come from the
baseline flux analysis (FVA and pFBA) and were not re-run on the new model.

## 5. Independent audit

A fresh read-only reviewer was given TASK.md, the code, the data, the diff, the logs and the
models. It recomputed SHAs, re-derived the diff, re-ran the solves and the focused tests, and
opened every cited source.

| Claim | Verdict |
| --- | --- |
| 1 Control build byte-identical | supported |
| 2 Diff confined to the targets | supported (its own diff also lists the derived metabolite-side change: m10[C_cy] leaves R1303 and m214[C_ex] joins R2274) |
| 3 Restored state = selection "before"; balanced | supported |
| 4 Guard fails closed, idempotent, opt-in | **partially supported**: idempotence failed after an SBML round trip |
| 5 pytest unchanged except new tests | supported |
| 6 Solves unchanged | supported, reproduced |
| 7 Sources support the assignments | supported, with a wording correction for R612 |
| 8 R1302/R1303 labelled annotation-only everywhere | supported |

**Fixed after the audit:**
- **Notes round trip.** SBML reload returns `&apos;` in notes, so a rerun on the saved v1
  candidate raised "conflicting restore notes". Note comparison now unescapes both sides, and a
  save-reload-rerun test was added. The new test was confirmed to hit the escaped notes.
- **R612 wording.** "A URA3 rescue mutant restores it" became "completely prevents the
  secretion", which is what the abstract states.
- **Sequence identity.** The UniParc link for the two deleted UniProt accessions was added.
  Rebuilt as v2.
- **TASK.md.** The control path, the list of data files and the scope of the expected-effect
  statement were corrected in an appended section.

**Not changed (accepted as by design or recorded as limits):**
- A reaction whose restore notes are present but whose definition is back in the cleared state
  is restored again.
- There is no rollback if the edit loop fails halfway. This is unreachable after the
  pre-validation, and the build aborts before writing.
- The default (E0) build path with the flag is untested.

## 6. Limits and findings outside this change

- **LIP2 knockout on triolein (inferred, structural).** R2274 is the only reaction that consumes
  triolein[ex], and the triolein exchange R2273 is export-only. Under a medium with triolein
  uptake, a LIP2 knockout would block all triolein use. Pignède et al. report that the lip2 KO
  still grows on triglycerides, so the model would call LIP2 essential where the strain is not.
  No triolein condition was run. Fixing this needs evidence for the other lipases, and is not
  part of this change.
- **Selection re-adds two known mis-annotations (verified, outside this scope).**
  - The same 2026-09-10 selection writes `YALI1E07744g or YALI1D02625g` on R765/R766 and
    `YALI1F23441g or YALI1E11370g` on R671.
  - The builder removes YALI1E07744g (trehalase) and YALI1E11370g (GatB) from these reactions as
    mis-annotations before the selection runs.
  - Both a56f3649 and the new candidate carry them.
  - The selection also turns R190 from an AND of two genes into an OR of three, and widens R159
    to the same OR.
- **Stale strain-profile text.** The PO1f profile copy under `reference_pipeline_restore_20260909`
  still says R612 has no GPR. That is a recorded file and was not edited.
- **Pre-existing log line.** The build log line `LibSBML error code -3` ("Setting notes on
  R_R304") appears in every build, including the control.
- **Opt-in.** The default build (`b4ce0974…`) and all earlier candidates are unchanged. Making the
  restore the default, or folding it into the selection as an exclusion (as done for R_NTP1 on
  2026-09-14), is a separate decision.
- **Interpreter.** The builds and tests ran with the main checkout's `.venv` interpreter (Python
  3.13.5, cobra 0.30.0, optlang 1.8.3).
- **Not published.** Nothing is committed or pushed.
