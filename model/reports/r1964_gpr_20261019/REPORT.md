# R1964 provisional GPR and D-xylulose cross-references: implementation report

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). The evidence review
behind the change is [r1964_gpr_review_20261019/REPORT.md](../r1964_gpr_review_20261019/REPORT.md).
Commands are in [commands.txt](commands.txt).

Tiers: **verified** (run or checked in this session, output cited), **provable** (follows from
verified facts), **inferred** (reasoned, no direct check), **recorded** (from an earlier record, not
re-run).

## 1. Summary

- **New candidate (verified).**
  `model/candidates/E5_coq9_resp_pkg_r132_r1964_gpr_20261019.xml`, SHA-256
  `db42974360491a33e3366e5e2a0379c6712143a836432dbb1dbf36facadebf78`.
  - It is a56f3649 plus the opt-in step `--r1964-gpr`.
  - It is a candidate, not a reference model. Uncommitted.
- **What changed (verified, whole-model diff, 5 differences):**
  - R1964 GPR: "" → `YALI1F03452g or YALI1E06616g`.
  - R1964 cross-references: L-xylulose ABTD1/MNXR154181 → the D-xylulose reaction (MNXR198271,
    Rhea 17921–17924, KEGG R05604, EC 1.1.1.11, BiGG/VMH ABTDG, MetaCyc).
  - R1964 notes: 10 new keys. `PROTEIN_CLASS 1.1.1.11` is kept.
  - New gene YALI1E06616g.
  - YALI1F03452g gains R1964 as a second reaction (R2190, R1964).

  Bounds, stoichiometry, species and every other reaction, metabolite and gene are unchanged. R1964
  stays balanced (residual {}, verified).
- **Control (verified).** The same build without the flag reproduces `a56f3649…` byte for byte.
- **LPs (verified; 5 of 5 authorized, all optimal):**
  - WT is 1.430805735352851 h⁻¹ on both models; the relative difference is 0.
  - KO/WT = 1.000 for YALI1F03452g on M0 and M1, and for YALI1E06616g on M1. No call changes at
    1/5/10/15 %.
  - R1964 and R2190 flux is 0 in all five solutions.
- **Audit (§4).** C1–C6 and C8–C10 are supported. C7 (guard) is partially supported: non-spec notes,
  part of YALI1F03452g's identity and YALI1E06616g's name were unchecked. The guard was tightened
  and the test strengthened. A verification rebuild of the candidate and the control is
  byte-identical (§3.1).
- **Tests (verified):**
  - the focused test passes (5 tests, 23 subtests, after the audit fix);
  - full suite, before and after the fix: unchanged test by test, apart from the 5 new passing tests
    (165 → 170 passed; the 11 failures, 14 errors and 31 skips are pre-existing and identical).
- **Evidence status of the assignment: provisional.**
  - In vitro activity of the purified proteins: verified from Wang et al. 2020.
  - The reaction products are assigned, not identified.
  - The C-4 (this reaction) versus C-2 (D-ribulose) reading is inferred and unsettled.
  - The GPR is not experimental validation of either gene's in vivo role.

## 2. Change and mechanism

| Item | Location |
| --- | --- |
| Curated data (decision, evidence, before/after state) | [`model/curation/gpr/r1964_gpr_assignment.json`](../../curation/gpr/r1964_gpr_assignment.json) (`R1964-ARADH-GPR-20261019`). SHA-256 `36b6a730…` after the audit fix; the first build used `5b86a2ba…` |
| Guarded step | `platform/scripts/gem_annotate/r1964_gpr.py`, `apply_r1964_gpr_assignment` |
| Wiring | `cli.py` flag `--r1964-gpr`; `main.py` calls the step after `--r132-mito-adh` and before the CoQ functional/C5/literature steps |
| Focused test | `platform/tests/test_r1964_gpr.py` |
| READMEs | `model/curation/README.md`, new section "Opt-in R1964 provisional GPR …"; `model/candidates/README.md`, new candidate row |
| Code diff | [logs/code_diff.patch](logs/code_diff.patch) |

**What the step checks before writing (fail closed, no partial edits):**
- the spec's schema, ID, status, reaction ID, the before GPR "" and the after GPR, and the exact gene
  sets;
- R1964's name, bounds, exact stoichiometry and species (formula, charge, compartment);
- that the (GPR, annotation) pair is exactly the before-state or the after-state;
- no conflicting notes, and no after-state without this curation's notes;
- R1964's notes outside this curation equal `before_notes` (`PROTEIN_CLASS 1.1.1.11`). Added after
  the audit;
- YALI1F03452g: its full annotation equals the recorded one (SBO, NCBI, UniProt, KEGG, RefSeq, EC),
  and its only other reaction is R2190. Before the audit, only NCBI, RefSeq and KEGG were checked;
- YALI1E06616g: absent, or present with exactly the spec identifiers, no other reaction, and no name
  other than the ID-derived default. The name check was added after the audit.

It then:
- sets the GPR, which makes cobra create YALI1E06616g;
- writes that gene's identifiers;
- replaces R1964's annotation dictionary with the spec's after-state;
- merges the notes.

A rerun returns `already_correct`.

**Why the step runs after the annotation stages.** It runs at `main.py:701`, after `annotate_reactions`
(`main.py:239`), `annotate_remaining_reactions` (`main.py:288`) and `backfill_reaction_xrefs`
(`main.py:318`), with line numbers as edited. The
gene-EC route therefore cannot re-attach YALI1F03452g's EC 1.1.1.250 (D-ribulose) identifiers to
R1964 (verified: the built R1964 carries only the spec's after-annotation).

**Gene name.** YALI1E06616g has no verified name, so the step sets none. The SBML writer stores the
SBML ID as `fbc:name` (`G_YALI1E06616g`), as for the 15 other unnamed genes in the model
(verified).

**Note text.** The first focused-test run failed its SBML-roundtrip assertion. The cause was the
apostrophes in one note, which the reader returns as `&apos;`. The note text was reworded before
the build; no apostrophes remain.

## 3. Verification

### 3.1 Builds (verified)

| Build | Command (from `platform/`) | Output SHA-256 |
| --- | --- | --- |
| Control | recorded a56f3649 options, no flag → `.scratch/r1964_gpr_20261019/control_no_flag.xml` | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` (= baseline) |
| Candidate | same options + `--r1964-gpr` → `model/candidates/E5_coq9_resp_pkg_r132_r1964_gpr_20261019.xml` | `db42974360491a33e3366e5e2a0379c6712143a836432dbb1dbf36facadebf78` |
| Control v2 (after the audit fix) | as Control → `.scratch/…/control_no_flag_v2.xml` | `a56f3649…` (= baseline) |
| Candidate v2 (after the audit fix) | as Candidate → `.scratch/…/candidate_v2.xml` | `db429743…` (= candidate) |

The TASK.md budget named 2 builds. The 2 v2 builds were added to verify the audit fix; they are
builds without solves, not LPs.

The candidate's `.build.json` (`6019fde5…`) records:
- `r1964_gpr: true`, HEAD `285eaee`, `input_unchanged: true`, `source_unchanged_during_build: true`;
- the step result `applied`, with `added_genes: ["YALI1E06616g"]`;
- counts of 2306 reactions, 1875 species and 1074 genes (a56f3649 has 1073 genes).

The `.coq9_genes.tsv` sidecar is identical to a56f3649's. Both builds print the libSBML
"Setting notes on … R304" message seen in earlier builds, such as the R2041 build log. It is not
related to this change.

**Provenance caveat.**
- The curation README section was written while the candidate build was running, after the build had
  recorded its input hashes. The `.build.json` therefore holds the pre-edit README hash
  (`1ff63ddf…`) and lists `M model/curation/README.md` in its end-of-build dirty status.
- The builder hashes every file under `model/curation/` into `data_sha256` (`main.py:761`), the
  README included, but does not use the README's content. The curated data file's recorded hash
  equals the file as it was at that build (`5b86a2ba…`). The audit fix later changed the data file
  (§4); the v2 rebuild records the new hash.

### 3.2 Whole-model diff (verified)

`scripts/model_diff.py` compares every reaction, metabolite and gene field, gene–reaction
membership, groups (name and members), compartments, the objective direction, each reaction's
objective coefficient, and model notes. Output: [model_diff.json](model_diff.json). The auditor's own
independent comparison also covered the objective expression, group kind and the medium, and found
the same five differences.

| Section | Item | Field | Before → after |
| --- | --- | --- | --- |
| reactions | R1964 | gpr | "" → `YALI1F03452g or YALI1E06616g` |
| reactions | R1964 | annotation | ABTD1 / MNXR154181 set → the D-xylulose set (TASK.md) |
| reactions | R1964 | notes | `PROTEIN_CLASS` → `PROTEIN_CLASS` + 10 curation keys |
| genes | YALI1E06616g | — | added (ncbigene 2912324, kegg.genes yli:2912324, refseq XP_503594.1, uniprot Q6C6W8, SBO:0000243) |
| genes | YALI1F03452g | reactions | [R2190] → [R1964, R2190] |

The raw XML diff ([logs/xml_diff.txt](logs/xml_diff.txt), 4 hunks) touches only R1964's notes,
annotation and `geneProductAssociation`, plus the new `geneProduct`.
[before.json](before.json) and [after.json](after.json) hold both states; the R1964 residual is {}
in both.

### 3.3 LPs (verified; [lp_results.json](lp_results.json), [scripts/run_lps.py](scripts/run_lps.py))

- **Context:** the protocol of `coq9_respiratory_package_20261008/run_controls.py`, through
  `coq_candidate_checks_20261005/manifest.json`:
  - SD-Leu medium with the `po1f_sd_leu_accrispr_v1` strain overlay;
  - Gurobi with Threads 1, Presolve 0, tolerances 1e-7 and a 60 s limit.
- **KO method:** `Gene.knock_out()` inside a model context.
- **Calls:** 5 of 5 budgeted.

The first column is TASK.md's LP number; the second is the execution order.

| TASK # | Run | LP | Status | Growth (h⁻¹) | KO/WT | Closed by the KO | R1964 / R2190 flux |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 1 | M0 WT | optimal | 1.430805735352851 | — | — | 0 / 0 |
| 3 | 2 | M0 KO YALI1F03452g | optimal | 1.430805735352848 | 1.000 | R2190 | 0 / 0 |
| 2 | 3 | M1 WT | optimal | 1.430805735352851 | — | — | 0 / 0 |
| 4 | 4 | M1 KO YALI1F03452g | optimal | 1.430805735352848 | 1.000 | R2190 (R1964 stays open via YALI1E06616g) | 0 / 0 |
| 5 | 5 | M1 KO YALI1E06616g | optimal | 1.430805735352848 | 1.000 | none | 0 / 0 |

- **WT:** M0 = M1 exactly. This is provable for a GPR-only diff and observed here.
- **KO/WT:** 0.999999999999998 (a 3e-15 h⁻¹ difference, solver round-off) in the three KO rows.
- **Stop conditions (audit Defect 4).** The script enforced only the non-optimal stop inside its
  loop. The WT-difference and call-change stops were evaluated after the run, from
  `lp_results.json`. Neither would have fired: the WT relative difference is 0, and `essential_at`
  is false at every cutoff for all three KOs. The script was not changed or rerun, because its SHA
  is recorded in `lp_results.json` and the LP budget is used up.
- **Benchmark:** neither gene is in `consensus_essential_genes.csv` (review §6.1). Both stay
  non-essential in the user-defined negative class, so the confusion matrix cannot change through
  these genes. That is inferred from the three KO rows; no full screen was run.
- **Not run:** the double KO (not authorized). By the OR rule it closes R1964 and R2190; R1964 is 0
  in every solution here.

### 3.4 Tests (verified)

- Focused: `tests/test_r1964_gpr.py`, 5 passed, 23 subtests after the audit fix (19 before). It
  covers:
  - 17 stale-input cases that fail closed without edits, each checked against its expected error
    message;
  - a check that no subcase leaks into the shared model;
  - 3 out-of-scope specs;
  - an exact in-memory change, idempotence, and an exact SBML roundtrip;
  - the OR logic under `knock_out` (no LP);
  - that the flag is opt-in.
- Full suite before / after ([logs/pytest_comparison.txt](logs/pytest_comparison.txt)):
  - before (HEAD 285eaee, no code change): 165 passed, 11 failed, 14 errors, 31 skipped,
    388 subtests;
  - after: 170 passed, 11 failed, 14 errors, 31 skipped, 407 subtests.
  - The only per-test changes are the 5 new `test_r1964_gpr` tests (absent → passed). Every
    pre-existing failure and error has the same outcome before and after.

## 4. Independent audit

A fresh read-only reviewer audited the deliverables against TASK.md. It had no conversation
history. It was not allowed to build, solve or edit; it could run the focused test and query
NCBI/UniProt/KEGG/Rhea/BiGG.

| # | Claim | Verdict |
| --- | --- | --- |
| C1 | Control build = a56f3649 byte for byte | supported |
| C2 | Candidate SHA and build record | supported |
| C3 | Independent whole-model comparison: only the 5 declared differences | supported |
| C4 | R1964 after-state equals the data file; residual empty | supported |
| C5 | YALI1E06616g identifiers and identity (2912324 = XP_503594.1 = Q6C6W8) | supported |
| C6 | Cross-references describe D-arabinitol ⇌ D-xylulose; ABTD1 uses L-xylulose | supported. VMH ABTDG was checked only through MNXref (VMH API unreachable); MetaCyc through MNXref and BiGG |
| C7 | Guard fails closed and is idempotent | **partially supported** (defects 1–3 below) |
| C8 | Flag opt-in; no other new code path | supported |
| C9 | 5 optimal LPs, protocol, budget | supported, with a protocol gap (defect 4) |
| C10 | pytest unchanged apart from the new tests | supported |
| C11 | Report statements consistent with the files | mostly supported (defects 5–8) |
| C12 | Only declared files changed | supported, but `model/candidates/README.md` was undeclared (defect 5) |

**Defects and what was done:**

| # | Defect | Action |
| --- | --- | --- |
| 1 | Notes outside the spec keys were not checked: a stale `PROTEIN_CLASS` or an extra note still gave `applied` | **Fixed.** `before_notes` added to the data file, a new guard check, and test cases `protein_class` and `extra_note` |
| 2 | YALI1F03452g identity was only partly checked (NCBI, RefSeq, KEGG) | **Fixed.** The full recorded annotation is now required to match; test case `existing_gene_uniprot` |
| 3 | The name of a pre-existing YALI1E06616g was not checked | **Fixed.** The name must be empty or `G_YALI1E06616g`; test case `new_gene_name` |
| 4 | The LP script did not enforce the WT-difference and call-change stops inside its loop | Reported (§3.3). Both were checked after the run and would not have fired; the LPs were not rerun |
| 5 | The `model/candidates/README.md` row was undeclared | Fixed in §2 |
| 6 | "The builder does not read the README" is inaccurate: the README is hashed | Fixed in §3.1 |
| 7 | The `model_diff.py` coverage was overstated (objective, groups) | Fixed in §3.2 |
| 8 | LP row numbers did not match TASK.md | Fixed in §3.3 |
| 9 | The STATE.md entry was still missing | Written after this audit |
| 10 | In the data file, the MNXR198271 URL was tagged `verified` although only the local tables were read | **Fixed.** The evidence entry now cites the local MNXref files and says the URL was not opened |
| 11 | The auditor created `platform/scripts/gem_annotate/__pycache__/` (4 `.pyc` files) | Left in place. It is gitignored and changes no SHA, test or record. Deleting files needs the user's OK |

**A further test defect found while fixing.**
- `Model.copy()` shares `Gene.annotation` dicts with the original (verified). The first version of
  the test mutated YALI1F03452g's annotation in place in the `existing_gene_identity` case. That
  edit leaked into the class-level model, so the later gene subcases may have raised for the wrong
  reason.
- The test now rebinds instead of mutating, asserts each subcase's expected error message, and
  checks that the shared model is unchanged.
- In-place `Gene.annotation` edits on model copies also appear in other tests, for example
  `test_r1025_gpr_assignment.py:55`, `test_r1026_gpr_assignment.py:55` and
  `test_r153_gpr_assignment.py:59` (grep). Whether they leak was not examined; that affects
  another line.

**Re-verification after the fixes (verified):**
- the focused test passes with 5 tests and 23 subtests;
- the v2 builds from the fixed code and data reproduce both earlier outputs byte for byte:
  - control `.scratch/r1964_gpr_20261019/control_no_flag_v2.xml` = `a56f3649…`;
  - candidate `.scratch/r1964_gpr_20261019/candidate_v2.xml` = `db429743…`;
- the v2 build record ([logs/candidate_v2.build.json](logs/candidate_v2.build.json)) holds the fixed
  data file `36b6a730…`, `r1964_gpr.py` `73835a03…` and the edited curation README `0bff5e41…`;
- the full suite after the fix ([logs/pytest_comparison_v2.txt](logs/pytest_comparison_v2.txt)):
  170 passed, 11 failed, 14 errors, 31 skipped, 411 subtests. The only per-test changes from the
  pre-change run are the 5 new tests (absent → passed);
- the post-fix code diff is [logs/code_diff_v2.patch](logs/code_diff_v2.patch). The audited version
  is `logs/code_diff.patch`.

Neither the guard change nor the data-file change alters the built model, so the candidate file
from the first build stands. Its `.build.json` records the pre-fix code and data hashes. The v2
build record (`logs/candidate_v2.build.json`) records the fixed ones.

## 5. Limits and what was not tested

- **Science (unchanged by this build).**
  - The products of YlAraDH1/YlAraDH2 were assigned, not identified.
  - D-ribulose was never tested.
  - The regiochemistry reading is inferred (review §5).
  - The GPR does not show which gene, if either, carries D-arabitol metabolism in vivo; Wang 2020
    found no deletion phenotype on glucose.
- **No flux consequence in this context (provable, review §6.2; observed).** On a56f3649, R1964
  carries flux only in the R803/R714/R2042/R2041/R2190 loop. After either loop fix (R2041
  forward-only, or R2190 → D-ribulose) it is blocked unless a medium opens D-xylose, D-arabinose or
  D-arabitol uptake. The new GPR matters for KO mapping and for future pentose media.
- **Cross-references.**
  - KEGG R05604 sits on MNXR198270 in MNXref 4.5, the α-D-xylulofuranose form; the model's
    MNXR198271 uses open-chain D-xylulose (`m1015`).
  - SEED rxn03883, also on MNXR198270, was not added because it was not in the approved list.
  - VMH ABTDG was added as the direct counterpart of the removed VMH ABTD1 (TASK.md).
- **Not done:**
  - no commit or push;
  - no merge with the R2041 (`ae54dd8`), GPR-fix (`claude/clever-kare-7d2e1a`) or R2190 work. All
    of them edit `main.py`/`cli.py` next to this flag; merge by hand;
  - no KO screen; no double KO; no builds without the E5 chain.
- **Colleague-owned audit.** The lipid-unlump `r1521_kegg_mnx_frontier_source_audit.json` still lists
  "R1964: D-/L-xylulose stereochemistry conflict". This candidate resolves the identifier side, but
  updating that file is the owner's decision.

## 6. Next gates

- The user decides whether this candidate becomes part of a reference chain, and how it is merged
  with the R2041 and GPR-fix candidates.
- R2190 plan: add the limit line from review §5.4. Option A remains independent of this change.
- Evidence that would upgrade or overturn the assignment: product identification of purified
  YlAraDH1/2 with D-arabitol + NAD⁺ (D-xylulose vs D-ribulose), the Liu 2026 *JAFC* full text, and
  Wang 2020's supplementary file.
