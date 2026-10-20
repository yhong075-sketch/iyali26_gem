# R765/R766 without YALI1D02625g: opt-in change set and candidate (2026-10-19)

- **Scope.** Authorization, baseline, acceptance criteria and the amended control gate:
  [TASK.md](TASK.md). Evidence for the decision:
  [../d02625g_identity_20261019/REPORT.md](../d02625g_identity_20261019/REPORT.md).
- **Not committed.** This change is uncommitted. Only the identity report (`dd406b4`) and the merge
  of `claude/clever-kare-7d2e1a` (`c091b78`) were committed, both at the user's instruction. Nothing
  was pushed.
- **No LP was run on any GEM** (baseline, control or candidate). The two authorized full pytest runs
  include existing toy-fixture solves.

Tiers:
- **verified**: computed or read in this session and shown in a linked file;
- **inferred**: reasoned from verified facts.

## 1. What changed and why

- **Science.** R765 "transketolase 1" and R766 "transketolase 2" go from
  `YALI1E07744g or YALI1D02625g` to `YALI1E07744g`.
  - YALI1D02625g is a transketolase N-terminal (PP-domain) fragment: no PYR or C-terminal domain, and
    no conserved GDG Mg²⁺/ThDP motif.
  - An OR member needs evidence that it catalyses alone (AGENTS.md §5). There is none. The tier of
    "it cannot catalyse alone" is **inferred** (identity report §8).
  - It entered only through the 2026-06-02 EC-match isozyme table. The source model has
    `YALI1E07744g` only.
  - YALI1D02625g is a real, expressed gene, so it stays in the model, unlinked.
- **Recall.** YALI1E07744g is an essentiality FN in the 2026-09-10 screen: label essential, KO/WT
  1.0000, FN at 10%
  ([essentiality_per_gene.tsv](../screen_test_metadata_trna_20260910/essentiality_per_gene.tsv)).
  No LP was run, so any change of that call is untested here. It is a consequence, not the reason.
- **Mechanism in the model** (verified, [static_checks.json](static_checks.json), knockout bounds only,
  no solve):

  | Knockout | 51c6fabd closes | Candidate closes |
  | --- | --- | --- |
  | YALI1E07744g | none | R765, R766 |
  | YALI1D02625g | none | none (the gene has no reactions) |
  | both | R765, R766 | R765, R766 |

## 2. Implementation (through the pipeline)

| File | Change |
| --- | --- |
| `model/curation/gpr/selection_tkl_fragment_gpr_removal.json` (new, sha256 `db656c60…92f5`) | See below |
| `platform/scripts/gem_annotate/selection_gpr_restore.py` | `TKL_FRAGMENT_GPR_RESTORE` constant, one `CHANGE_SETS` entry (targets R765, R766) and a two-line docstring bullet. The apply logic is unchanged |
| `platform/scripts/gem_annotate/cli.py` | Flag `--remove-selection-tkl-fragment-gpr` (`store_true`, default off) |
| `platform/scripts/gem_annotate/main.py` | Parameter, plus a call to `restore_after_selection(model, selection, TKL_FRAGMENT_GPR_RESTORE)` after the three existing restores. The result is stored as `selection_tkl_fragment_gpr_restore` in the build record |
| `platform/tests/test_selection_gpr_tkl_fragment.py` (new) | 7 in-memory tests (below) |
| `model/candidates/README.md` | One row for the new candidate, added after the builds |

The data file records:
- change ID `SELECTION-GPR-TKL-D02625G-20261019`;
- for R765/R766: cleared = `YALI1E07744g or YALI1D02625g`, restored = `YALI1E07744g`, with the
  equation and bounds unchanged;
- `unchanged_equation_residual` `{charge: 2.0}`;
- a gene guard on YALI1E07744g (name `COBRAProtein646`, sbo `SBO:0000243`);
- `removed_genes` `["YALI1D02625g"]`;
- the evidence text, tier "inferred" and sources.

**Tests**, as written:
- record scope;
- apply + idempotence (YALI1D02625g kept with no reactions; no gene or metabolite added);
- knockout propagation;
- an SBML round trip that is a no-op;
- 8 stale-input cases that fail closed without edits:
  - GPR `YALI1D02625g`;
  - GPR `YALI1E07744g and YALI1D02625g`;
  - an equation edit;
  - a bounds edit;
  - a species charge change;
  - a gene name change;
  - a conflicting note;
  - a missing reaction;
- the flag is off by default (argparse level) and requires a complete selection;
- on the real 51c6fabd model, applying the set changes only R765/R766 in `model_definition`, and gene
  names and annotations stay the same.

**Limits of the tests** (independent audit; not fixed, because editing tests is a code change):
- The 51c6fabd test does not assert the baseline SHA. It does not compare notes or annotations of
  other reactions, metabolite annotations, or other genes' reaction links. The whole-model diff in §4
  covers these on the built candidate.
- There is no stale case for the sbo/annotation guard, and none for a missing gene.
- The `main.py` wiring is evidenced only by the two build records, not by a unit test.
- Plan step 4's "the three existing change sets still reproduce exactly" was not written as a new
  test. It rests on the unchanged existing tests (§4) and the field-level control (§3.1).

**Guard strength** (audit finding, verified):
- `SBO:0000243` is carried by all 1,073 genes, and `COBRAProtein646` is an auto-generated name.
- So the guard catches an ID or name remap of YALI1E07744g. It does not check biological identity.
- The removed gene YALI1D02625g (`ncbigene` 2910520 in the model) is not guarded: `removed_genes`
  is a record field, and the apply code never reads it.

**Deviations from the approved plan:**
- Two were approved by the user before the code was written (TASK.md, authorization 3):
  - the gene guard uses name and sbo, because YALI1E07744g has no `ncbigene` annotation in the model.
    NCBI Gene 2912397 and Q6C6T4 appear in the evidence text and sources only;
  - the pre-existing +2 charge residual is recorded as `unchanged_equation_residual`. The equation
    is not edited.
- The third is the missing reproduction test above.

**Residual.**
- The 2026-09-10 metadata selection set these equations
  (`metadata_reaction_selection.json`, verified):
  - its "before" equations have no H⁺. They are charge-balanced but short of 2 H in mass;
  - its "after" equations add 2 H⁺ (R765) and 4 H⁺ (R766), which makes them mass-balanced with a +2
    charge residual.
- **Arithmetically (verified), the residual traces to two species:**
  - m889 R5P, stored as C5H11O8P at charge −2. At charge 0, R765's charge sum would be 0;
  - m166 F6P, stored as C6H11O9P at charge 0. At charge −2, R766's sum would be 0.
- That these two formula/charge pairs follow different protonation conventions (neutral-acid
  formula at −2, dianion formula at 0) is **inferred** and uncited. It is not a verdict that either
  species is wrong. A formula or charge edit would propagate to every reaction that uses it, and
  needs its own review. R671's recorded residual is +1.

## 3. Builds (offline, `--no-solve`; exact commands in [commands.txt](commands.txt))

| Build | Output | SHA-256 |
| --- | --- | --- |
| Control: baseline flags, no new flag | `.scratch/tkl_d02625g_gpr_20261019/control_no_new_flag.xml` (git-ignored). A copy and its build record are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/tkl_d02625g_gpr_20261019/control_build/` | `248ea2c0d8d69aa579a5a70a1a1860fb9434d46d927d2e2b5afdc797c0ae9995` |
| Candidate: baseline flags + `--remove-selection-tkl-fragment-gpr` | `model/candidates/E5_coq9_resp_pkg_r132_tkl_fragment_gpr_20261019.xml` | `52d3106bc6b598ade840ef56a12938aec9590efd5dc39419400175a59ca87f34` |
| Candidate build record | `…tkl_fragment_gpr_20261019.build.json` | `c1684e8e19304436d3d82ace852b1ef5a9bc53902eabd35eaf7e16688cd677af` |

- **Code identity.** Both builds ran at the same time, 04:29:59–04:32:30 UTC, on HEAD `c091b78`
  plus the uncommitted change in §2 ([logs/build.head](logs/build.head),
  [logs/build.status](logs/build.status)).
  - The audit confirmed that all 43 `source_sha256` entries in the build records equal the current
    tree.
  - For every file in 51c6fabd's build record, the code at `c091b78` is the code that built it.
- **Logs:** [logs/build_control_no_new_flag.log](logs/build_control_no_new_flag.log) and
  [logs/build_candidate.log](logs/build_candidate.log).
  - Each has one `LibSBML error code -3` line on the R_R304 notes. That line appears in every earlier
    build, as recorded in the 2026-10-19 fix report.

### 3.1 Control gate: not byte-identical; accepted on field-level identity by the user

- **Byte level (verified, [logs/control_byte_diff.txt](logs/control_byte_diff.txt)):** the control is
  `248ea2c0…`, not `51c6fabd…`. The whole-file diff is 4 lines, all in R159's `listOfReactants`. The
  species and coefficients are the same; only their order differs.

  | Model | R159 reactant order |
  | --- | --- |
  | 51c6fabd | m267, m199, m130, m32, m10, m141 |
  | control | m267, m141, m32, m199, m130, m10 |
  | candidate | m267, m32, m199, m141, m10, m130 |

  The product order (m35, m326, m143, m50) is the same in all three.
- **Field level (verified, [model_diff_control.json](model_diff_control.json)):** 0 differences over
  every reaction, metabolite, gene, group, compartment, objective and model-level field. The auditor
  re-derived this independently, including object order.
- **Cause:**
  - R159 is edited by the existing change set SELECTION-CPA-URA2-RESTORE-20261019, not by this
    change. Its restore adds 7 species to R159 (m10, m130, m141, m143, m199, m32 and m50). It also
    removes m325 and changes m35's coefficient.
  - `apply_selection_gpr_restore` builds its delta by iterating
    `set(current) | set(target["stoichiometry"])` (`selection_gpr_restore.py` lines 157–158), then
    calls `add_metabolites`. Cobra appends new species in iteration order, and the SBML writer keeps
    that order.
  - Python randomizes string hashing per process. **Inferred** from the code.
- **Static support (verified, re-run by me after the audit).**
  - The auditor's `seedsim.py` replays lines 157–158 on R159's pre-restore order (from a56f3649).
  - Under `PYTHONHASHSEED` 399, 317 and 188, it gives exactly the three orders in the table above.
  - It must run without `python -I`, which ignores `PYTHONHASHSEED`. My first re-run used `-I` and
    did not match.
  - The script and the commands are in commands.txt and in the artifacts folder. This is a
    simulation of the expression, not a rebuild.
- **Testable prediction** (needs two authorized builds, not run):
  - the `$BASE` build with `PYTHONHASHSEED=399` should give `51c6fabd…` byte for byte;
  - the candidate build with `PYTHONHASHSEED=188` should give `52d3106b…`.
- **Decision.** The user, 2026-10-19: "Accept the control gate on field-level identity
  (0 differences; byte difference only R159 speciesReference order from set iteration)."
  - Consequence: an unseeded rebuild of 51c6fabd or of this candidate reproduces its content but
    usually not its bytes. A seeded rebuild probably reproduces both (prediction above).
  - A sorted-iteration fix would change the byte layout of every build with
    `--restore-selection-cpa-ura2`. It is proposed as a separate task, not done here.
- **Sequence.**
  - The mismatch was known at 21:32:40 PDT.
  - Within the authorized budget and before asking the user, I ran the field-level diffs, the static
    checks and the after-change pytest run (21:34:08–21:39:05).
  - The candidate had already been built at the same time as the control.
  - The acceptance time was not recorded. Details are in the TASK.md amendment.
- **Controls not rebuilt:** a56f3649 (no restore flags) and c98961e9 (only
  `--restore-selection-cleared-gprs`). Both were outside the two authorized builds, and neither
  uses the CPA/URA2 flag.

## 4. Whole-model verification (verified)

- **Candidate vs 51c6fabd, field by field** ([model_diff.json](model_diff.json), `model_diff.py`
  from the 2026-10-19 fix report, unchanged). There are exactly 5 differences:
  - R765 `gpr`: `YALI1E07744g or YALI1D02625g` → `YALI1E07744g`;
  - R765 `notes`: the new keys `selection_gpr_restore`, `…_evidence`, `…_evidence_status` and
    `…_sources`. All existing note keys are unchanged;
  - R766 `gpr` and `notes`: the same;
  - gene YALI1D02625g `reactions`: [R765, R766] → [].
- **Notes:** the note values equal the data file after `html.unescape`. The reloaded text holds a
  literal `&apos;`, a known SBML reload quirk that the module handles.
- **Byte diff vs 51c6fabd** ([logs/candidate_byte_diff.txt](logs/candidate_byte_diff.txt)): the
  R765/R766 notes and `fbc:or` blocks, plus the R159 order noise of §3.1.
- **Static checks** ([static_checks.json](static_checks.json), the final inline script reproduced in
  commands.txt):
  - R765/R766: bounds [−1000, 1000], residual `{charge: 2.0}` in both models;
  - knockouts as in §1;
  - 1,073 genes and 2,306 reactions in both models;
  - build record:
    - the selection is complete;
    - the R765/R766 selection rows are "applied" (stoichiometry applied, gpr already correct, as in
      51c6fabd);
    - all four change sets are "applied".
- **pytest** ([logs/pytest_comparison.txt](logs/pytest_comparison.txt)):
  - before (at `c091b78`): 181 passed, 31 skipped, 11 failed, 14 errors, 431 subtests passed;
  - after: 188 passed, 443 subtests passed, otherwise identical;
  - every pre-existing test has the same outcome; the 7 new tests pass, including all 12 of their
    subtests.
  - The 11 failures and 14 errors are pre-existing:
    - failures: `test_coa_protonation_curation` ×6, `test_vacuole_supply` ×2,
      `test_vacuole_explicit_inputs`, `test_dfba_new_fn_essentiality`, `test_lipid_moiety_ledger`;
    - errors: `test_lp_sn12_candidate` ×14.
  - None of them mentions R765, R766, D02625, E07744, `selection_gpr_restore` or `CHANGE_SETS`
    (audit grep). The audit traced the failures to a missing `../iyali26_gem/model.xml`.

## 5. Not done, limits, next gates

- **No LP on any GEM.** The 6 regression LPs need a separate authorization: WT, KO YALI1E07744g and
  KO YALI1D02625g, on 51c6fabd and on the candidate. So does the benchmark at 1/5/10/15%.
  - Expected without solving: WT is unchanged, because R765/R766 stay catalysed by YALI1E07744g, and
    the medium and strain overlay are identical (inferred).
  - Whether the YALI1E07744g knockout falls below each cutoff is unknown.
- **Byte reproducibility** (§3.1): a seeded rebuild could confirm the cause (two builds). A
  sorted-iteration fix, followed by a documented re-baseline, is proposed as a separate task.
- **Test strengthening** (§2 limits): needs the user's go-ahead, because it is a code edit. Proposed:
  - assert 51c6fabd's SHA in the baseline test, and compare the full `model_diff.py` table;
  - add stale cases for the sbo annotation ("identity differs") and for a missing YALI1E07744g;
  - assert that YALI1D02625g carries `ncbigene` 2910520 before the change.
- **Not changed:**
  - Lines 22–23 of `gpr_isozyme_additions.csv` (data rows 21–22) still add YALI1D02625g, and the
    metadata selection still re-adds it.
  - This change set removes it after both, and only when its flag is set.
  - Making that the default (alternative (b) in the identity report) is a separate decision.
- **Evidence limits** (identity report): no *Y. lipolytica* enzymology for either gene; no nr
  ortholog result; no whole-genome tblastn.
- **Records:** the `model/STATE.md` entry is added with this report. Publishing needs a separate
  instruction.

## 6. Independent audit

- **Auditor.** A fresh read-only agent was given only TASK.md, the diff, the logs and the
  deliverables. It ran no build, no pytest and no LP; reading SBML initialized Gurobi, but nothing
  was solved.
- **Verdict on the model.** The candidate's content is correct and confined to the targets. It
  re-derived this independently. All SHAs, counts and pytest outcomes reproduce. The code change is
  confined to what §2 says.
- **Claims.** Of 31 claims checked, 25 were supported (one only by file times), 5 partially
  supported and 1 overstated. No defect touches the candidate's content.

| # | Finding | Disposition |
| --- | --- | --- |
| 1 | The "before" knockout claim cited a file without baseline knockouts | Fixed: baseline knockouts added to static_checks.json (§1) |
| 2 | Cause testable statically; "cannot be reproduced byte for byte" overstated | Fixed: seed simulation re-run and cited (§3.1); wording now "unseeded rebuild"; seeded rebuild listed as a prediction |
| 3 | Evidence generation not reproducible from the record | Fixed: commands.txt now has exact paths and the inline scripts verbatim. Copies are in the artifacts folder (`scripts/`). No new script file in the repository |
| 4 | Control evidence only in git-ignored `.scratch/` | Fixed: control XML, build record and genes TSV copied to the artifacts folder, with SHA256SUMS |
| 5 | Tests weaker than §2 implied | Disclosed (§2 limits); strengthening proposed for approval (§5). Not edited |
| 6 | Gene guard strength not disclosed; removed gene unguarded | Disclosed (§2); test assertion proposed (§5) |
| 7 | "No LP was run" imprecise | Fixed: "no LP on any GEM; pytest includes toy-fixture solves" |
| 8 | Stop condition and gate amendment not documented in sequence | Fixed: TASK.md amendment and §3.1 sequence |
| 9 | Imprecise R159 order description | Fixed: full orders tabulated (§3.1) |
| 10 | Residual attribution imprecise and untiered | Fixed: traced to m889/m166 arithmetically (verified); convention mismatch tagged inferred; no verdict (§2) |
| 11 | Minor omissions | Fixed: FN source linked, LibSBML −3 line, subtest counts, third plan deviation, STATE entry |

The fixes were not re-audited.
