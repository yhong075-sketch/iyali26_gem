# V-ATPase R794/R795 rule (A1 + V0 e) and F-subunit scope (2026-11-03; merged as pipeline only 2026-11-04)

**Status:** the opt-in builder step `--vatpase-a1-gpr`, its curated data and test, and the Scope-1 benchmark edit are
merged as pipeline only (§8).

**Candidate:** the candidate built in this round on `2bc263f2` is **not kept**. On 2026-11-04 the user asked to merge
only the pipeline part, because main has since taken a better candidate (`e584151f`, PR #47).
- The candidate was `…_ceramide_pool_vatpase_a1_20261103.xml`, SHA-256
  `ab320e18477a186e77773de2acde761c17dc8319f989f5e66b467c88ef3e408e`.
- Its files are archived in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/vatpase_gpr_a1_20261103/candidate_ab320e18_not_kept/`.

The results in §1–§6 are for that candidate. §8 checks the step on main's line.

The scope, the user's authorization and the acceptance criteria are in [TASK.md](TASK.md).

## 1. What changed and why

**Before.** R794 (Golgi V-ATPase, [0, 1000]) and R795 (vacuolar V-ATPase, [0, 0]) both carried a rule inherited from
iYli21:

`(V1 A–H, V0 c, c′, d and a1 E12482) or (a2 F38820 and c″ F13017)`

The second branch has no ATP-hydrolysing head and no c/c′ ring, so it cannot pump. It kept both reactions open
under every single knockout. That made all 11 consensus-essential rule genes FN for a rule reason.

**After.** Both reactions get the same rule:

`A11258 and E14125 and D00581 and E10492 and A09766 and F20965 and YALI0E16192g and F31854 and F21690 and E37063 and
E32332 and F13017 and F25663 and (E12482 or F38820)` (all YALI1…g except the F subunit)

| Part | Change |
| --- | --- |
| V1 A–H and V0 c, c′, d | unchanged; YALI0E16192g is the V1 F subunit |
| c″ F13017 | joins the core |
| V0 e YALI1F25663g | new model gene, added to the AND |
| a1, a2 | the only OR in the rule |

**The new gene.** YALI1F25663g carries only cross-references verified locally in the NCBI feature table:
- NCBI Gene 90949792;
- RefSeq XP_065950526.1;
- SBO:0000243.

It has no name, UniProt or KEGG ID.

**Notes.** The reaction notes record:
- the decision;
- the standard used: it meets AGENTS.md §5 by homology and deliberately does not meet the 2026-09-21 same-compartment
  test;
- the evidence with its counter-evidence;
- the tier: a homology-based candidate rule, not a verified native rule.

**Unchanged:** stoichiometry, bounds (R795 stays closed), every other reaction and gene. That includes R_NTP1, which
still reads `YALI1A09766g or YALI1E41893g`.

**Scope-1.** YALI0E16192g is added to `model/expected/untargeted_genes_out_of_scope.json`, and the `benchmark.md`
sentence now names it. No library guide lies in its W29 locus.

**Decisions** (the user, 2026-11-03, on the decision write-up and the plan):

| Code | Decision |
| --- | --- |
| Std-§5 | AGENTS.md §5 decides the a-subunit OR |
| Fgate-lift | YALI0E16192g may be a required AND member |
| A1 | The same rule on R794 and R795 |
| Label-1 | Label the structural FN in this report and STATE only |
| Scope-1 | YALI0E16192g goes on the untargeted list |
| ID-keep | Keep the CLIB122 ID |
| E-1 | Add V0 e |
| Verify-1 | Build, diff, WT and 15 target knockouts |

Score effects were not a reason. The only effects are FN +1 and TN −1, both of which lower the score.

## 2. Evidence

Tags:
- **verified** = checked in this round or in the 2026-11-01/02 BLAST round by this session;
- **agent-verified** = a workflow agent did it, and it was not repeated;
- **recalled** = from an earlier report.

### Subunit identity

**13 of the 15 rule genes are reciprocal best hits of *S. cerevisiae* subunits** (agent-verified, BLAST 2026-11-01).
H E10492 is the weakest, at 26.8 %. The two exceptions:
- **a1 E12482:** its best *S. cerevisiae* hit is VPH1, 876 vs 811 bits over STV1, so its Stv1-type assignment is not
  reciprocal.
- **V1 F:** it has no W29 protein. CLIB122 Q6C5Q2's best hit is VMA7, one-way.

My first curation note said "all 15". The audit caught it, and the note was corrected (§5).

### The rule structure

**c″ belongs in the core.** The *S. cerevisiae* c″ null loses all activity
([Hirata 1997](https://europepmc.org/article/MED/9030535), abstract; verified 2026-11-01).

**The §5 OR test is met by homology.** Both *S. cerevisiae* a isoforms form active pumps:
- [Manolson 1994](https://europepmc.org/article/MED/7514599);
- [Kawasaki-Nishi 2001](https://europepmc.org/article/MED/11592965);
- [Perzov 2002](https://europepmc.org/article/MED/11948198).

These were verified from the downloaded copies on 2026-11-01. Which *Y. lipolytica* gene is the Stv1 type and which
the Vph1 type has only moderate support (tree support at most 67/100 across analyses; agent-verified).

**The 2026-09-21 same-compartment test is not met.** Substitution within one compartment is partial and one-way: the
Vph1 type covers for the Stv1 type in the early Golgi, but not the reverse
([Deschamps 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7002768/), full text, verified). A1 is recorded as meeting
§5 and not meeting that test, as the user chose (Std-§5).

### V0 e (YALI1F25663g)

**For adding it:**
- It is the reciprocal best hit of Vma9, at 61.4 % (agent-verified).
- The e null fails V1–V0 assembly ([Compton 2006](https://europepmc.org/article/MED/16569636), abstract, verified).

**Against:**
- Once the pump is assembled, e is not needed for coupled pumping in vitro
  ([Bueler 2015](https://europepmc.org/article/MED/25546637), abstract, verified).
- It is in no other *Yarrowia* model.

**Screens (E-2, verified 2026-11-03).** Read from the acCRISPR per-guide file MOESM6:

| Screen | Guides in the gene | Cutting scores (CS) | Depletion | Gene call |
| --- | --- | --- | --- | --- |
| Cas9 | 6 | 4 cut well (5.33–7.15) | strong; counts fall to about 0 by day 4–6 while no-Cas9 controls stay in the hundreds | FS −5.23, essential |
| Cas12a | 2, overlapping in exon 1 | −1.60 and −0.36: neither cut | none | FS +0.45, non-essential |

The Cas12a call is therefore a non-measurement (inferred). The 1.5 CS cutoff is recalled from the acCRISPR method, not
re-read. Evidence: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/vatpase_gpr_blast_20261101/e_guides/`.

### F subunit (YALI0E16192g)

- **The W29 "pseudogene, frameshift" call** comes from a missed 2-nt first exon plus a canonical 56-nt intron. Spliced,
  the locus gives the intact 122-aa CLIB122 protein. The splice was verified on 2026-11-01. That the gene is functional
  is inferred: no RNA evidence was checked. This is the same evidence class as COX9 (2026-10-30).
- **No library guide targets it:** 0 of 46,868 Cas9 and 57,014 Cas12a guides in NC_090774.1:1935676-1936100, and no
  guide ID naming E19360 or E16192 (verified 2026-11-02; the audit re-read the TSV).

## 3. Results

### Builds

All builds ran from `platform/` with `/opt/miniconda3/bin/python`, using 2bc263f2's argv
([scripts/args_2bc263f2.json](scripts/args_2bc263f2.json)) plus `--vatpase-a1-gpr`.

| Build | Seed | SHA-256 | Meaning |
| --- | --- | --- | --- |
| Control (flag off) | 0 | `2bc263f2…a062` | reproduces the baseline byte for byte (verified twice; the audit rebuilt it independently) |
| Candidate | 0, 1 | `ab320e18…408e` | seed-identical |

The exact commands are in [logs/build_commands.txt](logs/build_commands.txt).

Every build log, the control included, has a libSBML "Setting notes on R_R304" error. It is pre-existing: it is also
in the 2026-11-01 ceramide build. It was not investigated here.

### Whole-model diff

[model_diff_vs_2bc263f2.json](model_diff_vs_2bc263f2.json) has **5 differences**:
- R794 rule and notes;
- R795 rule and notes;
- gene YALI1F25663g added.

The other 14 genes keep their reaction links (A09766 also stays on R_NTP1). The audit's raw XML diff confirms that
nothing else changed. Mass and charge balance are unaffected, because no stoichiometry or species changed.

### Verify-1

16 LPs, from [verify1_ko_targets.json](verify1_ko_targets.json):
- Task A medium + PO1f;
- Gurobi with feasibility and optimality tolerance 1e-9;
- a reset before each solve;
- the targeted check.

| Point | Status | Growth | Reactions closed |
| --- | --- | --- | --- |
| WT | optimal, passed | 0.3120535436 (= baseline) | R794 flux 0 in this optimum |
| 13 core + e knockouts | optimal, passed | 0.3120535436 (KO/WT = 1 − 1.6e-14) | R794 (R795 already closed) |
| a1 E12482, a2 F38820 | optimal, passed | WT | none (R795 already closed) |

So closing R794 does not lower growth in this context: the 13 core knockouts close it and grow at the WT value. This
is the "pure ATP drain" argument, now observed in these 13 LPs. One optimum does not show that no other route
exists.

### Scores

[score_recount.json](score_recount.json) (`scripts/recount.py`):
- rule B (certified-infeasible knockouts score as growth 0);
- mtDNA and untargeted genes out of scope;
- negatives = model genes not on the consensus list (a user-defined class);
- a development reference, not independent validation.

| Cutoff | Baseline 2bc263f2 (old scope) | Candidate ab320e18 (new scope) |
| --- | --- | --- |
| 1 % | 114 / 212 / 37 / 701 | 114 / 213 / 37 / 700 |
| 5 % | 121 / 205 / 39 / 699 | 121 / 206 / 39 / 698 |
| 10 % | 125 / 201 / 44 / 694 | 125 / 202 / 44 / 693 |
| **15 %** | **142 / 184 / 44 / 694** | **142 / 185 / 44 / 693** |

Columns are TP / FN / FP / TN.

- **Changes:** F25663 enters as an FN, and YALI0E16192g leaves the TN count.
- **Positive coverage:** in-model positives 326 → 327, so out-of-model positives drop by one.
- **No existing gene changes class.**
- **Recall at 15 %:** 43.6 % → 43.4 %.

**How the candidate numbers were obtained.** The 15 rule genes use today's measured knockouts. The other 1,067 genes
(1,050 in scope) carry their class over from the baseline sweep. That carry-over is sound but not measured: the diff
shows the two LPs are the same optimisation problem for any knockout outside the rule. The audit confirmed this
reasoning and spot-checked 2 genes on both models with identical results. A full sweep (Verify-3) was not authorized.

### The 12 FN rule genes (Label-1)

The 12 FN rule genes are the 11 earlier ones plus V0 e:
- A09766, A11258, E10492, E12482, E14125, E32332, E37063, F13017, F20965, F21690, F31854, F25663 (all YALI1…g).

| Genes | Label | Cause |
| --- | --- | --- |
| The 11 genes other than E12482 | **Structural FN (growth in the model does not require this pump)** | The knockout closes only R794 (R795 is already closed), and no reaction or biomass term needs a Golgi or vacuolar proton gradient. Organelle acidification is not represented: there is no measured acidification demand to add, and adding an unmeasured one would be match-rate tuning. |
| E12482 (a1) | **Isoform redundancy** | The knockout closes nothing, because F38820 covers it in the OR. That is the same class as the other FN kept non-essential by an OR rule. |

These labels are backed by the 13 Verify-1 knockouts above. They stay counted as FN, and no benchmark list records
the labels.

## 4. Independent audit

A fresh read-only reviewer was given TASK.md, the diffs, the logs and the deliverables. It ran 3 builds and 11 LPs of
its own. Its verdicts on the 10 claims:

| Claim | Verdict |
| --- | --- |
| Builds | supported |
| Diff | supported |
| Boolean rule | supported |
| Apply step | supported, small gaps |
| Verify-1 | supported |
| Scores | supported, count off by one |
| Scope-1 | supported |
| E-2 guides | supported |
| Notes | partially supported |
| AGENTS.md | no contradiction |

**What was fixed after the audit:**
- **The notes overclaimed.** They said "all 15 are reciprocal best hits" (13 are), stated the F splice without
  "inferred", and left out the Bueler 2015 counter-evidence for e. All three are corrected in the curation data and
  TASK.md.
  - Because the notes are part of the model, this needed a rebuild: 9771e44a → ab320e18.
  - The diff between the two builds is the R794/R795 notes only.
  - Verify-1 was re-run on ab320e18 with identical results.
  - The first build and its logs are kept: model files in
    `$IYALI26_RESEARCH_ROOT/artifacts/tasks/vatpase_gpr_a1_20261103/superseded_9771e44a/`, logs in
    [logs/v1_9771e44a/](logs/v1_9771e44a/).
- **The count** "1,066 carried-over genes" was 1,067. It is fixed in TASK.md and `recount.py`.
- **Missing tests** for the applied-state checks (new gene with a changed annotation or name, rule gene with a new
  reaction) are added.

**Left as found:**
- The post-edit consistency check raises after editing an in-memory model. This is the same as the R171 step; the
  builder then aborts and writes nothing.
- The rule is duplicated as module constants that must equal the data file. This is the same convention as R171.
- The "rule applied but new gene missing" branch is not tested: cobra creates a rule's genes automatically, so the
  state is hard to construct.

## 5. Tests

| Run | Tree | Result |
| --- | --- | --- |
| `test_vatpase_a1_gpr.py` alone | this branch | 7 passed, 54 subtests |
| Full suite, before | origin/main `b57272d3` by `git archive` | 12 failed / 418 passed / 24 skipped / 56 errors ([logs/pytest_before_b57272d3.txt](logs/pytest_before_b57272d3.txt)) |
| Full suite, after | this branch | 14 failed / 465 passed / 24 skipped / 14 errors ([logs/pytest_after.txt](logs/pytest_after.txt)) |
| W29 audit test after regenerating the audit | this branch | 12 passed ([logs/pytest_after_w29_audit_rerun.txt](logs/pytest_after_w29_audit_rerun.txt)) |
| The 8 files that errored only before, re-run at `b57272d3` from a short path | `/tmp/claude-501/bt` | 42 passed ([logs/pytest_before_shortpath_8files.txt](logs/pytest_before_shortpath_8files.txt)) |

Both full runs used the `workflow-automation-review-6765da/platform/.venv` interpreter, which has no scipy. The
before tree lives under the session scratchpad, not in a git worktree.

**Test-by-test comparison** of the FAILED and ERROR lines:

- **Errors only in the before run (42).** These are setup errors in 8 files (cofactor biomass, heme iron, lipoylation on
  2c9cebd0, R1726, R2135, R2193, thiamine proton lock, TPP requirement). Each is cobra's `OSError` on a model path of
  260+ characters: the before tree lives under the longer scratchpad path. The same 8 files pass at `b57272d3` from a
  short path (42 passed). They are a path artefact, not a code difference.
- **Failures only in the after run (2).** Both are in `test_w29_tag_protein_audit.py`
  (`test_audit_covers_exactly_the_model_genes`, `test_raw_inputs_reproduce_the_committed_audit`). The new gene
  YALI1F25663g was not in the committed W29 tag-protein audit. I regenerated the audit as the test instructs
  (`python -m tools.audit_w29_tag_protein`). It gained one row, `YALI1F25663g,YALI1_F25663g,,mRNA,YALI1F:2566174-2566753-,
  XP_065950526.1,YALI1F:2566174-2566753-,ok`, plus the sidecar counts. The flag is `ok`, so no decision row is needed.
  The file then passes (12 passed). The full suite was not re-run after the regeneration.
- **Everything else is identical:** the 12 failures and 14 errors in both runs are the same tests, and the 7 new tests
  pass.

## 6. Limits and what was not tested

- **Only the 15 rule-gene knockouts were solved.** Verify-2 and Verify-3 were not chosen. The full sweep is inferred
  unchanged (§3).
- **The rule is a homology-based candidate:**
  - there are no *Y. lipolytica* localisation, knockout or complex data;
  - the a1/a2 isoform assignment has only moderate support;
  - the §5 OR test rests on *S. cerevisiae* data.
- **R795 stays closed** (UM2 Option A0), so its rule has no LP effect.
- **Still open, outside this round:**
  - the R_NTP1 rule (V1 A or an acid-phosphatase-like protein);
  - assembly factors VMA21 (B30163) and RAV1 (C23873);
  - coupling the pump to growth;
  - the 2 H+/ATP coefficient.
- **Disclosure (2026-11-01, BLAST round):** a workflow agent submitted two EBI Clustal Omega jobs (public sequences
  only) using the user's e-mail address as EBI's required contact field, without asking first. It is recorded in
  `vatpase_gpr_blast_20261101/aiso/runs.md`.

## 7. Files

**Changed:**
- `model/curation/gpr/vatpase_a1_gpr.json`
- `platform/scripts/gem_annotate/vatpase_a1_gpr.py`
- `cli.py` and `main.py` (flag `--vatpase-a1-gpr`, refused together with `--vatpase-gpr-hypothesis`)
- `platform/tests/test_vatpase_a1_gpr.py`
- `model/expected/untargeted_genes_out_of_scope.json`
- `model/expected/benchmark.md`
- this folder
- `model/STATE.md`

**Not merged (2026-11-04):**
- the candidate XML with its `.build.json` and `.coq9_genes.tsv` (archived);
- its `candidates/README.md` row;
- the W29 audit row for YALI1F25663g. The audit covers only genes in committed models, so it stays as on main. It
  must be regenerated when a candidate built with `--vatpase-a1-gpr` is kept.

**Scripts in this folder:**

| Script | Purpose |
| --- | --- |
| `before_state.py` | records the baseline state |
| `write_curation.py` | writes the curation data file |
| `run_builds.zsh` | runs the control and candidate builds |
| `model_diff.py` | whole-model diff (code unchanged from the ceramide round) |
| `ko_targets.py` | the Verify-1 LPs |
| `recount.py` | the score recount |
| `run_pytest.zsh` | the full pytest runs |

## 8. Pipeline-only merge on main (2026-11-04)

**Authorization.** The user, 2026-11-04: "we merge a better candidate than this, can you just merge the pipeline
part?". The user then chose:
- "Code + Scope-1 + docs" (no candidate XML, no README row, W29 audit as on main);
- "Commit and open a PR".

**Branch.** `claude/vatpase-gpr-a1-20261103`, rebased onto origin/main `9b6eafd8` (PR #47). The conflicts in
`cli.py`, `main.py` and `STATE.md` were additions on both sides, and both were kept. The step runs after the R24
step, before PR #47's mcnaugh_pass and BiGG-crosswalk steps.

**Checks** (verified 2026-11-04):

| Check | Result |
| --- | --- |
| In-memory apply on main's candidate `e584151f` | `applied` |
| Control build, `e584151f`'s argv with the flag off ([scripts/args_e584151f.json](scripts/args_e584151f.json), [scripts/run_builds_main.zsh](scripts/run_builds_main.zsh)) | reproduces `e584151f…39d2` byte for byte |
| The same argv plus `--vatpase-a1-gpr`, scratch only, not kept | `e8991a2c123fffbaa373a2ef4c4dd57eb84e62c8a58fee139494ee85b3f28233` |
| Its whole-model diff against `e584151f` ([logs/model_diff_main_vatpase_a1_vs_e584151f.json](logs/model_diff_main_vatpase_a1_vs_e584151f.json)) | the same 5 differences as on 2bc263f2: R794/R795 rule and notes, and gene YALI1F25663g |
| Focused pytest: this step, mcnaugh_pass, BiGG crosswalk, BiGG export, W29 audit, ceramide pool, 2026-09-11 V-ATPase hypothesis ([logs/pytest_main_merge_focused.txt](logs/pytest_main_merge_focused.txt)) | 52 passed |

**Not run on main's line:**
- no LP;
- no knockout sweep;
- no full pytest.

No score on `e584151f` + the flag is claimed. The Verify-1 results in §3 are for ab320e18 on 2bc263f2.

**Scope-1 on main's scoring.** It removes YALI0E16192g from the TN count of any later sweep. This was not recounted for
`e584151f`.
