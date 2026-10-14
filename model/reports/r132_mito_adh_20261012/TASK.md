# TASK B: R132 "mitochondrial alcohol dehydrogenase" (2026-10-12/13)

## Authorization

The user's prompt file `claude_code_prompt_w29_tasks_20261012.md`, sections 0–2 plus Task B
("read and do the task B"), with these decisions given in this conversation on 2026-10-13:

- **Predictors:** submit the 9 public RefSeq sequences to the official TargetP-2.0, DeepLoc (DTU)
  and MitoFates (AIST) web forms through the built-in browser. No account, no email.
- **Task A medium:** copy the uncommitted Task A file byte-identical into this worktree.
- **Medium for steps 4–5:** the user answered "PO1f + CSM-Leu-". Read as: Task A's medium
  unchanged (PO1f runtime overlay plus CSM-Leu supplements, glucose 2.414). Glycerol and oleate
  runs use the same file with only the carbon source swapped.
- **Code plan:** approved, run through.

**Scope limits** (prompt §0):
- No push, merge or branch switch.
- No edit of any existing model XML, medium, evidence or benchmark file.
- No deletions of files.
- No change made to improve a match rate.
- One curation change, built to a new file in `model/candidates/`. Everything stays uncommitted.

## Starting state (verified 2026-10-13)

| Item | Value |
|---|---|
| Worktree / branch | `claude/task-b-e1ff88` @ `7dc9511` (merge of PR #3). `git diff 2a34475 HEAD` is empty. |
| Baseline model | `model/candidates/E5_coq9_respiratory_package_20261008.xml`, SHA-256 `0ac2addfa0cb6de2d89393e433766fb04325ffe42c2d131eaefaca3bcec29dba` |
| Medium (screen and W29 runs) | `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, SHA-256 `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45`, copied from worktree `competent-turing-bc6e9b` (Task A, uncommitted there) |
| PO1f profile | `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json`, `35307853…` |
| Experimental list | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv`, `1e887f5a…` |
| W29 data | `$IYALI26_RESEARCH_ROOT/reference/w29_worland2024/`, as recorded by Task A (`growth_parameters.csv` `0264fa47…`, `glucose_central_fluxes.csv` `844615d0…`). Link to Worland et al. 2024: **inferred, pending user confirmation**. |
| Research root | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` (from `~/.zshrc`; exported per command) |
| Main checkout | `restructure/model-platform` @ `2a34475`, untracked `.codex/` and `.mcp.json` (not ours; untouched) |

**Target R132 as built in the baseline:**
- `m278[C_mi] + m30[C_mi] <=> m279[C_mi] + m27[C_mi]`, i.e. acetaldehyde + NADH ⇌ ethanol + NAD⁺
  in the matrix. Bounds (−1000, 1000), no notes.
- GPR: `YALI1D02041g or YALI1C08017g or YALI1E15452g or YALI1A15227g or YALI1F32330g or
  YALI1F11481g or YALI1B11402g or YALI1E23809g or YALI1F37194g`.
- Matrix ethanol `m279[C_mi]` is used only by R132 and R1114. Deleting R132 would leave R1114 as a
  dead end.
- 7 of the 9 genes are used only by R132. The two exceptions are YALI1A15227g (also R130, R1377,
  R1864, R1866, R1868, R1885–R1887, R2100, R_ALCDH_nadp_hi) and YALI1E15452g (also R257, R513,
  R801).
- The builder keeps genes without reactions; the baseline already has 15. So the benchmark gene
  count stays 1,073 either way.

Full baseline state of every target: [before.json](before.json).

## Decision rule (fixed 2026-10-13T19:10Z, before any predictor was run)

**Predictors** (official web servers, full RefSeq sequences, no trimming):
- TargetP-2.0, organism group Non-Plant;
- MitoFates, organism group Fungi;
- DeepLoc 2.x, whichever version DTU serves; recorded, but not part of the rule.

**A gene counts as predicted mitochondrial** if TargetP-2.0 mTP ≥ 0.5 **or** MitoFates
probability ≥ 0.385. This is the complement of the prompt's "none" condition (mTP < 0.5 and
MitoFates < 0.385 for every gene).

**Branch 1. No gene is predicted mitochondrial.**
- Propose removing R132.
- Evidence of absence:
  - no Adh3-type ortholog with a presequence;
  - Kerscher et al. 2001;
  - Gatter et al. 2016;
  - each citation opened before it is written into evidence.
- R1114 and R997 have no GPR (diffusion) and stay.
- Genes and metabolites are kept.

**Branch 2. Some genes are predicted mitochondrial.**
- Propose a GPR restricted to those genes, with per-gene citations, and flag the capacity as
  unknown.
- If no gene-level citation supports any predicted gene as an alcohol dehydrogenase, stop and
  report instead of building.

The rule is applied exactly. A score that sits on a threshold is reported as it is, without
rounding.

## Other evidence, gathered alongside

- **Local Smith-Waterman** (BLOSUM62, gap open 11 / extend 1). The 12 Yarrowia proteins (9 R132
  genes plus the Gatter ADH1/2/3 genes in R130) are aligned against:
  - *S. cerevisiae* Adh1–5, Sfa1, Adh6, Adh7 and Bdh1;
  - *K. lactis* ADH3 and ADH4, the mitochondrial ADHs of a Crabtree-negative yeast, for the
    prompt's yeast-orthology rule.

  Sequences and outputs: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r132_mito_adh_20261012/`.
- **Citations to open before use:**
  - Kerscher et al. 2001 J Cell Sci 114:3915–3921;
  - Gatter et al. 2016 Fungal Genet Biol 95:30–38;
  - Bakker et al. 2000 J Bacteriol;
  - the UniProt Adh3 transit-peptide annotation.

## Implementation (Branch 1 or 2, one curation change)

- **Curated data:** `model/curation/r132_mito_adh/r132_mito_adh.json`, holding R132's exact before
  state, the predictor scores and the evidence.
- **Apply step:** `platform/scripts/gem_annotate/r132_mito_adh.py`, `apply_r132_mito_adh()`.
  - Fails closed unless R132, R1114, R997, m279 usage and the gene RefSeqs match exactly.
  - A rerun is a no-op.
- **Wiring:** `cli.py` / `main.py` flag `--r132-mito-adh`. It requires
  `--coq9-respiratory-package` and runs right after it.
- **Test:** `platform/tests/test_r132_mito_adh.py`, on an in-memory model. Then the full `pytest`
  run, compared test by test with a clean run.
- **Builds:**
  - Control: rebuild the baseline into `.scratch/` and require SHA `0ac2addf…`.
  - Candidate: `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`.
- **Whole-model diff:** only R132 may differ, plus, under Branch 1, the genes whose reaction lists
  shrink as a direct result.

## Computation (declared budget)

Solver: Gurobi, 1 thread, Presolve 0, tolerances 1e-7. The same contexts are used for baseline and
candidate.

1. **Controls** on the Task A medium plus PO1f profile, both models:
   - WT;
   - KO/WT for R2062 (complex I), R305 (complex III), R304 (complex IV) and R171 (ATP synthase).
2. **Full gene benchmark,** both models: 2 × 1,074 LPs. Scoring follows `benchmark.md`:
   - 15 % primary, plus 1 %, 5 % and 10 %;
   - mtDNA genes out of scope;
   - negatives are the user-defined class;
   - a knockout with no solution is unresolved, with the "counted as lethal" view as a labelled
     second column, as in Task A.
3. **Route analysis** (WT and complex I KO, before and after): FVA on the shuttle and the
   replacement route, plus one pFBA vector per panel for the maps. FVA, not single pFBA vectors,
   is used for judgments.
4. **Step 5, diagnostic only.** The W29 13C upper bounds go in a new ASSUMPTION file,
   `model/conditions/assumptions/w29_13c_upper_caps_20261012.json`:
   - PC R690 ≤ 0.058 × uptake;
   - net cytosol → matrix malate ≤ 0.061 × uptake, summed over every carrier that moves malate
     between the two compartments;
   - glycine cleavage ≤ 0.002 × uptake. The earlier 0.009 is also reported: it was the
     synthesis-side upper bound of the MFA row "R51 net".

   Runs: glucose, glycerol and oleate at their measured uptake, each as supplied and with R132
   off plus the caps. Glycerol and oleate use their own MFA columns for the caps.

   New media files: the Task A medium with glucose closed, plus glycerol 6.2666 or oleate 0.7630.
5. **Budget:** at most 3,000 LPs and 30 min wall time, builds excluded.

## Stop conditions

- The predictors cannot be run.
- The whole-model diff shows anything besides R132 (and, under Branch 1, the direct gene-list
  consequence).
- The baseline rebuild does not reproduce `0ac2addf…`.
- WT is infeasible.
- A recalled value differs by more than 2 %: stop, report, do not tune. Exception: the WT 0.297
  "W29-type defined medium" value, which the user replaced with the Task A medium. The difference
  is reported, not treated as a stop.

## Addendum 2026-10-13: predictor results, Branch 2 selected, stop fired

**Runs.** All three ran on the official servers with the full RefSeq sequences from
[r132_genes.fasta](r132_genes.fasta) (`3c437523…`).

| Tool | Settings | Job ID | Result date |
|---|---|---|---|
| TargetP-2.0 | Non-Plant, long output | `6ACE82AE0001C6E02ECCAE2B` | 2026-10-13 21:12 UTC (server timestamp) |
| MitoFates | fungi | `20261014.040043.60854` | server-local time, JST |
| DeepLoc 2.0 | High-quality (Slow), ProtT5, long output | `6ACE830B0001C7579B8DCFFC` | 2026-10-13 |

The raw outputs are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r132_mito_adh_20261012/predictors/`.
[collect_predictions.py](collect_predictions.py) writes [localization_predictions.tsv](localization_predictions.tsv)
and [decision_rule_result.json](decision_rule_result.json).

**Rule outcome.** Only YALI1F11481g counts as predicted mitochondrial, so the rule selects
**Branch 2**. The other 8 genes all have TargetP mTP ≤ 0.0104, MitoFates ≤ 0.006 and DeepLoc
mitochondrion ≤ 0.4434, with DeepLoc never predicting them mitochondrial.

| YALI1F11481g (XP_505149.2, 401 aa) | Value | Threshold |
|---|---|---|
| TargetP-2.0 mTP | **0.502811** (CS 27–28, Pr 0.2265) | 0.5, crossed by 0.0028 |
| MitoFates | **0.197**, "No mitochondrial presequence" | 0.385, not crossed |
| DeepLoc 2.0 Mitochondrion | **0.6242**, mitochondrial transit peptide signal | 0.6220 (DeepLoc's own), crossed by 0.0022; not part of the rule |

**Branch 2 sub-clause: no gene-level ADH citation exists for YALI1F11481g → stop before any code.**
Every check below was run this session.

- **Annotation:** RefSeq XP_505149.1 and KEGG yli:2908523 call the protein "NADPH-dependent
  alkenal/one oxidoreductase". UniProt Q6C2G3 lists Pfam ADH_N (PF08240) and ADH_zinc_N_2 (PF13602),
  i.e. the MDR superfamily, with no function comment. The W29 UniProt entry A0A1D8NMH1 named in the
  model annotation is now **inactive**.
- **Literature:** Europe PMC finds 3 papers for YALI0F08129g / YALI1F11481g. Two open-access full
  texts were read:
  - PMC6201493 lists YALI0F08129g among 11 "putative alcohol dehydrogenase genes" found by BLAST.
    Its deletion gave no 1-butanol phenotype. This is not activity evidence.
  - PMC13640256 calls YALI1_F11481g "a mitochondrial inner membrane protease", with no experimental
    support given; this conflicts with the MDR annotation.
  - PMC5092062 (J Vis Exp 2016) is not open access and was not read.
- **Gene-model dependence:** XP_505149.1 (340 aa, CLIB122 model) is exactly residues 62–401 of
  XP_505149.2 (401 aa, CLIB89 model). TargetP's cleavage site (27–28) and MitoFates' amphipathic
  helix (18–27) both lie inside the 61-aa extension, which only the newer gene model has. Which
  start codon is real was not checked.
- **Homology:** YALI1F11481g is 22–26 % identical to Sc Adh3/Adh5 and KlADH3/KlADH4 over local
  alignments that start at its residue 119 or 121. The four Adh1-like Yarrowia proteins (A15227g,
  and the R130 genes D33961g, E21128g, A16292g) are 68–71 % identical to the mitochondrial Sc Adh3,
  KlADH3 and KlADH4. They align from mitochondrial-ADH residue 26–30, i.e. none carries the
  presequence. Table: `artifacts/.../sw_results.tsv`.

**Not done** because the stop fired: no curated file, code, test, build, solve, map, audit,
REPORT.md or STATE.md entry. The research-root artifacts and this folder's files are kept.
Decision requested from the user.

## Addendum 2026-10-13: start-codon check (user: "Check start codon first")

**Authorization.** The user chose to check the start codon before deciding: submit the 340-aa
XP_505149.1 once to each predictor (same settings), and look for transcript support of the
61-aa extension. Only those three submissions were made. Transcript support was looked up in
records only; no BLAST or other submission.

**Predictors on both gene models** ([start_codon_check.tsv](start_codon_check.tsv); raw outputs
under `artifacts/.../predictors/*_XP_505149.1/`):

| Protein | Length | TargetP mTP | MitoFates | DeepLoc 2.0 mito (threshold 0.6220) |
|---|---|---|---|---|
| XP_505149.2 (CLIB89 model) | 401 aa | 0.502811 | 0.197 | 0.6242, called Mitochondrion |
| XP_505149.1 (CLIB122 model) | 340 aa | 0.000138 | 0.001 | 0.3923, called Cytoplasm \| ER |

The mitochondrial call exists only with the 61-aa extension.

**Genomic DNA** (NCBI NC_090775.1 CLIB89 chr 1F and NC_006072.1 CLIB122 chr F, fetched this session):
- The 1,023-nt CDS, the 183 nt that encode the extension and 600 nt further upstream are
  identical in both strains (0 mismatches). The two RefSeq versions differ only in which start
  was annotated, not in sequence.
- The 401-aa ATG is the most upstream start possible: the nearest in-frame stop is 51 nt
  upstream.
- Start context:
  - 401-aa ATG: `TGCGCTTCC ATG ACT`, with T at −3;
  - 340-aa ATG: `CAATACACC ATG AAA`, with A at −3.

  That a purine at −3 favours start use is general eukaryotic knowledge, **established / not
  cited**, and it was not tested here for *Y. lipolytica*.

**Transcript support:** none specific to this gene was found.
- XM_505149.3 is a PROVISIONAL RefSeq "derived from an annotated genomic sequence (NC_090775)",
  "incomplete on both ends", with no RNA-seq evidence block.
- Magnan et al. 2016 PLoS ONE 11:e0162363 (PMC5014426, opened this session) built YALI1 models
  by merging three pipelines and selected CDSs "based on size and RNA-Seq support". This is a
  genome-wide rule that favours the longest ORF, not gene-level evidence for this start. Their
  RNA-seq data "will be published elsewhere".
- The W29 UniProt entry A0A1D8NMH1 (401 aa) was deleted as "Not part of a reference proteome",
  an administrative reason. The CLIB122 entry Q6C2G3 (340 aa) has only a rule-based (ARBA)
  "Lipid droplet" location.

**File note.** `artifacts/.../start_codon/clib122_locus_minus_strand_783up_plus_cds.fasta` came
from a failed window search (CDS not found, wrong slice). It is **invalid** and was kept under
the no-delete rule; the valid file is `clib122_locus_gene_strand_783up_cds.fasta`. Earlier, two
of this task's own just-written outputs were removed from the same artifacts folder: a 150-byte
Europe PMC error response and the misplaced CLIB122 window FASTA.

Decision requested from the user.

## Addendum 2026-10-13: user decision, remove R132

The user chose "Remove R132" after seeing all predictor scores and the start-codon check. This is
a **user amendment made after the scores**, not the outcome of the fixed rule, which had selected
Branch 2.

**Recorded reason:**
- The only mitochondrial call (YALI1F11481g) exists only with an upstream start codon that has no
  gene-specific transcript support. With the CLIB122 start, all three predictors give non-mitochondrial
  (TargetP 0.0001, MitoFates 0.001, DeepLoc 0.392).
- The protein has no alcohol dehydrogenase evidence.
- No R132 gene with ADH homology carries a presequence.

**Implementation:** as planned for Branch 1.
- R132 is removed. R1114 and R997 stay; all genes and metabolites are kept.
- Evidence is cited only after each source has been opened.

## Addendum 2026-10-13: build done, runs done, recalled-value stop fired

**Done.** All of this is verified this session.

- **Curated data:** [r132_mito_adh.json](../../curation/r132_mito_adh/r132_mito_adh.json).
- **Apply step:** `apply_r132_mito_adh` in `platform/scripts/gem_annotate/r132_mito_adh.py`, behind
  the `--r132-mito-adh` flag.
- **Test:** `platform/tests/test_r132_mito_adh.py`, 4 tests passed.
- **Full pytest:** identical to a clean HEAD copy test by test, apart from the 4 new tests. The
  11 failures and 14 errors are pre-existing.
- **Baseline rebuild:** reproduced `0ac2addf…` exactly.
- **Candidate:** `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`, SHA-256
  `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`.
- **Whole-model diff:** only R132 is removed. The other differences are its own pathway-group
  membership, the 9 genes' reaction lists and the usage lists of its 4 metabolites.
- **Runs:** [run.json](run.json) and [net_malate_fva.json](net_malate_fva.json), 2,911 LPs in total,
  within the 3,000 budget. The second script's docstring says 8 LPs; it used 12, including 4 growth
  solves.

**Stop fired: recalled values differ by more than 2 % in the user's chosen medium.** Details:

- The CI KO shuttle requirement is ≥ 257.4 per 100 glucose against a recalled ~270.
- After removal, net malate import is −2.45 in the CI KO and [86.5, 142.4] in WT, against a recalled
  ~48.
- Step-5 complex I is 186.2 against ~177.
- Step-5 predicted/measured growth is 0.979 / 1.107 / 1.091 against 0.93 / 1.07 / 0.99.
- "As supplied" over-prediction is +17 / +33 / +23 % against +12 / +29 / +18 %.

On an amino-acid-free glucose medium (labelled reconciliation), the recalled glucose values are
reproduced: +11.9 % as supplied, 0.927 with the R361 cap. **Inferred:** the differences come from
the medium substitution the user chose, not from the model change. This was not tested for
glycerol, oleate or the route analysis.

**New finding** (verified, step-5 runs):
- R361 re-oxidizes the dihydrolipoamide of 2-oxoglutarate dehydrogenase. The R361 cap "as specified"
  forces R614 (2-OGDH) to 0.2 per 100 glucose, against an MFA estimate of 16 (95 % CI 0–29).
- None of the three glycine-cleavage entry steps carries flux in any variant.
- So the cap's effect (1.069 → 0.979 of measured growth on glucose) comes from throttling 2-OGDH, not
  from glycine cleavage.

Nothing was tuned. Decision requested from the user before the audit, REPORT.md and STATE.md entry.

## Addendum 2026-10-13: user decision, reconcile in a W29-type defined medium (written before these runs)

**Authorization.** The user chose Option 2 and raised the budget to about 3,300 LPs in total, i.e.
at most 389 more LPs after the 2,911 already used.

**W29-type defined medium,** as defined by the user, applied in memory:
- Start from Task A's medium: `po1f_csm_leu_w29uptake_20261012.csv`, or the glycerol / oleate files
  derived from it, with the PO1f profile loaded.
- Close adenine and the 11 amino acids: R1003, R1202, R1204, R1215, R1217, R1220, R1222, R1223,
  R1231, R1232, R1233, R1234.
- Restore R612 (URA3) to (0, 1000).
- Fix substrate uptake at the measured value (lower = upper bound): glucose R1070 2.414, glycerol
  R1141 6.2666, oleate R1440 0.7630.
- Two variants:
  - **V2 (primary):** uracil R1354 closed.
  - **V1:** uracil open at 1000.
- The profile's R45 LEU2 GPR swap has no effect here, because no gene is knocked out in these runs.

**Variants:**
- *Holistic overlay:* the baseline plus the prompt-header diagnostic overlay, reused verbatim from
  Task A's `run_task_a.py` (`f28cef16…`, uncommitted in worktree `competent-turing-bc6e9b`).
- *Holistic + shuttles restricted:* the candidate (R132 removed) plus the overlay plus the 13C caps
  (PC, net malate import, and glycine cleavage in three forms).
- Context runs: "as supplied" (baseline, no overlay) and "candidate + overlay, no caps".

**Expected predicted/measured growth** (user-provided; glucose / glycerol / oleate):

| | V2 | V1 (reproduces the recalled values) |
|---|---|---|
| holistic overlay | 1.103 / 1.279 / 1.138 | 1.116 / 1.293 / 1.151 |
| holistic + shuttles restricted | 0.916 / 1.054 / 0.983 | 0.928 / 1.068 / 0.994 |

**Stop rule, V1 only.**
- It fires if any V1 ratio differs by more than 2 % (relative) from its expected value.
- The "shuttles restricted" comparison uses the R361 cap at the **looser** fraction, which is the
  setting of the earlier run being reproduced.
- The R361 0.002 form and the glycine-cleavage-module form are reported next to it, without a stop.
- The V2 and PO1f + CSM-Leu comparisons are informational.

**Also run:**
- the step-4 route FVA (baseline vs candidate, WT and complex I KO, glucose) in V2 and V1, plus the
  net-malate-import range;
- the PO1f + CSM-Leu medium as a third column, labelled **sensitivity (adds amino acids)**, with the
  same overlay variants. Its uptake is "≤ measured", as in that medium file.

**Planned LPs** (about 350 of the 389 allowed):
- step-5 growth: 51;
- V2 glucose FVA on R570, R2062, R614 and R690, plus pFBA for the maps: about 66;
- route V2: 4 × 30;
- route V1: 4 × 28.

**After this:** if the stop does not fire, run the audit, write REPORT.md and add the STATE.md entry.

## Addendum 2026-10-13: W29 runs done; corrections after the independent audit

**W29-type defined runs.** [w29_defined_run.json](w29_defined_run.json), 337 LPs.
- The V1 stop did not fire: every V1 value is within 0.03 % of the user-expected values.
- V2 also matches its expected values.
- Total 3,248 LPs (2,899 + 12 + 337), within the raised budget.

**Corrections.** The earlier addenda are kept as written; these entries supersede them.

1. **TargetP timestamp.** "2026-10-13 21:12 UTC (server timestamp)" is server local time, CEST
   (UTC+2), i.e. **19:12:47 UTC**. The second TargetP job was 21:25:51 CEST = 19:25:51 UTC.
2. **MitoFates job IDs.** The IDs embed timestamps that read, as JST, 19:00:43 and 19:13:14 UTC.
   - The submitted FASTA was created at 19:11:52 UTC and the XP_505149.1 FASTA at 19:24:44 UTC (file
     birth times).
   - The results were saved at 19:14:52 and 19:28:18 UTC.
   - So the MitoFates server clock reads about 12–14 min behind, consistently for both jobs. No job
     ran before the rule.
   - Ordering evidence (filesystem birth times, UTC): TASK.md, with the decision rule, at 19:11:24;
     before.json 19:11:43; r132_genes.fasta 19:11:52; TargetP output saved 19:14:34; DeepLoc
     19:16:48.
   - TASK.md is untracked and was appended later, so this ordering rests on file metadata, not on
     version control.
3. **"TargetP mTP ≤ 0.0104"** for the other 8 genes is precisely ≤ 0.010418.
4. **Deletions: three files were removed, not two.** The earlier addenda undercount them. None was
   authorized by the user, so the removals conflict with this task's "no deletions" scope. All three
   had been written by this task minutes earlier in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r132_mito_adh_20261012/`:
   - `evidence/epmc_30359264_fulltext.xml`: a 150-byte Europe PMC error response, from querying the
     full text with a wrong PMC ID;
   - `start_codon/clib122_NC_006072.1_1180425-1182230_minus.fasta`: an NCBI efetch window at the KEGG
     coordinates, which did not contain the gene;
   - `start_codon/clib122_NC_006072.1_1150000-1210000.fasta`: a 60-kb efetch window that did not
     contain the gene.

   Both windows can be re-fetched with the same efetch queries (db=nuccore, id=NC_006072.1, the
   stated coordinates, strand=2 for the first). The valid locus file is
   `clib122_locus_gene_strand_783up_cds.fasta`.
5. **"R361 looser" is looser only for glucose.** The UB of the MFA row "R51 net" is the
   glycine-synthesis side, not a cleavage bound. It is *tighter* than the cleavage cap for glycerol
   (0.000949 vs 0.001422) and oleate (0.015284 vs 0.108917). The variant name is kept because the
   recorded outputs use it.
6. **Solves outside the LP tally.** The 3,248 LPs count the task's run scripts only. Not counted:
   - the 5 optimizations of an exploratory medium probe (scratch `probe_medium.py`, run before the
     budget was declared, to find which medium gives WT ≈ 0.297);
   - the optimization tests inside the two full `pytest` runs (repository test suite, clean HEAD and
     worktree).
7. **Glycine-cleavage "zero" claim is narrower than written.** R357, R360 and R362 = 0 is verified
   by FVA only in the 15 PO1f + CSM-Leu step-5 variants of run.json, which have no overlay.
   - In the 4 V2-glucose variants it rests on pFBA vectors. Each is an optimal solution of the growth
     LP, so a cap that is slack there does not bind: for an LP, removing a constraint that is
     inactive at an optimum leaves the optimum unchanged.
   - It was not computed for V1, V2 glycerol/oleate, or the CSM-Leu overlay variants.
8. **w29_defined_run.json has no inputs block.** The script asserts the SHA-256 of both models, all
   media, the profile, the experimental list, the mtDNA list, the W29 flux file and the overlay
   source. The caps file has no recorded digest there; its SHA (`df1e33a8…`) matches run.json and was
   not modified after either run.
