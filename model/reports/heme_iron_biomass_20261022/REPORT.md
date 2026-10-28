# Heme a and iron(2+) in biomass_C, Yeast8 values (2026-10-22)

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). Exact commands are in
[commands.txt](commands.txt). Evidence tags: **verified** = checked this session with the output kept
here; **inferred** = reasoned, not run; **recalled** = from memory or an unopened source.

## 1. What changed and why

- **Trigger.** A Claude Science test ("Gap 5") showed that no biomass term contains iron. As a
  result, iron uptake and the whole heme route carried zero flux. Reproduced on `cf0825b8`
  (verified; TASK.md table).
- **Change.** `biomass_C` now also consumes:
  - heme a `m1635[C_cy]` at 1e-6 mmol/gDW;
  - iron `m948[C_cy]` at 3.04e-5 mmol/gDW.

  Seven `heme_iron_biomass_*` notes record the curation. Nothing else changed: the medium, the
  R1189 iron cap (1.61e-5), heme b and the iron chemistry are all as before.
- **Why these values.**
  - They are the *S. cerevisiae* consensus-model values (yeast-GEM v9.1.1).
  - Heme a at 1e-6 is also what the source model's own closed biomass R1710 carries.
  - The test's 5e-5 heme b + heme a has no source found, and it would make the cap bind
    (growth 0.161).
- **New candidate.**
  [E5_coq9_resp_pkg_r132_combined_gpi_coproducts_heme_iron_20261022.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_coproducts_heme_iron_20261022.xml),
  SHA-256 `7f836c3aabc0cc581f6ed2ed7e54d2d4e008ee5b2bdfaf24aab368781916539d`. It is `cf0825b8` plus
  `--heme-iron-biomass`.
  - Curated data: [heme_iron_biomass.json](../../curation/biomass/heme_iron_biomass.json).
  - Guarded step: `platform/scripts/gem_annotate/heme_iron_biomass.py`.
  - Test: `platform/tests/test_heme_iron_biomass.py`.

## 2. Evidence for the values

| Term | Value | Source | Tag |
| --- | --- | --- | --- |
| heme a [c] | 1e-6 mmol/gDW | yeast-GEM v9.1.1 `r_4598` cofactor pseudoreaction, `s_3714`; source table cites "Forster2003" | verified (model file and table read; archived) |
| iron(2+) [c] | 3.04e-5 mmol/gDW | yeast-GEM v9.1.1 `r_4599` ion pseudoreaction, `s_0924`; source table cites DOI 10.22161/ijeab/2.2.2 | verified (same) |
| heme b | none | no heme b in any yeast-GEM biomass pseudoreaction | verified |
| Primary measurements behind both values | — | Forster2003; DOI 10.22161/ijeab/2.2.2 | recalled (not opened) |
| Internal precedent | heme a `m1635[C_cy]` 1e-6 | source model `R1710` ("yeast 8 biomass pseudoreaction", closed 0/0) | verified ([before.json](before.json)) |

- **Links.**
  - [yeast-GEM.yml at 2d594ae1](https://github.com/SysBioChalmers/yeast-GEM/blob/2d594ae1c4a2d550ccef120d96a58c7bbf586255/model/yeast-GEM.yml)
  - [biomassComposition_Cofactor_Ion.tsv](https://github.com/SysBioChalmers/yeast-GEM/blob/2d594ae1c4a2d550ccef120d96a58c7bbf586255/data/physiology/biomassComposition_Cofactor_Ion.tsv)
  - [PR #174](https://github.com/SysBioChalmers/yeast-GEM/pull/174)
- **Archived copies.** `$IYALI26_RESEARCH_ROOT/artifacts/tasks/heme_iron_biomass_20261022/evidence/`
  (with `SHA256SUMS`). The raw files were fetched from `main`, and the GitHub API reported `main` at
  `2d594ae1` at that time.
- **Organism.** These are *S. cerevisiae* model values, not *Y. lipolytica* measurements.

## 3. Build and whole-model verification (verified)

- **Builds** (argv of `cf0825b8` from its `.build.json`):
  - without the new flag: `cf0825b8…`, identical bytes to the baseline;
  - with `--heme-iron-biomass`: `7f836c3a…` under both `PYTHONHASHSEED=0` and `1` (`cmp` identical).
- **Whole-model diff** ([model_diff_new_vs_cf0825b8.json](model_diff_new_vs_cf0825b8.json),
  `combined_candidate_20261020/scripts/model_diff.py`): 2 differences, both on `biomass_C`.
  - Stoichiometry: `m1635[C_cy]` −1e-6 and `m948[C_cy]` −3.04e-5 added.
  - Notes: 7 `heme_iron_biomass_*` keys added; the 5 existing notes unchanged.
- **Balance.** Only the pseudo-reaction `biomass_C` changed. Its mass balance is not computable by
  cobra (for example, chitin `m401` has no elements), so no balance claim is made. No metabolic
  reaction was touched.

## 4. Tests (verified)

- **Focused test with the final guard:** 4 passed, 22 subtests.
  - 16 stale-input cases fail closed with no edits, including:
    - a missing R1189;
    - a closed R863;
    - an altered R383;
    - an objective moved off `biomass_C`.
  - 6 scope cases are rejected.
  - The exact edit survives an SBML roundtrip and is idempotent.
  - The flag requires `--coq9-biomass-dilution`.
  - The first version had 17 subtests; the guard was extended after the audit (§8).
- **Full suite, test by test:**
  - [before](logs/pytest_before.log): 221 passed, 11 failed, 17 errors, 31 skipped.
  - [after, first guard](logs/pytest_after.log): the same outcomes plus the 4 new passes (225
    passed). This run overlapped with the builds and solves.
  - [after, final guard, run alone](logs/pytest_after_final_guard.log): 225 passed, 11 failed, 17
    errors, 31 skipped, 558 subtests. Identical test by test to the first after-run.
  - The 11 failures and 17 errors are pre-existing and unchanged. `pydantic_ai` is absent, so some
    gem_review files fail to collect.
- **Final-guard rebuilds:** seeds 0 and 1, into `.scratch`, both give `7f836c3a…`. The guard does
  not change the model. The final-code provenance is in
  [build_final_guard_seed0.build.json](build_final_guard_seed0.build.json). The candidate's own
  `.build.json` records the first guard version.

## 5. Solves (Task A medium + PO1f, GLPK; [solves.json](solves.json), 11 LPs; verified)

| | cf0825b8 | 7f836c3a |
| --- | --- | --- |
| WT growth (h⁻¹) | 0.3124267 | 0.3124263 (Δ −4.5e-7) |
| R1189 iron uptake | 0 | −9.81e-6 (61 % of the 1.61e-5 cap) |
| R2193 iron, ex → cy | 0 | 9.81e-6 |
| R863 iron, cy → mi | 0 | 3.12e-7 |
| R302 ferrochelatase, R384 heme o, R383 heme a, R1709 heme a export | 0 | 3.12e-7 each |
| R674 protoporphyrinogen oxidase | 0 | 1.56e-7 |

- With R1189, R2193, R302, R383 or R1709 closed in turn, the new model is optimal at zero growth
  (|growth| < 1e-17).
- **Mechanism (verified by the closures and the knockouts below):**
  - `biomass_C` now needs heme a in the cytosol, made in the mitochondrion by
    R61 → R667 → R409 → R792 → R790 → R240 → R674 → R302 → R384 → R383, then exported by R1709.
  - Ferrochelatase takes iron from R1189 → R2193 → R863.
  - Every step has one gene (or one AND rule), and there is no alternative route, so each step is
    required.
- About 97 % of the iron demand is the free-iron term. Heme a needs 3.1e-7 of the 9.8e-6.
- **Compartment convention.** Heme a is drawn from `m1635[C_cy]`, as in yeast-GEM ([c]) and the
  source model's R1710. This makes the gene-less export R1709 required. Biologically, heme a stays in
  mitochondrial cytochrome c oxidase (established, uncited), so `m609[C_mi]` was also an option.
  `biomass_C` already takes CoQ9 from `C_mi`. No gene result depends on this choice.

## 6. Single-gene knockouts and benchmark scores

The comparison ran on both models, 2 × (1,074 + 1) LPs, Task A medium + PO1f, GLPK with a 60 s
limit per LP. Outputs: [ko_compare.tsv](ko_compare.tsv) and
[ko_compare_summary.json](ko_compare_summary.json). Two fresh-solve corrections were then applied
(§7), giving [scores_corrected.json](scores_corrected.json).

### Genes that change class

There are 13, identical at 1, 5, 10 and 15 %. All go from non-essential to essential (KO growth 0).
No other gene changes class or growth by more than 1e-6. Identities are KEGG names for the
`kegg.genes` ID that the model annotates on each gene ([gene_identities.tsv](gene_identities.tsv),
read 2026-10-22). These are database identities, not experimental validation in *Y. lipolytica*. The
YALI1 → `yli:` mapping is the model's own annotation and was not re-checked.

| Gene | Reaction | KEGG name / ortholog | Experimental list |
| --- | --- | --- | --- |
| YALI1C07816g | R61 | 5-aminolevulinate synthase, K00643 | positive |
| YALI1F27686g | R667 | porphobilinogen synthase, K01698 | positive |
| YALI1F34051g | R409 | hydroxymethylbilane synthase, K01749 | positive (sweep `time_limit`; fresh solve 0) |
| YALI1F13039g | R792 | uroporphyrinogen-III synthase, K01719 | positive (sweep wrong; fresh solve 0) |
| YALI1C02391g | R790 | uroporphyrinogen decarboxylase, K01599 | positive |
| YALI1E35410g | R240 | coproporphyrinogen III oxidase, K00228 | **negative** |
| YALI1D30196g | R674 | protoporphyrinogen oxidase, K00231 | **negative** |
| YALI1F25924g | R302 | ferrochelatase, K01772 | positive |
| YALI1F31135g | R384 | heme o synthase, K02257 | positive |
| YALI1F32002g | R383 (AND) | heme a synthase COX15, K02259 | positive |
| YALI1B19490g | R383 (AND) | adrenodoxin-NADP+ reductase, K18914 | positive |
| YALI1B03314g | R383 (AND) | mitochondrial ferredoxin-2, K22071 | **negative** |
| YALI1F07747g | R2193 | "Phosphorus acquisition-controlling protein" (bHLH), no KO | **negative** |

### Scores

Negatives are the user-defined class: model genes not on the consensus list. The 16 mitochondrially
encoded genes are out of scope, and 1,290 listed positives have no model gene. The list was used in
model development, so these scores are a development reference, not independent validation.

| Cutoff | cf0825b8 TP / FN / FP / TN / unresolved, MCC | 7f836c3a (corrected) |
| --- | --- | --- |
| 1 % | 65 / 246 / 26 / 705 / 16, 0.281 | 74 / 237 / 30 / 701 / 16, 0.301 |
| 5 % | 71 / 240 / 28 / 703 / 16, 0.297 | 80 / 231 / 32 / 699 / 16, 0.315 |
| 10 % | 75 / 236 / 33 / 698 / 16, 0.294 | 84 / 227 / 37 / 694 / 16, 0.314 |
| 15 % | 83 / 228 / 33 / 698 / 16, 0.323 | 92 / 219 / 37 / 694 / 16, 0.341 |

At 15 %:
- recall goes from 0.267 to 0.296;
- specificity from 0.955 to 0.949;
- precision from 0.716 to 0.713.

Scores were not a reason for this change (AGENTS.md §2). A gained true positive does not validate
the mechanism.
- YALI1B19490g (adrenodoxin reductase) and YALI1B03314g (ferredoxin) are essential in the model
  only through heme a.
- In *S. cerevisiae*, the essentiality of Arh1/Yah1 is usually attributed to Fe–S cluster
  biogenesis (recalled, uncited). The model does not contain Fe–S cluster biogenesis.

**What the four new false positives mean** (inferred, not tested):
- **YALI1F07747g.** The model's only iron-import gene, but its annotated ID is a transcription factor
  in KEGG. Real cells have several iron-uptake systems (recalled, uncited). So this FP most likely
  reflects the model's single-gene iron import, not biology.
- **YALI1E35410g (HEM13), YALI1D30196g (HEM14), YALI1B03314g (ferredoxin).** Their single-step GPRs
  make them lethal in the model. Possible reasons they are not on the list:
  - screen coverage;
  - other enzymes or routes;
  - ID mapping.

  None of these was checked.

## 7. A solver artifact found in the sweep

This change introduces the artifact. It is the main limit of this result (also in §9).

- **First run (observed).** The sequential GLPK sweep hung for 34 min inside one simplex solve,
  confirmed by a stack sample showing GLPK factorisation routines. The sample is archived at
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/heme_iron_biomass_20261022/logs/ko_compare_run1_stack_sample.txt`.
  - The run was stopped ([log](logs/ko_compare_run1_stopped.log)), and the script gained a 60 s
    limit per LP and a progress log.
  - The first run printed no per-gene lines, so which LP hung is unknown. It was probably
    YALI1F34051g, the LP that later timed out (inferred).
- **YALI1F34051g (R409).** The rerun hit `time_limit` on the new model. A fresh GLPK solve is
  optimal at 0 in 0.2 s, and Gurobi agrees: 0 on the new model, 0.31243 on the baseline
  ([ko_fresh_resolve.json](ko_fresh_resolve.json)).
- **YALI1F13039g (R792, HEM4).** The sweep reported the full WT growth.
  - **Verified:** a fresh GLPK pFBA gives optimal 0 ([hem4_bypass.json](hem4_bypass.json)).
  - **Verified:** a warm-started GLPK solve (WT first, then the knockout in the same context)
    reproduces the WT value. R792 is at 3.12e-7 while its bounds are (0, 0), and the status is still
    "optimal". The audit found this to be the only bound violation above 1e-12 in that solution.
  - **Verified:** fresh Gurobi gives 0 at tolerance 1e-7 and 1e-9
    ([hem4_tolerance_check.json](hem4_tolerance_check.json)). That script's docstring wrongly gives
    Gurobi's default as 1e-6; the recorded value is 1e-7. Its residual metric covers |S v| only,
    not bounds.
  - **Inferred:** GLPK kept the WT basis because the violation is within its scaled primal tolerance.
- **Why this change exposes it (inferred).** Heme a at 1e-6 × growth 0.31 puts every heme-route flux
  at 1.6–3.1e-7. That is at the edge of GLPK's primal tolerance (default 1e-7, recalled).
- **Risk to other genes.**
  - Accepting a stale basis gives growth at or above the true optimum. That can be WT growth or an
    intermediate value (inferred).
  - The other 11 heme/iron genes read 0 in the sweep, which is the true optimum.
  - The audit listed every reaction with 1e-9 < |v| < 1e-5 in the WT and pFBA optima of both models.
    The only ones with genes are the 11 heme/iron reactions and YALI1D22124g's acyl-CoA ligases
    (R1445, R1451, R1452). That gene is essential in both sweeps. This is a screen, not a proof: each
    sweep solve starts from the previous knockout's basis, not from WT.
  - The baseline sweep and earlier GLPK sweeps of other candidates are exposed in the same way and
    were not re-examined.
- **How to rule it out in future screens** (proposed, untested):
  - after every knockout solve, check that all variables respect their bounds, and in particular
    that knocked-out reactions are ≤ 1e-9; treat any violator as unresolved and re-solve it fresh
    (the Gurobi screen script already computes `violation`);
  - or solve each knockout from a fresh basis;
  - or use Gurobi with FeasibilityTol ≤ 1e-9.

  To test, rerun every gene called non-essential, in both models, with one of these methods and
  compare.
- **Other tools (inferred).** Any FVA or blocked-reaction tool that treats |v| < 1e-6 as zero will
  misread this route.

## 8. Independent audit

The full audit is in [AUDIT.md](AUDIT.md). A fresh read-only reviewer used 18 LPs.
- **Verdicts:**
  - C1–C4 and C6–C9, C11: supported.
  - C5 (guard): partially supported.
  - C10: supported for HEM4/HEM3; whether other genes are affected is only partly settled.

Responses:

| Issue | Response |
| --- | --- |
| 1. Numerical fragility | Accepted. It is now the first limit in §9 and is stated in §7. Changing the coefficient or the solver protocol is a user decision. |
| 2. Budget statement | Accepted. Corrected in §9 and in TASK.md: the stopped first run went beyond the sweep budget, and the rerun was started without a separate check-in. |
| 3. Guard gaps | Fixed in part. The guard now checks that `biomass_C` is the only objective, and pins R1189, R863, R302, R384 and R383 besides R1709 and R2193. The focused test gained 5 cases (22 subtests). Final-guard rebuilds still give `7f836c3a` with seeds 0 and 1. **Left open:** values of existing `biomass_C` notes (their HTML escaping differs between the builder's in-memory model and a model read from file, so only keys are pinned); R1710's bounds; the heme steps upstream of R302; enforcing `baseline_model_sha256`. |
| 4. Evidence not kept | Fixed. Saved: [scripts/reproduce_gap5.py](scripts/reproduce_gap5.py) and [scripts/plan_a_in_memory.py](scripts/plan_a_in_memory.py), with their transcribed output in [logs/reproduce_gap5_output.txt](logs/reproduce_gap5_output.txt); the stack sample and the 13 raw KEGG entries in the research workspace; the `.outcomes` command in commands.txt. |
| 5. Concurrency | Recorded in commands.txt and §4. The final-guard suite ran alone. |
| 6. `hem4_tolerance_check.py` | Noted in §7. The script was left as run. The bound violation is shown by the recorded R792 flux and the audit's check. |
| 7. Compartment convention | Stated in §5 and TASK.md. |
| 8. Score gain is not validation | Stated in §6. |
| 9. TASK.md wording | Corrected in a section added to TASK.md; the original text is kept. |
| 10. §7 overstated | Rewritten, with verified and inferred statements labelled. |
| 11. Open items | §8 and the STATE.md entry are written. The authorization quote is from the user's message in this session. |

**New finding from the audit:** the UniProt cross-reference on YALI1F07747g, A0A1D8NM31, is deleted
in UniProt, and KEGG's entry for 2908631 links to Q6C2U3. This is added to the R2193 finding.

## 9. Limits, untested items and next gates

- **Numerical fragility (main limit).**
  - Heme a at 1e-6 puts the heme-route fluxes at 1.6–3.1e-7, at LP tolerance.
  - In this evaluation, 2 of the 13 heme/iron knockouts came out wrong or unresolved in the GLPK
    sweep and needed fresh solves (§7).
  - Any later sweep, FVA or blocked-reaction analysis of this candidate needs one of the checks in
    §7.
- **Condition.** Only the Task A medium + PO1f with GLPK. SD-Leu, the overlay conditions and other
  media were not run.
- **Not comparable.** These GLPK scores are not directly comparable with the Gurobi screens in
  earlier entries.
- **Organism.** The values are yeast-GEM (*S. cerevisiae*) consensus values; no *Y. lipolytica*
  measurement was found.
- **Iron headroom.** The total iron demand of 3.14e-5 mmol/gDW uses 61 % of the cap. Any later iron
  term (heme b, Fe–S clusters) that brings total cellular iron above about 5.15e-5 mmol/gDW would
  make the cap bind.
- **Findings recorded, not changed** (TASK.md):
  - **The YALI1F07747g GPR on R2193.** KEGG names its annotated ID a bHLH regulator. Its UniProt
    cross-reference A0A1D8NM31 is deleted, and KEGG links that gene to Q6C2U3 (audit). Followed up
    in [r2193_iron_import_gpr_20261022](../r2193_iron_import_gpr_20261022/REPORT.md): the gene is a
    PHO4-like bHLH TF, and the iron-uptake genes are two FET3/FTR1 pairs.
  - **Iron chemistry.** The iron species have no formula and charge 0, and R302, R727, R384 and R383
    are imbalanced.
  - **Siroheme** has no consumer.
  - **Heme b** is absent.
- **Budget (not fully kept).**
  - Builds: 3 authorized; 5 run. The 2 extra are the scratch-only final-guard rebuilds after the
    audit.
  - Diagnostic LPs: 17 (11 + 2 + 4), within the ~20 budget.
  - Fresh re-solves: 5.
  - Knockout sweep: 2,150 LPs as authorized, **plus** the stopped first run. Its LP count was not
    recorded (up to about 2,100, inferred), so it went beyond the authorized sweep.
  - The rerun with a 60 s limit per LP was started without a separate check-in.
  - The planning turn used 10 LPs before authorization.
- **Next gates (user decisions):**
  - adopt this candidate into the working line or not;
  - the R2193 gene;
  - a solver protocol for near-tolerance fluxes;
  - heme b;
  - the iron chemistry.
- Nothing is committed or pushed.

### What kind of statement each result is

- **Statically verified:**
  - the diff;
  - byte identity;
  - the guard behaviour (tests);
  - the coefficients matching yeast-GEM.
- **Observed in runs:**
  - WT growth;
  - fluxes;
  - closures;
  - the knockout classes;
  - the HEM4 warm-start artifact.
- **Provable from structure:** with the new term, any blocked heme-route step forces zero growth,
  because heme a has no other source (inferred from the producer lists; consistent with every solve).
- **Untested:**
  - other media;
  - the solver-protocol fixes in §7;
  - the biology behind the four false positives.
