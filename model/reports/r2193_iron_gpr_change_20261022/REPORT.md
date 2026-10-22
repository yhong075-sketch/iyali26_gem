# R2193 iron import: GPR to the two FET3/FTR1 pairs (2026-10-22)

Scope, authorization, budget and acceptance criteria are in [TASK.md](TASK.md). It includes the stops
during the run, the user's decision on the knockout check, and the corrections after the audit. Exact
commands are in [commands.txt](commands.txt). The evidence is in the investigation report
[r2193_iron_import_gpr_20261022](../r2193_iron_import_gpr_20261022/REPORT.md); see its correction
note.

Evidence tags:
- **verified** = checked this session with the output kept here;
- **inferred** = reasoned;
- **recalled** = not opened.

## 1. What changed and why

- **Before.** R2193 (iron(2+) import, `m1273[C_ex]` → `m948[C_cy]`) had the GPR `YALI1F07747g`. That
  gene encodes a PHO4-like bHLH transcription factor, not a transporter. Since the heme/iron biomass
  term (`7f836c3a`), iron import is required, so this regulator gene was lethal in the model.
- **After.** R2193's GPR is `(YALI1D09277g and YALI1D09308g) or (YALI1D08684g and YALI1D08564g)`.
  - These are the *Y. lipolytica* FET3/FTR1 pairs, D09 and D08.
  - YALI1F07747g stays in the model with no reaction.
  - The four new genes carry their NCBI Gene, RefSeq, UniProt and KEGG cross-references.
  - R2193's name, bounds, stoichiometry and annotation are unchanged. R863 is unchanged.
- **Candidate (current).**
  [E5_coq9_resp_pkg_r132_combined_gpi_coproducts_heme_iron_r2193_gpr_v3_20261022.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_coproducts_heme_iron_r2193_gpr_v3_20261022.xml),
  SHA-256 `e1586027cd5e3727a75b237097d3578d2eba976b9162294f0d6ae2b8be1a127b`. It is `7f836c3a` plus
  `--r2193-iron-gpr`.
  - Data: [r2193_iron_gpr.json](../../curation/gpr/r2193_iron_gpr.json).
  - Step: `platform/scripts/gem_annotate/r2193_iron_gpr.py`.
  - Test: `platform/tests/test_r2193_iron_gpr.py`.
- **Superseded builds, deleted at the user's request** (2026-10-22, with their `.build.json` and
  `.coq9_genes.tsv`; the diffs below are kept):
  - `…_r2193_gpr_v2_20261022.xml` (`2cf19c69`): the same model with the uncorrected evidence notes;
  - `…_r2193_gpr_20261022.xml` (`4d738403`): the v2 model with a hash-seed-dependent gene order.

## 2. Evidence (investigation, corrected after the audit)

| Gene | Role | Identity evidence | Essential list |
| --- | --- | --- | --- |
| YALI1D09277g (Q6C9Y8) | FET3-like ferroxidase, D09 pair | reciprocal best hit of FET3 (553 bits); Pfam multicopper oxidase; KEGG K19791; detected 24/24 | positive |
| YALI1D09308g (Q6C9Y7) | FTR1-like permease, D09 pair | reciprocal best hit of FTR1 (392 bits); PF03239; K07243; detected 15/24 | positive |
| YALI1D08684g (Q6CA12) | FET3-like ferroxidase, D08 pair | paralog; one-way best hit of FET3 (536 bits; FET3's best hit is Q6C9Y8); K19791; detected 23/24 | not listed |
| YALI1D08564g (Q6CA15) | FTR1-like permease, D08 pair | paralog; one-way best hit of FTR1 (346 bits; FTR1's best hit is Q6C9Y7); K07243; detected 24/24; co-expressed with VHb, raising intracellular iron (Ye et al. 2026, primers matched; no Ftr1-only control) | not listed |
| YALI1F07747g (Q6C2U3) | PHO4-like bHLH transcription factor | HLH domains only; PHO4 reciprocal best hit over the bHLH region; Pomraning 2016 | not listed |

- **BLAST inputs.** BLAST used the CLIB122 UniProt proteins. The W29 D08684g protein (644 aa, against
  589 aa for Q6CA12) was not searched.
- **AND within a pair** follows the *S. cerevisiae* Fet3–Ftr1 complex (Stearman et al. 1996). That
  the *Y. lipolytica* pairs work the same way is inferred.
- **OR between pairs.** It rests on both pairs being full FET3/FTR1 homologues with Pfam, KEGG and
  protein-detection support. It is weaker for D08 than for D09, because the D08 genes are paralogs,
  not reciprocal best hits.

## 3. Build and whole-model verification (verified)

- **First builds** (`4d738403` / `b3d66568`): bytes differed by seed. This was a stop condition. Only
  the order of the four new `geneProduct` entries differed, because cobra creates genes new to a GPR
  in set order.
- **v2.** The step registered the genes in sorted order: `2cf19c69` under both seeds.
- **v3** (after the audit):
  - corrected evidence notes;
  - the genes are created by growing the rule one gene at a time, through cobra's own context-aware
    path;
  - notes restored on context exit.
- **v3 builds:**
  - control without the flag: `7f836c3a` (also for the first and v2 controls);
  - with the flag: `e1586027` under both seeds (`cmp` identical).
- **Whole-model diffs:**
  - [v3 against `7f836c3a`](model_diff_v3_vs_7f836c3a.json): 7 differences, namely R2193's GPR and
    notes, the four added genes, and YALI1F07747g's reaction link;
  - [v3 against v2](model_diff_v3_vs_v2.json): only two R2193 note texts;
  - [v2 against the first build](model_diff_v2_vs_first_build.json): 0 fields.
- **Balance.** Stoichiometry is unchanged. R2193 stays uncheckable, because the iron species have no
  formula.
- **Provenance.** `PYTHONHASHSEED` is recorded in commands.txt, not in build.json.

## 4. Tests (verified)

- **Focused test (v3 code):** 6 passed, 27 subtests, under seeds 0, 1 and 7.
  - 15 stale-input cases fail closed with no edits.
  - 5 scope cases are rejected.
  - The exact edit and sorted gene order survive an SBML roundtrip; the step is idempotent.
  - 7 knockout combinations follow the pair logic.
  - Exiting a `with model:` block restores the baseline (new).
  - The flag is opt-in.
- **Full suite after v2** ([pytest_after.log](logs/pytest_after.log), run alone): 230 passed, 11
  failed, 17 errors, 31 skipped. Identical test by test to the previous run
  ([heme final-guard](../heme_iron_biomass_20261022/logs/pytest_after_final_guard.log)) plus the 5
  new passes.
- **Full suite after v3** ([pytest_after_v3.log](logs/pytest_after_v3.log), run alone): 231 passed,
  11 failed, 17 errors, 31 skipped. Identical test by test to the v2 run, plus the new context-exit
  test. The failures and errors are pre-existing.

## 5. Knockouts and scores (verified)

The sweep ran on `7f836c3a` and v2: 1,074 and 1,078 genes, 2,154 LPs, Task A medium + PO1f, GLPK, 60 s
per LP.
- **Method:** a basis reset before every knockout, then the targeted check (user decision): a point
  fails if a reaction the knockout closed carries |v| > 1e-9, or any bound or |S v| violation
  exceeds 1e-6.
- **Outputs:** [ko_compare.tsv](ko_compare.tsv) and [ko_compare_summary.json](ko_compare_summary.json).
- **v3 inherits these results.** v3 differs from v2 only in note text, which does not enter the LP.
  v3's WT is 0.3124262956, the same as v2.
- **Experimental list:** `consensus_essential_genes.csv`, SHA-256
  `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b`, 1,612 positives.
  - 322 → 324 are model genes. 1,290 → 1,288 have no model gene (out of scope).
  - 11 in-model positives are unresolved (infeasible knockouts) in both models.
  - 16 mtDNA genes are out of scope.
  - Negatives are the user-defined class.
  - Development reference, not independent validation.
- **The check found nothing on either model.**
  - Every closed reaction carried exactly 0.
  - The largest general violation was 4.0e-7 and the largest |S v| 2.5e-7. No point failed, and
    there were no fresh re-solves.
  - The same 20 knockouts are infeasible in both models. They are unresolved without a re-solve (a
    narrowed rule, §8).
- **WT:** 0.3124262956 on both models.
- **Cross-check.** The baseline sweep reproduces the heme report's corrected scores at all four
  cutoffs without fresh re-solves.
  [method_check_targeted.json](method_check_targeted.json) shows the reset fixing the HEM4 artifact.
  Under the current check, the warm-start point is flagged (R792, closed, at 3.12e-7) and the reset
  point passes at 0.

**Changes** (identical at 1, 5, 10 and 15 %):

| Gene | 7f836c3a | v2 / v3 | Experimental list |
| --- | --- | --- | --- |
| YALI1F07747g | essential (KO growth 0) | non-essential (KO = WT) | negative: FP → TN |
| YALI1D09277g, YALI1D09308g | not in model | non-essential | positive: new FN ×2 |
| YALI1D08684g, YALI1D08564g | not in model | non-essential | negative: new TN ×2 |

No other gene changes growth by more than 1e-6.

**Scores** (TP / FN / FP / TN / unresolved, MCC):

| Cutoff | 7f836c3a | v2 / v3 |
| --- | --- | --- |
| 1 % | 74 / 237 / 30 / 701 / 16, 0.301 | 74 / 239 / 29 / 704 / 16, 0.303 |
| 5 % | 80 / 231 / 32 / 699 / 16, 0.315 | 80 / 233 / 31 / 702 / 16, 0.317 |
| 10 % | 84 / 227 / 37 / 694 / 16, 0.314 | 84 / 229 / 36 / 697 / 16, 0.315 |
| 15 % | 92 / 219 / 37 / 694 / 16, 0.341 | 92 / 221 / 36 / 697 / 16, 0.342 |

- **The two new false negatives are the stated limit of the OR rule.** Experimentally the D08
  paralogs do not compensate for the D09 pair, and a GPR cannot express that. Choosing the D09 pair
  alone would have matched the list, but was rejected as fitting to the benchmark.
- Scores were not a reason for the change.

## 6. What the results mean

- **Shown (verified):**
  - The model no longer makes a transcription factor's knockout lethal through iron import.
  - R2193 closes only when each pair loses at least one gene.
  - Growth and every other knockout are unchanged.
- **Inferred, not tested:**
  - that each pair can import iron on its own in *Y. lipolytica* (weaker for the D08 paralogs);
  - why the D08 pair does not compensate.
- **Method:**
  - The basis reset removed the heme report's stale-point artifact without fresh re-solves. The
    targeted check confirmed it: every closed reaction is at 0.
  - The violations the first (1e-9) check flagged are within GLPK's tolerance but not random. Most
    |S v| residuals equal 1e-6 × growth on a heme-route row: GLPK drops the 1e-6 heme a demand
    within tolerance (audit, verified with 13 LPs; Gurobi at 1e-9 gives the same growth). So the 1e-6
    general threshold has a blind spot. A heme-step knockout "rescued" this way would pass the check.
  - In this sweep every heme-route knockout returned 0, so no call is affected. Before reusing
    `scripts/checked_solve.py`, scale the residual tolerance to growth × the smallest biomass
    coefficient, or cross-check points with non-zero residuals on Gurobi (proposed, untested).

## 7. Limits, untested items and next gates

- **Condition.** Only the Task A medium + PO1f with GLPK.
- **Identity.** Orthology and protein detection, not *Y. lipolytica* biochemistry. For the D08 pair,
  the orthology is one-way.
- **Mechanism.** R2193 is a lumped Fe²⁺ uniport. Fet3/Ftr1 is an O₂- and Cu-dependent
  oxidase–permease that moves Fe³⁺ (established, uncited). The GPR follows the yeast-GEM convention;
  the mechanism is not modelled.
- **Unchanged and still open:**
  - iron species with no formula and charge 0;
  - R863 is gene-less (MRS3/4 ortholog YALI1_B28215g not added);
  - no FRE reduction step;
  - YALI1F07747g's deleted UniProt cross-reference.
- **Acceptance criterion 4** ("within 1e-9 of 0.3124263") fails literally by 4.4e-8, because the
  stated value was rounded. The intended criterion, WT unchanged between the models, holds exactly.
- **Budget overruns** (TASK.md):
  - 9 builds against 3;
  - the stopped first sweep (about 202 LPs plus 7 fresh re-solves);
  - 7 diagnosis LPs, 4 method-check rerun LPs and 1 v3 WT LP beyond the plan.
- **Not committed or pushed.** The superseded v1 and v2 builds were deleted at the user's request.

## 8. Independent audit and responses

The full audit is in [AUDIT.md](AUDIT.md); it covered v2.
- **Supported:** C1–C4, C6 and C8.
- **Partially supported:** C5, C7, C9 and C10.

Responses:

| Issue | Response |
| --- | --- |
| 1. Heme-skip blind spot | Accepted. Stated in §6, TASK.md and `checked_solve.py`. "Noise" replaced by "within GLPK tolerance, systematic". No call changes here; the fix is proposed, not applied (a method change for the user to decide). |
| 2. "Reciprocal best hit" wrong for D08 | Fixed. The data file, the R2193 notes (rebuilt as v3), REPORT §2 and the investigation report/table are corrected. |
| 3. Ye 2026 overstated | Fixed in the v3 notes, data and REPORT: co-expression with VHb, no Ftr1-only control. |
| 4. Method-check trap | Fixed. `knocked` is now required in `solve()`. `method_check.py` passes the closed set. [method_check_targeted.json](method_check_targeted.json) names R792 and flags the warm-start point. |
| 5. Narrowed re-solve rule | Accepted. It was not approved separately: a deviation from acceptance criterion 5, reported here and in TASK.md. No infeasible knockout was re-solved; the same 20 genes are infeasible in both models. |
| 6. Guard gaps | Fixed in part. Genes are created through cobra's context-aware path, notes are restored on context exit, and a new test covers it. Not changed: the reviewed SHA is not enforced (as in other steps), the old gene's auto-name is not pinned, and the post-edit sanity check stays. |
| 7. Provenance gaps | The full SHA of the experimental list is now in §5; `PYTHONHASHSEED` stays in commands.txt; the `flag_diagnosis.py` WT denominator (cf0825b8's, 1.4e-6 relative) is noted in TASK.md; diagnosis and sweep points differ numerically (vertex choice); noted. |
| 8. Wording | Criterion 4 is explained in §7. Positive coverage and the 11 unresolved positives are in §5. The names differ from the plan (`r2193_iron_gpr*`); noted. |
| 9. Budget | Updated in §7. |
| 10. Process state | The STATE.md entry is written. The full suite reran after the v3 code (below). |
| 11. Biology | Added to §7. |

**Full suite after v3:** 231 passed; identical test by test to the v2 run plus the new context test
(§4).
