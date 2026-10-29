# Cytosolic threonyl- and seryl-tRNA synthetase rules split from their mitochondrial paralogs (2026-10-29)

Scope and authorization: [TASK.md](TASK.md) (user, 2026-10-29: "Thr/Ser split first", then "Approve as written").
Evidence review: [../or_rule_review_r171_aars_20261029/REPORT.md](../or_rule_review_r171_aars_20261029/REPORT.md).
Audit: [AUDIT.md](AUDIT.md). Nothing committed or pushed.

## 1. What changed and why

| Reaction | Before (iYli21) | After |
| --- | --- | --- |
| `R759` threonyl-tRNA synthetase, cytosol | `YALI1A09121g or YALI1F14882g` | `YALI1F14882g` |
| `R760` threonyl-tRNA synthetase, mitochondrion | `YALI1A09121g or YALI1F14882g` | `YALI1A09121g` |
| `R723` seryl-tRNA synthetase, cytosol | `YALI1B07957g or YALI1F03908g` | `YALI1F03908g` |

`YALI1B07957g` stays in the model with no reaction. `R673` (prolyl) is unchanged. Each reaction also gets six
`aars_thr_ser_split_*` notes (curation id, status, authorization, previous rule, evidence, evidence tier).

**Why.** Each rule OR'd a cytosolic enzyme with its mitochondrial-type paralog, so neither gene's knockout could close
the cytosolic reaction. An OR member needs evidence that it catalyses this step alone in this compartment (AGENTS §5);
the mitochondrial-type paralogs have none for the cytosol. The recall gain is a consequence, not the reason.

**Evidence** (from the review; tagged there):
- *Verified:* UniProt names and locations: Q6C234 (`YALI1F14882g`) "Threonine–tRNA ligase, cytoplasmic", cytoplasm;
  Q6CHF9 (`YALI1A09121g`) mitochondrial matrix; Q6C350 (`YALI1F03908g`) "Serine–tRNA ligase, cytoplasmic", in the SES1
  UniRef50 cluster; Q6CFL2 (`YALI1B07957g`) DIA4-like, own cluster. Two of the locations are automatic (ARBA).
- *Inferred:* the compartment of each enzyme, from homology, length and N-terminal composition. No *Y. lipolytica*
  experimental localisation exists in the sources searched.

## 2. Code and data

| Item | File | SHA-256 |
| --- | --- | --- |
| Decision (data) | `model/curation/gpr/aars_thr_ser_split.json` | `0ed217f6…` |
| Guarded step | `platform/scripts/gem_annotate/aars_thr_ser_split.py`, `apply_aars_thr_ser_split` | `7dcb9633…` |
| Test | `platform/tests/test_aars_thr_ser_split.py` (6 tests, 21 subtests) | `34a765c4…` |
| Wiring | `cli.py` flag `--aars-thr-ser-split`; `main.py` call after the R2171 step | diff in git |

The step checks name, bounds, stoichiometry, species, annotation, rule and notes of all three reactions, and the
annotation and reaction links of all four genes, before any edit. It rejects a partly applied state, is idempotent,
and restores notes on cobra context exit. No gene is added, so the written order does not depend on the hash seed.
The flag has no dependency on another flag; the guards pin the before-state of the targets.

## 3. Build (verified)

Command (from `platform/`, worktree venv Python 3.13.5), full argv in [logs/build_commands.txt](logs/build_commands.txt):
the 28 opt-in flags of `965ab947` plus `--aars-thr-ser-split`.

| Build | Output SHA-256 |
| --- | --- |
| Control, without the flag | `965ab947f2aadc2c69443e7a794358c502acfdf0892e78582278ff07a37429d7`, byte-identical to the baseline |
| Candidate | `6abb6b3873d826e7e5b3105853ac2a1144d4be8e9b0b6ce8fbb69611d3c7c4e4`, [candidates/E5_…_r2171_1to1_aars_thr_ser_split_20261029.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_r2171_1to1_aars_thr_ser_split_20261029.xml) |

Both build logs contain the same libSBML error on setting the notes of `R304`. It is also in the R2171 task's
`logs/build_control_no_flag.log` and `logs/build_candidate_seed0.log` (grep, this session), so it predates this change
(not investigated here). The control's `.build.json` is copied to [logs/control_965ab947.build.json](logs/control_965ab947.build.json);
the control XML itself was a scratch file. `logs/build_commands.txt` gives the candidate as the shared argv plus the
extra flag; the full argv is in the candidate's `.build.json`.

## 4. Whole-model diff (verified)

[model_diff_vs_965ab947.json](model_diff_vs_965ab947.json), from the R2171 task's `model_diff.py` (copied unchanged):
**9 differences**, all intended: the rule and notes of `R759`, `R760`, `R723`, and the reaction links of
`YALI1A09121g` (loses `R759`), `YALI1F14882g` (loses `R760`) and `YALI1B07957g` (loses `R723`). No stoichiometry,
bound, species, annotation or other rule changed, so mass and charge balance are untouched.

## 5. Knockouts (verified)

Task A medium + PO1f, Gurobi 13.0.2 at 1e-9, reset per knockout, targeted check (`ko_compare_two.py` and
`checked_solve.py` from earlier tasks, unchanged); 2,164 LPs; [ko_compare_summary.json](ko_compare_summary.json),
[ko_compare.tsv](ko_compare.tsv). A first attempt with the worktree venv stopped before any solve (no scipy);
[logs/ko_compare_attempt1_noscipy.log](logs/ko_compare_attempt1_noscipy.log).

- WT 0.3116608216 on both models (residual 5e-14).
- Growth changes above 1e-6: `YALI1F14882g` and `YALI1F03908g`, 0.31166 → 0. No other gene.
- Class changes: those two genes, FN → TP, at every cutoff. No failing point; the same 20 unresolved knockouts on both
  models: 4 mtDNA genes and 16 nuclear genes whose knockout LP is infeasible, the same set on both models.

| Cutoff | 965ab947 TP / FN / FP / TN | 6abb6b38 TP / FN / FP / TN | Recall | Specificity | Precision | MCC |
| --- | --- | --- | --- | --- | --- | --- |
| 1 % | 94 / 221 / 35 / 699 | 96 / 219 / 35 / 699 | 0.298 → 0.305 | 0.952 → 0.952 | 0.729 → 0.733 | 0.350 → 0.356 |
| 5 % | 100 / 215 / 37 / 697 | 102 / 213 / 37 / 697 | 0.317 → 0.324 | 0.950 → 0.950 | 0.730 → 0.734 | 0.363 → 0.370 |
| 10 % | 104 / 211 / 42 / 692 | 106 / 209 / 42 / 692 | 0.330 → 0.337 | 0.943 → 0.943 | 0.712 → 0.716 | 0.361 → 0.368 |
| 15 % | 112 / 203 / 42 / 692 | 114 / 201 / 42 / 692 | 0.356 → 0.362 | 0.943 → 0.943 | 0.727 → 0.731 | 0.386 → 0.393 |

16 unresolved at each cutoff. Inputs: medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`), profile
`po1f_sd_leu.json` (`35307853…`), list `consensus_essential_genes.csv` (`1e887f5a…`). Positive coverage: 326 of the
1,612 listed genes are model genes; 315 are scored (TP + FN) and 11 are among the 16 unresolved knockouts; 1,286
listed positives have no model gene and are out of scope, not FN. Negatives are the user-defined class (model genes
not on the list); mtDNA genes are out of scope; the list was used in development, so these are development-reference
scores, not independent validation.

`YALI1A09121g` and `YALI1B07957g` stay FN, as expected: `R760` carries no flux at the optimum (dead-end mitochondrial
tRNA) and the model has no mitochondrial seryl-tRNA reaction. That is the category-d problem of the 2026-10-27 study,
not a rule problem.

## 6. pytest

Full suite, worktree venv, before vs after ([logs/pytest_compare.txt](logs/pytest_compare.txt)): 318 → 324 passed,
the same 19 failed and 20 errors, subtests 837 → 858. Test by test, the only differences are the 6 new tests, all
passing; no existing test changed outcome. The before run finished before any code was installed. The failures and
errors are pre-existing; the R2171 task recorded 19 failures and 14 errors, and the error count has since grown
without this change (not investigated).

## 7. Acceptance (TASK.md §4)

| # | Criterion | Result |
| ---: | --- | --- |
| 1 | Control reproduces `965ab947` byte for byte | met |
| 2 | Diff confined to the three rules, their notes and the four genes' links | met (9 differences; `YALI1F03908g` keeps `R723`, so only three genes change links) |
| 3 | Step fails closed, idempotent, context-clean | met (focused test) |
| 4 | Full pytest unchanged apart from the new tests | met (only the 6 new tests differ) |
| 5 | Sweep: exactly the two expected class changes at every cutoff, no failing point | met |

## 8. Verification levels and limits

- **Verified this session:** build identity, diff, WT, sweep and scores; the focused test.
- **Inferred:** the compartment assignment of each enzyme (no *Y. lipolytica* experiment).
- **Not tested:** media other than Task A; the prolyl pair; a combined build with the UM7 line (`5f54eb22`); the
  `R304` libSBML notes error.
- **Next gates:** the prolyl split after a presequence prediction; the R171 OR branches (review candidate 2); an
  annotation round for `YALI1A09121g`, which has no NCBI or UniProt cross-reference in the model.
- Large outputs: none above 1 MB. Departure from TASK.md §5: the plan put the sweep TSV in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/`; `ko_compare.tsv` (276 KB) stays in this folder instead, under the 1 MB rule.
- `scripts/ko_compare_two.py` and `scripts/model_diff.py` are unchanged copies from the R2171 task (SHA-256 `bb3118a6…`,
  `8a464ff8…`); the knockout script's docstring still names that task's models.
- Audit coverage limit: the stale-input test cases mutate R723 and `YALI1B07957g`; R759 / R760 field-level staleness
  goes through the same loop but is not exercised case by case.
