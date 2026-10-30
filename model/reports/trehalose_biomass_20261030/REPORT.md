# Trehalose removed from `biomass_C` (2026-10-30)

**Outcome.** Candidate
[`E5_…_aco_split_treh_biomass_v2_20261030.xml`](../../candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_v2_20261030.xml),
SHA-256 `164542fb6056f01e257ecfaed1b456961e8a07e46e9e4e1942f79489d5014beb` = `e1e7c868` + `--trehalose-biomass-removal`.
`biomass_C` no longer drains trehalose `m294[C_cy]` (was 0.00686 mmol/gDW); nothing else in the model changed.
The three genes on the trehalose-6-phosphate synthase reaction R143 go from essential to non-essential at every
cutoff (FP → TN); no other gene changes class. Not a designated reference. Uncommitted.

The first build, `d2e4c0c4b02c7e48…` (`…_treh_biomass_20261030.xml`), was swept, drawn and audited; after the audit
two evidence notes were reworded (§6) and the model rebuilt as v2. v2 differs from `d2e4c0c4` only in those two
`biomass_C` note texts ([diff](model_diff_v2_vs_d2e4c0c4.json)), so the sweep and the figure (run on `d2e4c0c4`)
apply to v2 unchanged. `d2e4c0c4` is kept until the user decides.

**Authorization** (TASK.md §0): the user, 2026-10-30, "Approve C on e1e7c868" (option C as written in TASK.md §9,
baseline `e1e7c868`, §9 budget; option D deferred; no commit or push).

## 1. What changed and why

| | Before (`e1e7c868`) | After (`d2e4c0c4`) |
| --- | --- | --- |
| `biomass_C` trehalose `m294[C_cy]` | −0.00686 mmol/gDW (6.86 nmol/mg DW, 2.35 mg/gDW) | absent |
| `biomass_C` notes | 31 keys | + 7 `trehalose_biomass_*` keys (id, status, authorization, previous term, evidence, tier, scope) |

**Why** (evidence tiers: **verified** = read this session; **inferred**; **recalled**):
- **Measured content** (verified, full text of Flores et al. 2011, *PLoS One*,
  [PMID 21931609](https://europepmc.org/article/MED/21931609), PMC3171402): wild-type *Y. lipolytica* (CJM 645)
  contains < 1 nmol trehalose/mg dry weight on 0.17 % YNB without amino acids and ammonium sulfate + 0.1 % glutamate,
  2 % glucose, 30 °C (Table 4, three cultures; enzymatic assay after boiling-water extraction), and < 1 "grown in
  glucose up to stationary phase or in glycerol". The model's 6.86 nmol/mg is more than 6.8 times that upper bound.
- **Phenotype** (verified, same paper): the *tps1* disruptant grows on glucose at least as fast as the WT (doubling
  time 139 ± 4 vs 149 ± 4 min, four experiments), and no trehalose is detected in `nth1 tps1` or in heat-shocked
  `tps1` (the single mutant's content at 30 °C is not reported). So trehalose is not an obligate biomass component in
  this condition. The paper reports "< 1" and "not detected" without stating a detection limit.
- **The old value has no source** (verified locally): 0.00686 sits in the immutable `model/source/iyali26.xml`; it is
  hard-coded without a comment in the Debiaggi 2025 eciYali5-GEM script (`xCARBOHYDRATE`, line 56), whose four
  carbohydrate numbers are iYali4's `xBIOMASS` terms (trehalose 0.0032) scaled 2.134–2.144×. Where iYali4's 0.0032
  came from is unknown.
- **Identity** (verified): `YALI1E17616g` = `YALI0E14685g` = UniProt [O74932](https://rest.uniprot.org/uniprotkb/O74932)
  TPS1, EMBL AJ011032 (the cDNA deposited by Flores 2011).

The change was chosen on the measurement and the phenotype, not on the match rate (TASK.md §6 compares keep / lower /
remove; lowering to any positive value would keep the three genes essential).

## 2. Implementation (the pipeline)

- **Data:** [`model/curation/biomass/trehalose_biomass.json`](../../curation/biomass/trehalose_biomass.json),
  written by [scripts/write_curation.py](scripts/write_curation.py) from the SHA-checked baseline: the term, the
  species identity and reaction context before/after, the full accepted `biomass_C` composition without trehalose
  (61 entries), its bounds and 31 note keys, the notes to add, evidence with URLs and tiers, reviewer and the
  authorization.
- **Step:** `platform/scripts/gem_annotate/trehalose_biomass.py`, `apply_trehalose_biomass`. Before any edit it
  requires: a valid spec; `biomass_C` the only objective (coefficient 1), maximised; `m294[C_cy]` with the recorded
  name, formula, charge and compartment; every other `biomass_C` coefficient equal to the reviewed composition
  (`math.isclose`, rel 1e-12) and the same bounds; the term at −0.00686 (same comparison) when present; the species' reaction
  set equal to the before or after list; no note outside this curation added or removed (keys only: the values of
  other notes are not pinned, as in the dolichol step); no conflicting value of this curation's notes; and either term present without the notes, or term absent with all notes. Then it sets the coefficient to 0
  with `add_metabolites(..., combine=False)` (cobra drops zero entries and registers the undo in a model context),
  restores the notes on context exit, adds the notes, and re-checks that only the term left. A rerun, also on a saved
  and reloaded model, returns `already_correct`.
- **Wiring:** flag `--trehalose-biomass-removal` (`cli.py`), called in `main.py` right after the aconitase step,
  which is after every other step that edits `biomass_C` (TPP, heme/iron, cofactor terms, dolichol route) and after
  the E5 energy candidate; nothing after it edits `biomass_C`.
- **Test:** `platform/tests/test_trehalose_biomass.py` (5 tests, 19 subtests, no LPs): 13 stale or partial states
  (the objective-direction case added after the audit)
  fail closed without edits; 6 scope edits of the spec are rejected; the applied model equals the baseline with only
  the term and notes changed, survives an SBML round trip and stays idempotent; context exit restores the baseline;
  the flag is opt-in. Two test-side fixes during development: an expected message containing `[C_cy]` is passed
  through `re.escape` (it was read as a character class), and the solver is flushed before the snapshot in the
  "biomass_C removed" case (optlang refreshes the objective lazily; the step's objective check triggered the refresh,
  not an edit).

## 3. Verification

| # | Acceptance criterion (TASK.md §8) | Result |
| --- | --- | --- |
| 1 | Control (flag off) reproduces `e1e7c868` | met twice: before and after the audit fixes, `e1e7c8686c1fbbbb…` byte for byte ([build_results](logs/build_results.txt), [v2](logs/build_results_v2.txt), [commands](logs/build_commands.txt)) |
| 2 | Seeds 0 and 1 identical | met: `d2e4c0c4…` (first build) and `164542fb…` (v2) at both seeds |
| 3 | Whole-model diff only the term and the notes | met: v2 vs `e1e7c868` 2 differences, `biomass_C` stoichiometry (`m294[C_cy]` −0.00686 removed) and notes (7 keys added, none changed or removed) ([model_diff_v2_vs_e1e7c868.json](model_diff_v2_vs_e1e7c868.json); first build [model_diff_vs_e1e7c868.json](model_diff_vs_e1e7c868.json)); the auditor confirmed it from a raw XML diff (8 lines) |
| 4 | Every other `biomass_C` coefficient kept; no other reaction touched | met (same diff) |
| 5 | Checked sweep, four cutoffs | met on `d2e4c0c4` (§4), the same optimisation problem as v2 |
| 6 | Full pytest: only the new tests differ | met for the first version: 409 → 414 passed, the same 19 failed and 20 errors (identical IDs), subtests 960 → 978 ([pytest_compare.txt](logs/pytest_compare.txt)). After the audit fixes, again: 409 → 414 passed, the same 19 failed and 20 errors (identical IDs), subtests 960 → 979 (+19 = the new subtests) ([pytest_compare2.txt](logs/pytest_compare2.txt)) |
| 7 | Independent audit | no unsupported claim (§6) |

Tree state: `7aa4baf` + this round's uncommitted files ([build_head.txt](logs/build_head.txt) for the first
builds, [build_head_v2.txt](logs/build_head_v2.txt) for v2).

## 4. Knockout sweep

Task A medium (`6eca70fd…`) + PO1f (`35307853…`), Gurobi FeasibilityTol/OptimalityTol 1e-9 with the real HOME, basis
reset per knockout, knocked-reaction flux check; SHA-checked short copies of both models; 2,164 LPs, no fresh
re-solves needed, no point failing a check
([ko_compare.tsv](ko_compare.tsv), [ko_compare_summary.json](ko_compare_summary.json),
[script](scripts/run_sweep.zsh)). The same 20 genes are non-optimal in both models (16 counted as unresolved in the
scores); none of them is a trehalose gene.

- **WT growth:** 0.3116608 → 0.3122964 (+0.20 %): the 0.00686 mmol/gDW of trehalose (two glucose units each, via
  G6P and UDP-glucose) no longer has to be made.
- **Class changes** (all three are user-defined negatives, i.e. not on the consensus list):

  | Gene | KO growth before → after | 1 % | 5 % | 10 % | 15 % |
  | --- | --- | --- | --- | --- | --- |
  | `YALI1E17616g` (TPS1) | 0.0 → 0.31230 | FP → TN | FP → TN | FP → TN | FP → TN |
  | `YALI1D17620g` ("TPS2", also R767) | 0.0 → 0.31230 | FP → TN | FP → TN | FP → TN | FP → TN |
  | `YALI1E36568g` ("TPS3") | 0.0 → 0.31230 | FP → TN | FP → TN | FP → TN | FP → TN |

- **No other gene changes class at any cutoff.** The KO/WT ratios closest to a cutoff moved by less than 0.001 and
  did not cross: 0.15946 → 0.15913 (15 %, 9 genes), 0.10970 → 0.10962 (10 %, 12 genes), 0.05493 → 0.05482 (5 %,
  `YALI1E34269g`).
- **Scores** (negatives user-defined; development reference, not independent validation):

  | Cutoff | Before TP / FN / FP / TN | After | Precision | Specificity | MCC |
  | --- | --- | --- | --- | --- | --- |
  | 1 % | 97 / 218 / 35 / 699 | 97 / 218 / 32 / 702 | 0.735 → 0.752 | 0.952 → 0.956 | 0.360 → 0.369 |
  | 5 % | 104 / 211 / 37 / 697 | 104 / 211 / 34 / 700 | 0.738 → 0.754 | 0.950 → 0.954 | 0.376 → 0.385 |
  | 10 % | 108 / 207 / 42 / 692 | 108 / 207 / 39 / 695 | 0.720 → 0.735 | 0.943 → 0.947 | 0.374 → 0.383 |
  | 15 % | 125 / 190 / 42 / 692 | 125 / 190 / 39 / 695 | 0.749 → 0.762 | 0.943 → 0.947 | 0.425 → 0.434 |

  Recall is unchanged (0.308 / 0.330 / 0.343 / 0.397); 326 in-model positives, 16 unresolved, in both models.

## 5. Figure

[escher/README.md](escher/README.md): WT before vs after
([PNG](escher/trehalose_wt_before_vs_after.png)) shows 0.002138 = 0.00686 × 0.3116608 through R143 and R767 into the
trehalose term before, and no flux after; the TPS1 knockout ([PNG](escher/trehalose_tps1_ko_before_vs_after.png))
grows 0 before and at WT after. 12 LPs.

## 6. Independent audit

A fresh read-only reviewer got TASK.md, the code, data, test, models, logs and figures (no conversation history) and
checked 9 claims ([AUDIT.md](AUDIT.md)): 7 supported (hashes and argv; the diff, also from a raw XML diff; the test
and the two test fixes; pytest; the sweep, recomputed from the TSV; the justification and limits; the figure
numbers), 2 partly supported (the guard; the evidence wording). No unsupported claim. What was fixed:

- **Evidence wording, in the curation notes (rebuilt as v2):** "below the detection limit" → "reported as below
  1 nmol/mg DW; the paper states no detection limit"; "grows at a wild-type rate and makes no detectable trehalose"
  → "grows at least as fast as the WT; no trehalose detected in tps1 with nth1 disrupted or after heat shock (the
  single mutant at 30 °C not reported)"; the medium attributed to Methods, not Table 4; "stationary phase" quoted as
  the paper says it.
- **Evidence wording, in TASK.md:** the TEF1-YlTPS1 activity rise was in a strain that keeps TPS2/TPS3 (CJM 613, a
  different background from CJM 645), so "alone" was dropped and it no longer supports option D.
- **Guard:** objective direction now checked (and tested); the coefficient uses the same `isclose` as the
  composition; the docstring and this report say that other notes are pinned by key only.
- **Figure:** the drain label no longer says "TERM REMOVED" on the before panels; the README's UDP-glucose numbers
  corrected (R638/R793 −0.0015, R4 +0.0006).

Noted, not changed: UniProt O74932 is 469 aa while Flores reports a 486-aa ORF from the cDNA; the identity rests on
the EMBL AJ011032 cross-reference, and the length difference is unexplained (not checked). The auditor wrote three
scratch files in the session scratchpad (outside the repo) although asked not to.

## 7. Limits and what was not tested

- **The measurement is an upper bound**, not a point value; the real content in our condition is between 0 and
  about 1 nmol/mg DW (≤ 0.34 mg/gDW), so removing the term drops at most that real mass. `biomass_C` is not
  re-normalised (earlier terms were also added without re-normalisation); the 2.35 mg/gDW removed is 0.23 % of DW.
- **Strain and medium differ:** Flores used CJM 645 (MATB met6) on YNB-glutamate glucose; the benchmark uses PO1f on
  the Task A medium. That PO1f contains as little trehalose is assumed, not measured.
- **Condition dependence:** heat shock (40 °C) raises trehalose to 7–20 nmol/mg DW, and *tps1* grows poorly at 35 °C;
  neither condition is modelled.
- **R143 rule (option D) untouched:** still `TPS2 and TPS1 and TPS3`; with no trehalose demand it no longer affects
  the score. Its own round, if wanted.
- **Non-objective pseudoreactions** `xBIOMASS`, `newBiom`, `R1372` keep their trehalose terms (out of scope).
- **Cas12a-only essentiality of TPS1:** the *TFC1* promoter overlap (Flores 2011) is a candidate explanation
  (inferred; guide positions not checked).
- Scores use the development consensus list, which was used for curation; they are not independent validation.

## 8. Math and method notes

- The term is a pure growth-dilution drain: with one route (R143 → R767), any positive coefficient makes every gene
  of that route essential and the KO ratio exactly 0. That is why option B (lowering) could not reconcile the model
  with the growing *tps1* mutant, whatever value is chosen.
- A firmer test of option C would be trehalose measured in PO1f on the Task A medium (or the fermentor runs), and a
  `tps1` strain grown there; neither is available.
