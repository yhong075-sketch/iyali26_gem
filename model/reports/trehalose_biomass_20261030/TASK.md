# Trehalose in `biomass_C`: the 0.00686 mmol/gDW term against the measured *Y. lipolytica* content (2026-10-30)

**Status: approved 2026-10-30, option C on baseline `e1e7c868`** (§0). Written before approval: this file, [before.json](before.json) (read-only dump
of the baseline, [scripts/before_state.py](scripts/before_state.py)) and [evidence/](evidence/) (UniProt rows and the
screen calls, with [SHA256SUMS.txt](evidence/SHA256SUMS.txt)). The full text of Flores 2011 as fetched is in
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/trehalose_biomass_20261030/flores2011/` (Europe PMC XML, SHA-256
`c5529d32…`, and the extracted text).

## 0. Authorization

- **Brief** from the session that started this one (2026-10-30): run one `curation-change` round on the trehalose
  false-positive group up to a TASK.md draft; lay out keep / lower / remove for the biomass term, with the R143
  AND-rule kept as a separate option; then stop for the user's choice. "Do not change the model to improve the match
  rate; the change must rest on the measured composition."
- **User, 2026-10-30**, answering the three questions of §7 with the next-step button:

  > "Approve C on e1e7c868"

  The button stood for: implement option C as written in §9 on baseline `e1e7c868`, with the §9 compute budget;
  option D deferred to its own round; no commit or push.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Worktree / branch / HEAD | `claude/trehalose-biomass-24e3c1` @ `7aa4baf` (= `origin/main` after `git fetch` today, the PR #39 merge); clean before this folder |
| Earlier trehalose work | none: no `trehal`/`R143`/`R767`/`m294` entry in `model/STATE.md`, no matching branch in `git branch -a`, no curation file. Only related history: `9ab337a` removed `YALI1E07744g` mis-annotated as a trehalase (it is a transketolase), which does not touch these reactions |
| **Proposed baseline** | `candidates/E5_…_r673_prolyl_split_aco_split_20261029.xml`, SHA-256 `e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897` (re-hashed). Newest candidate on `origin/main` now that PR #39 is merged. **Not a designated reference** (`model/candidates/README.md`) |
| Brief's baseline | `059316ef34417c88aadba61cfd5606134d27879a842b3835d2cc7df917a9248d` (re-hashed). `e1e7c868` = `059316ef` + `--aconitase-aco2-split`; the before-state of every trehalose target is identical in the two files (both dumps compared this session) |
| Baseline build argv | in `e1e7c868`'s `.build.json` (36 flags ending `--r673-prolyl-split --aconitase-aco2-split`) |
| Experimental list | `consensus_essential_genes.csv`, SHA-256 `1e887f5a…` (re-hashed). None of the four trehalose genes is listed → all are user-defined negatives |
| Medium / profile | Task A `6eca70fd…` + PO1f `35307853…` (per `model/expected/benchmark.md`) |

The baseline choice is the user's (§7, question 1).

## 2. Model state of the targets ([before.json](before.json))

| Reaction | Name in model | Bounds | Rule | Role |
| --- | --- | --- | --- | --- |
| `biomass_C` (objective) | biomass_C | (0, 1000) | — | consumes `m294[C_cy]` trehalose at **0.00686 mmol/gDW** = 6.86 nmol/mg DW = 2.35 mg/gDW (0.23 % of DW) |
| R143 | α,α-trehalose-phosphate synthase (UDP-forming) | (0, 1000) | `YALI1D17620g and YALI1E17616g and YALI1E36568g` | UDP-glucose + G6P → UDP + T6P; **sole producer of T6P `m297`** |
| R767 | trehalose-phosphatase | (0, 1000) | `YALI1D17620g` | T6P + H2O → trehalose + Pi; **sole producer of cytosolic trehalose** |
| R142 | α,α-trehalase | (0, 1000) | `YALI1D19108g` | trehalose + H2O → 2 glucose (context, not changed) |
| R1349 | trehalose vacuolar transport | (0, 0) | — | closed (context) |
| R1013 | α,α-trehalose exchange | (0, 1000) | — | export only; no reaction imports trehalose (context) |

Directional producers/consumers (static, from bounds and signs): `m294[C_cy]` is produced only by R767; `m297` only by
R143 and consumed only by R767. Trehalose is also consumed by the non-objective pseudoreactions `xBIOMASS` (0.0032),
`newBiom` (0.02814) and `R1372` (7.1e-7); the yeast 6/8 pseudoreactions R1387/R1710 (0.0234) are closed. **These are
out of scope and stay unchanged.** `biomass_C` has no note about trehalose.

So in the baseline, knocking out any one of the three R143 genes stops trehalose synthesis and therefore growth.
**Current knockout class** (recalled from the aconitase round's Gurobi 1e-9 sweep on both `059316ef` and
`e1e7c868`, [ko_compare.tsv](../aconitase_aco2_split_20261029/ko_compare.tsv); not re-solved): KO growth 0.0 for
`YALI1E17616g`, `YALI1D17620g`, `YALI1E36568g` → essential at every cutoff → **FP**; `YALI1D19108g` KO growth = WT
→ TN. 15 % on `e1e7c868`: 125 TP / 190 FN / 42 FP / 692 TN.

## 3. Gene identity

| Model gene | CLIB122 locus (locus map, verified) | UniProt (YALI0 entry, fetched this session) | Verified name | Function / evidence level |
| --- | --- | --- | --- | --- |
| `YALI1E17616g` | `YALI0E14685g` | [O74932](https://rest.uniprot.org/uniprotkb/O74932), TPS1, 469 aa, EMBL **AJ011032**, cites PMID 21931609 | **YlTPS1** (Flores 2011 deposited cDNA AJ011032) | T6P synthase; direct *Y. lipolytica* genetics and enzymology (Flores 2011) |
| `YALI1D17620g` | `YALI0D14476g` | [Q6C939](https://rest.uniprot.org/uniprotkb/Q6C939), uncharacterised, 800 aa | none verified; called "YlTPS2" in Flores 2011 by homology only | putative T6P phosphatase; Génolevures annotation "similar to *S. cerevisiae* TPS2"; mRNA measured, no mutant |
| `YALI1E36568g` | `YALI0E31086g` | [Q6C3Z3](https://rest.uniprot.org/uniprotkb/Q6C3Z3), uncharacterised, 1049 aa | none verified; "YlTPS3" in Flores 2011 | putative regulatory subunit (closest to Sc TPS3/TSL1); disruption mutant made (Flores 2011) |
| `YALI1D19108g` | `YALI0D15598g` | [Q6C8Z4](https://rest.uniprot.org/uniprotkb/Q6C8Z4), trehalase (automatic) | none verified; "YlNTH1" in Flores 2011 | putative neutral trehalase; disruption raises trehalose (Flores 2011) |

Screen calls ([evidence/screen_calls.csv](evidence/screen_calls.csv), verified): `YALI1E17616g` Cas9 FS +0.57
non-essential, Cas12a FS −1.95, q 0.025 **essential**; `YALI1D17620g` and `YALI1E36568g` non-essential in both;
`YALI1D19108g` non-essential in both. None is on the consensus list.

## 4. Flores et al. 2011 (read in full this session, **verified**)

Flores CL, Gancedo C, Petit T. *PLoS One* 2011, [PMID 21931609](https://europepmc.org/article/MED/21931609),
[PMC3171402](https://europepmc.org/article/PMC/PMC3171402), doi 10.1371/journal.pone.0023695. Full text fetched from
the Europe PMC REST full-text XML.

- **Medium and growth:** synthetic medium, 0.17 % YNB without amino acids and ammonium sulfate, 0.1 % glutamate,
  pH 6, auxotrophic supplements 20 µg/ml, 2 % glucose ("generally"), shaken at 30 °C (Methods). Wild type for
  Table 4 is CJM 645 (MATB met6).
- **Measurement:** boiling-water extraction, Megazyme enzymatic trehalose kit; unit **nmol/mg dry weight**
  (= µmol/gDW; 1 nmol/mg = 0.001 mmol/gDW). Table 4, three independent cultures, glucose medium: WT **< 1**;
  `TEF1p-YlTPS1` 8.9 ± 2.1; `nth1` 18 ± 2; `nth1 tps1` < 1. Text: WT "grown in glucose up to stationary phase or in
  glycerol" < 1. Heat shock (40 °C, 2 h) raised trehalose to 7–20 depending on strain. The growth phase of the
  Table 4 cultures is not stated (only "grown in minimal medium glucose").
- **tps1 phenotype:** the disruptant (cassette at nt 710 of the ORF; 5′ insertions failed, attributed to the
  overlapping, likely essential *YlTFC1* promoter) **grew on glucose**, doubling time 139 ± 4 min vs WT 149 ± 4 min
  (four experiments); no trehalose in `nth1 tps1` or in heat-shocked `tps1`. Growth at 35 °C was drastically slowed;
  sporulation reduced.
- **TPS2 / TPS3 and T6P synthesis:** YlTPS1 is "the only one active in trehalose synthesis when expressed in
  *S. cerevisiae*" among the Tps1-like sequences; TEF1-driven YlTPS1 (CJM 613, which keeps its own TPS2/TPS3) showed
  55 mU/mg protein of T6P synthase activity vs < 5 in CJM 645, a different background (corrected after the audit: an
  earlier draft said "alone"). The YlTPS3 cDNA did not complement *S. cerevisiae tps1*; `tps3` disruption abolished trehalose
  accumulation **upon heat shock** (basal content of `tps3` not reported); the authors suggest Tps3 stabilises the
  synthase complex during heat shock. YlTPS2 is named by homology only; no *tps2* mutant and no phosphatase assay.
  Nothing in the paper shows that Tps2 or Tps3 is required for T6P synthesis at 30 °C.

## 5. Where the 0.00686 came from (verified locally this session)

- It is in the immutable build input `model/source/iyali26.xml` (`biomass_C`, line 34795) and in `iyli21.xml`; no
  builder step touches it (`grep` of `platform/scripts/gem_annotate/`: no `m294` or trehalose edit of `biomass_C`).
- It matches the carbohydrate pseudoreaction hard-coded in the Debiaggi 2025 eciYali5-GEM build script
  (`code/debiaggi2025/generate_iYali5GEM.m` line 56, local copy under
  `$IYALI26_RESEARCH_ROOT/reference/external_models/chalmers_yarrowia/eciYali5-GEM/`, upstream
  [bioengtaltech/eciYali5-GEM](https://github.com/bioengtaltech/eciYali5-GEM), not opened online):
  `xCARBOHYDRATE` = m294 0.00686, chitin m401 0.868358, glucan m1123 0.943397, mannan m1324 0.235849 — the same four
  numbers as `biomass_C`. The script gives no source; it later rescales the block to close biomass mass (the final
  iYali5.yml has m294 0.00862).
- Those four numbers are iYali4's `xBIOMASS` carbohydrate terms (trehalose 0.0032, glucan 0.4415, mannan 0.1104,
  chitin 0.4068, local `sourceModels/iYali4.xml`, note "Modeling reaction" only) multiplied by 2.134–2.144×
  (computed this session). So the trehalose share is **inherited from iYali4, scaled, never measured for this model**.
  iYali4 is Kerkhoven et al. 2016 (*npj Syst Biol Appl*, doi 10.1038/npjsba.2016.5, **recalled, not opened**); where
  iYali4's 0.0032 came from is **unknown**.
- For comparison, the *S. cerevisiae* yeast 6/8 pseudoreactions in this model carry 0.0234 (verified in the model).

So the model term (6.86 nmol/mg DW) is more than 6.8× the measured *Y. lipolytica* upper bound (< 1 nmol/mg DW) on
glucose, and the only source behind it is an unexplained scaffold value.

## 6. Options

### Biomass term (pick one)

| | Change | Evidence for | Evidence against / limits | Expected effect (to be measured) |
| --- | --- | --- | --- | --- |
| **A. Keep** | none | none found for 0.00686 in *Y. lipolytica* | contradicts Table 4 (< 1) by more than 6.8×; makes TPS1 essential although the *tps1* mutant grows on glucose (Flores 2011) | no change; 3 FP stay |
| **B. Lower** to a sourced value | coefficient → ≤ 0.001 mmol/gDW | Table 4 gives an upper bound, < 1 nmol/mg DW | no measured point value exists: 0.001 is the reported bound (the paper states no detection limit), any number below it is a choice. Any positive coefficient keeps trehalose an obligate product, so the model would still kill the *tps1* mutant that Flores shows growing | the three genes stay essential (any positive drain with a single route); WT growth rises slightly |
| **C. Remove** (recommended) | delete `m294[C_cy]` from `biomass_C`; add curation notes on `biomass_C`; nothing else | (1) WT content < 1 nmol/mg DW on glucose, in stationary phase and on glycerol (Table 4, text); (2) the *tps1* disruptant, which makes no detectable trehalose with *nth1* also disrupted or after heat shock, grows on glucose at 30 °C at least as fast as the WT (139 vs 149 min), so trehalose is not an obligate biomass component in this condition (wording corrected after the audit); (3) the 0.00686 has no source (§5) | the term removed is 2.35 mg/gDW; the real trehalose it stood for is below the reported bound (< 1 nmol/mg DW = 0.34 mg/gDW; the paper states no detection limit), so a small real mass is dropped; `biomass_C` is not re-normalised (later terms were added without re-normalisation either). Flores used CJM 645 (MATB met6, from A1-5 × CJM 364), not PO1f, and YNB + glutamate, not Task A; heat stress (not modelled) raises trehalose to 7–20 | the R143 and R767 reactions no longer need flux; the three genes are expected to go essential → non-essential (FP → TN) at every cutoff; WT growth rises slightly; other genes near a cutoff could move with WT and must be checked |

The C sub-variant "remove and re-scale glucan/mannan/chitin to keep the mass" is **not** proposed: no *Y. lipolytica*
carbohydrate measurement backs a particular re-scaling, and it would change three more coefficients.

### R143 rule (separate decision, separate round)

**D. R143 → `YALI1E17616g` alone**; R767 keeps `YALI1D17620g`; `YALI1E36568g` stays in the model with no reaction.
- For: YlTPS1 is the only Tps1-like *Y. lipolytica* sequence active in trehalose synthesis when expressed in
  *S. cerevisiae*, and YlTPS3 does not complement *S. cerevisiae tps1* (Flores 2011, verified). (The TEF1-YlTPS1
  activity rise was in a strain that keeps TPS2/TPS3, so it does not show that Tps1 works without them; corrected
  after the audit.)
  AGENTS.md §5 asks for evidence of a complex or co-requirement for each AND member; for the T6P step at 30 °C there
  is none for Tps2 or Tps3. That *S. cerevisiae* Tps1 makes T6P without Tps2 and that Tps1–Tps2–Tps3/Tsl1 form a
  complex is **established/uncited** here (not opened this session).
- Against / limits: `tps3` abolished heat-shock accumulation, so Tps3 has a real role in one condition; the AND rule
  encodes complex membership, which this model uses elsewhere for complexes; Tps2's phosphatase role is homology only.
- Effect: after C, none on the essentiality score (trehalose is no longer required). Without C, `YALI1E36568g` would
  go FP → TN, `YALI1D17620g` would stay essential through R767 and `YALI1E17616g` would stay essential.
- One change per run (skill rule), so D is not part of this round. If approved, it gets its own TASK.md, before C or
  after it.

## 7. Questions for the user

1. **Baseline:** `e1e7c868` (proposed: newest on `origin/main`; same trehalose state) or `059316ef` (the brief's)?
2. **Biomass option:** A, B (with which value) or C (recommended)?
3. **R143 option D:** defer to its own round (proposed), or drop?

## 8. Acceptance criteria (option C)

1. Control build with the `e1e7c868` argv (flag off), `PYTHONHASHSEED=0`, `/opt/miniconda3/bin/python`, reproduces
   `e1e7c868…` byte for byte.
2. Candidate seeds 0 and 1 give identical bytes.
3. Whole-model diff vs the baseline: **only** the `m294[C_cy]` entry of `biomass_C` and the new `biomass_C` notes.
   No other reaction, metabolite, gene, rule, bound, formula, charge or annotation differs.
4. `biomass_C` keeps every other coefficient exactly; mass and charge balance of no other reaction changes (no other
   reaction touched).
5. Checked sweep, Task A medium + PO1f, Gurobi FeasibilityTol/OptimalityTol 1e-9 with the real HOME, basis reset per
   knockout, knocked-reaction flux check: WT status optimal, new WT growth recorded; the three R143 genes' KO/WT
   measured; every other gene's class compared at **1, 5, 10 and 15 %**, with any change listed and explained; no
   failed or non-finite solve, or each reported separately. Scores with recall, specificity, precision and positive
   coverage; negatives user-defined; development reference, not independent validation.
6. Full pytest: only the new tests differ from the before run.
7. Independent read-only audit: no unsupported claim left unfixed or unreported.

## 9. Plan (after approval of C)

1. **pytest before** at `7aa4baf`, full run in the background (about 13 min).
2. **Data:** `model/curation/biomass/trehalose_biomass.json`: target `biomass_C` / `m294[C_cy]`, previous coefficient
   −0.00686, new state "absent", origin trace (§5), evidence (Flores 2011 URL, table and text locations, tier
   `direct_Y_lipolytica_measurement_upper_bound` plus phenotype), reviewer, the user's authorization quote.
3. **Step:** `platform/scripts/gem_annotate/trehalose_biomass.py`, `apply_trehalose_biomass`, loaded with
   `MODEL.curation_file("trehalose_biomass.json")`. Fails closed unless `biomass_C` exists, is the objective and holds
   `m294[C_cy]` at exactly the data's previous coefficient; touches only that entry and `biomass_C` notes;
   idempotent (re-run on an applied model is a no-op that still checks the notes). Flag `--trehalose-biomass-removal`
   in `cli.py`, called in `main.py` after the aconitase step (every earlier `biomass_C` guard has run by then; this
   guard pins only its own term).
4. **Test:** `platform/tests/test_trehalose_biomass.py` on an in-memory model: applies; idempotent; a stale
   coefficient, a missing `m294[C_cy]` and a non-objective `biomass_C` each fail closed; nothing outside the term
   and notes changes.
5. **Builds** (from `platform/`, `PYTHONHASHSEED=0`): control to the scratchpad; candidate
   `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_20261030.xml`
   (158-byte name, 165 with `.build.json`, under the 255-byte limit) = baseline argv + `--trehalose-biomass-removal`;
   seed-1 rebuild to the scratchpad.
6. **Diff:** `model_diff.py` from the aconitase round (copied unchanged).
7. **Sweep:** `ko_compare_two.py` and `checked_solve.py` (aconitase round copies) on SHA-checked short-path copies of
   both models, conditions as in §8.5.
8. **Audit** by a fresh read-only reviewer given only this file, the diff, the logs and the deliverables; then
   REPORT.md, a `model/STATE.md` entry and a `model/candidates/README.md` row.
9. **Figure:** Escher before/after (skill `gem-route-escher-maps`) of the trehalose branch (G6P / UDP-glucose → R143
   → R767 → `biomass_C`, with R142) for WT and the `YALI1E17616g` knockout, pFBA on both models.

**Budget:** about 2,164 LPs for the two-model sweep (1,081 genes + WT each) + about 20 fresh re-solves + about 10 for
the figure; about 45 minutes of compute; no cluster job. **Stop** on any failed acceptance criterion and report it.
No commit or push without the user's explicit instruction.

## 10. Not in this round

- Option D (R143 rule), see §6.
- Trehalose terms in the non-objective `xBIOMASS`, `newBiom` and `R1372`.
- Condition-dependent trehalose (heat shock, 35 °C growth): not modelled.
- Why `YALI1E17616g` is essential in Cas12a only: the *TFC1* promoter overlap (Flores 2011) is a candidate
  explanation (**inferred**; guide positions not checked).
