# TASK: open the cholinephosphotransferase R1726 (Kennedy CDP-choline step), keep R1725 closed (2026-10-27)

## 0. Authorization

The user, 2026-10-27, after reading the Stage 1 draft of this plan (shared notes folder,
`kennedy_closures_TASK_draft_20261027.md`; its section 12 gives the full approval sentence):

> Approve Option A

"Option A" is section 6 of the draft: R1726 (0, 0) → (0, 1000) with its GPR unchanged, R1725 kept closed
with the reason recorded, through a guarded `--r1726-cpt-open` step after `--cofactor-biomass-terms` on
00e8e397, with the test, builds, diff, Gurobi knockout sweep, re-census, review, audit, REPORT and STATE
entry of sections 7–10. No commit or push. Sections 1–12 below are the draft unchanged apart from this
section, the status line and the paths of this folder.

**Status: first run stopped (STOP.md); re-run authorized and completed 2026-10-27 (section 13, REPORT.md, AUDIT.md).** Stage 1 outputs: [stage1_kennedy.json](stage1_kennedy.json),
[stage1_controls.json](stage1_controls.json), scripts in [scripts/](scripts/), the pFBA difference (417 kB)
in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/kennedy_closures_20261027/stage1_pfba_diff.json`.

## 1. Starting state (read 2026-10-27)

- **Branch / HEAD:** `claude/kennedy-closures-864c91` at `31e7aa4`, clean. **PR #24 is merged:**
  `origin/main` = `e26ac81` (merge of `31e7aa4`, 2026-10-27 16:18), and its tree is identical to this HEAD
  (`git diff HEAD origin/main` is empty). So the baseline below is on `main`.
- **Baseline model:** `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_20261027.xml`,
  SHA-256 `00e8e39775b64cfeb79cfafff65ef82eaf4760eff1b64fcddc3ce93b40919592` (re-hashed; the context
  loader reports the same canonical SHA). Built with the a1237821 argv plus `--cofactor-biomass-terms`
  (`model/reports/cofactor_biomass_terms_20261027/logs/build_argv.txt`). Not a designated reference.
- **Reference model:** none is designated.
- **Condition:** Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`,
  no choline and no ethanolamine) + `state/strain_profiles/po1f_sd_leu.json` (`35307853…`); Gurobi
  13.0.1 at FeasibilityTol/OptimalityTol 1e-9 through
  `model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py`.

## 2. Targets as they are

| Reaction | Equation (compartment `C_em`) | Bounds | GPR | Notes / annotations |
| --- | --- | --- | --- | --- |
| R1725 "choline/ethanolaminephosphotransferase" | CDP-ethanolamine m1602 + DAG m1649 → CMP m1601 + PE m1633 + H+ m377 | (0, 0) | `YALI1C15420g` | `PROTEIN_CLASS: 2.7.8.1`; ec-code 2.7.8.2; MNXR179146; seed rxn30109 |
| R1726 "cholinephosphotransferase" | CDP-choline m1603 + DAG m1649 → CMP m1601 + PC m1632 + H+ m377 | (0, 0) | `YALI1C15420g` | `PROTEIN_CLASS: 2.7.8.2`; ec-code 2.7.8.2 |

Both are declared irreversible. Their only other connection in `C_em` is the CDP-aminoalcohol transport
(R1602, R1603) from the cytosolic Kennedy branch (kinases R212 / R293 / R2080, cytidylyltransferases R213 /
R635).

## 3. Why they are closed: provenance

**Unknown.** No rationale was found anywhere searched:
- **Models:** (0, 0) in the immutable source `model/source/iyali26.xml` (`5c8c199e…`), in
  `model/source/iyli21.xml`, in `model/candidates/legacy/model.xml`, in all five lipid-unlump candidates, and
  in the baseline. Same GPR everywhere. The SBML notes carry only `PROTEIN_CLASS`.
- **Git:** `git log -S` on R1725 / R1726 hits only commits that rebuilt candidate XMLs. The first
  `iyli21.xml` commit (`05ae6e1`, 2026-03-19) already has the closure. The repository's root commit holds an
  *E. coli* core model, so it is not a predecessor.
- **Text:** `git grep` (non-XML) and a grep of `$IYALI26_RESEARCH_ROOT/state` and `reference/locus_map`
  find R1725 / R1726 only in the 2026-10-26 census, the lipid-unlump Stage 0 formula back-solve
  (`unlump_stage0_plan.csv`, `fill_candidates.csv`: PE / PC formulas solved from these reactions) and an
  EC-name check (`state/curation/data/ec_verification.csv`: R1725's name does not match EC 2.7.8.2).
- **STATE.md:** no entry before the census.

The closure is inherited from iYli21. Whether its authors closed it on purpose is not recorded; no intent is
inferred here.

## 4. Biology and gene identity

| # | Claim | Source | Tier |
| --- | --- | --- | --- |
| B1 | `YALI1C15420g` = YALI0C10989g (W29 / CLIB122 map). | `$IYALI26_RESEARCH_ROOT/reference/locus_map/s2_metabolic_genes.csv` (row `YALI1_C15420g,YALI0C10989g,2909135,…,2.7.8.2`) | verified (file read) |
| B2 | YALI0C10989g (KEGG yli:2909135) is RefSeq "Cholinephosphotransferase 1", KO **K00994** (CHPT1/CPT1, EC 2.7.8.2), module M00090 (choline → PC). | [KEGG yli:2909135](https://www.kegg.jp/entry/yli:2909135); UniProt [Q6CCA9](https://rest.uniprot.org/uniprotkb/Q6CCA9) "diacylglycerol cholinephosphotransferase (EC 2.7.8.2)", 394 aa, annotation score 2, automatic (ARBA) annotation only | verified (opened); orthology and automatic annotation, not experiment |
| B3 | A **separate** EPT1 ortholog exists: YALI0E26565g = `YALI1E31498g` (KEGG yli:2912382, RefSeq "Choline/ethanolaminephosphotransferase 1", KO **K00993** EPT1, EC 2.7.8.1, module M00092 ethanolamine → PE). **It is not in the model** (not in the baseline, source or iyli21 XML; S2 `in_model = no`). | [KEGG yli:2912382](https://www.kegg.jp/entry/yli:2912382); UniProt [Q6C4I2](https://rest.uniprot.org/uniprotkb/Q6C4I2) (unnamed, score 1); `s2_metabolic_genes.csv` | verified (opened); orthology only |
| B4 | *S. cerevisiae* CPT1 (K00994 ortholog) has "relative specificity for CDP-choline and -dimethylethanolamine". EPT1 uses CDP-ethanolamine, -monomethylethanolamine, -dimethylethanolamine and -choline. The cpt1 ept1 double null is viable. | Hjelmstad & Bell 1991, [PMID 1847919](https://pubmed.ncbi.nlm.nih.gov/1847919/) (abstract read); UniProt [P17898](https://rest.uniprot.org/uniprotkb/P17898) lists CDP-choline and CDP-DME activities only; [P22140](https://rest.uniprot.org/uniprotkb/P22140) lists CDP-ethanolamine and CDP-choline | verified (abstract and entries read) |
| B5 | In vivo, CPT1 makes 95 % of PC via CDP-choline and EPT1 5 %. CPT1 is implicated in PC synthesis from choline released by PC turnover. | McMaster & Bell 1994, [PMID 7961735](https://pubmed.ncbi.nlm.nih.gov/7961735/) (abstract read) | verified |
| B6 | The CDP-choline route adds net PC only when choline is supplied; otherwise it recycles choline from PC turnover. | UniProt [P17898](https://rest.uniprot.org/uniprotkb/P17898) function text | verified (curated entry; *S. cerevisiae*) |
| B7 | *S. cerevisiae* psd1 psd2 double mutants are **ethanolamine auxotrophs**; < 5 % of [3H]serine reaches PE. | Trotter et al. 1995, [PMID 7890740](https://pubmed.ncbi.nlm.nih.gov/7890740/); Trotter & Voelker 1995, [PMID 7890739](https://pubmed.ncbi.nlm.nih.gov/7890739/) (abstracts read) | verified (*S. cerevisiae*) |
| B8 | *Y. lipolytica* phospholipid synthesis "is still poorly understood". No Y. lipolytica CPT1 / EPT1 deletion or enzyme assay was found. The SPO14 (PLD) gene YALI0E18898g (= `YALI1E22633g`, model PLD R1743 / R1842, on the consensus essential list) was deleted in W29-derived strain PS05, giving the viable strain PS07. | Szczepańska et al. 2022, [PMID 36142650](https://pubmed.ncbi.nlm.nih.gov/36142650/), full text [PMC9505966](https://pmc.ncbi.nlm.nih.gov/articles/PMC9505966/) (passage read) | verified for what the paper states; the absence of CPT1/EPT1 work is from one standard web search, not exhaustive |

**What follows (inferred):**
- R1726's GPR (CPT1 ortholog → CDP-choline + DAG → PC) is consistent with B2, B4 and B5.
- R1725's GPR is not. By B4 the yeast CPT1 homologue is not an ethanolaminephosphotransferase. The
  reaction's own name ("choline/ethanolaminephosphotransferase") and its `PROTEIN_CLASS 2.7.8.1` describe
  EPT1, whose ortholog (B3) is absent from the model.
- None of this is *Y. lipolytica* experiment. B2 and B3 are orthology and automatic annotation; B4–B7 are
  *S. cerevisiae*.

## 5. In-memory LP results (Task A + PO1f, Gurobi 1e-9, checked; about 170 LPs)

**Growth and FVA.** Every WT point is optimal and passes the check.

| Variant | Growth | FVA R1725 at 100 % / 99.9 % / ≥ 0 | FVA R1726 at 100 % / 99.9 % / ≥ 0 |
| --- | ---: | --- | --- |
| closed (as built) | 0.3120049202 | 0 / 0 / 0 | 0 / 0 / 0 |
| R1725 open | 0.3120197782 (+1.5e-5) | 1.819e-4 fixed / [0, 1.146e-3] / [0, 0.797] | — |
| R1726 open | 0.3120558110 (+5.1e-5, +0.016 %) | — | 6.147e-4 fixed / [0, 0.02235] / [0, 21.33] |
| both open | 0.3120558110 | [0, 0] / [0, 1.048e-3] / [0, 0.797] | 6.147e-4 fixed / [0, 0.02235] / [0, 21.33] |

**Corrected energy check.** Every exchange is set to (0, 1000) and every positive lower bound is relaxed to
0; the only one is `xMAINTENANCE` 7.8625. R_NTP1 is then maximised. The LP is **optimal at 0** in all four
variants, so there is no energy-generating cycle. This agrees with the census re-check on a1237821
(REPORT section 6).

**What the opened routes do at the optimum** (pFBA difference vs closed; `stage1_pfba_diff.json`, artifacts folder). pFBA
optima are not unique: hundreds of degenerate swaps appear, such as MDH ↔ ADH, so only the lipid part is read.
- **R1726 runs a choline-neutral PC turnover cycle.** CDP-choline + DAG → PC (R1726) → PC to `C_en` (R1806)
  → PLD (R1842) → PA + choline → choline back to the cytosol (R1042) → choline kinase R212 → CCT R213 →
  CDP-choline (R1603), all at 6.147e-4. Net effect: DAG + ATP + CTP → PA. It replaces the lipid-particle
  DAG-lipase / MAG-lipase / acyl-CoA re-ligation detour (R1763, R1765, R1456, xPOOL_AC_LP / FA_LP → 0).
  That detour is used because the **DAG kinase R1720 (`YALI1F25444g`) is itself closed (0, 0)** in the
  model (one of the 69 inherited closures; census UM80). There is no net PC synthesis, because no choline
  enters on Task A. The cycle has a biological counterpart (B5, B6), but the +0.016 % growth gain is an LP
  side effect of the closed DGK, **not evidence** for the change. With PLD closed, R1726 still raises growth
  slightly (0.3120258). The route was not traced; phospholipase B / GPC turnover is a guess (inferred).
- **R1725 makes PE from sphingoid-base lyase phosphoethanolamine.** Serine palmitoyltransferase R722 →
  R37 → LCB kinase R734 → lyase SPHPL → P-ethanolamine → ECT R635 → CDP-ethanolamine → R1725, all at
  1.819e-4. This replaces the whole PS synthase → PSD1 route (R1841, R629 → ~0).

**Targeted controls** ([stage1_controls.json](stage1_controls.json)):

| Control | closed | R1725 open | R1726 open | both open |
| --- | ---: | ---: | ---: | ---: |
| PSD1 + PSD2 double KO (`YALI1D27274g`, `YALI1D04432g`) | 0 | **0.3120198 (WT)** | 0 | **0.3120556** |
| same, with lyase reactions SPHPL / R730 / R663 closed | 0 | 0 | 0 | 0 |
| PLD reactions R1743 + R1842 closed | WT | WT | 0.3120258 | 0.3120345 |
| PLD gene `YALI1E22633g` KO | WT | WT | 0.3120258 | 0.3120345 |

With R1725 open, the model predicts that a PSD-null strain grows at wild-type rate without ethanolamine,
because the sphingoid-base lyase supplies all PE. The verified *S. cerevisiae* phenotype is ethanolamine
auxotrophy (B7). The real limit is lyase capacity, which a stoichiometric model does not carry. This is the
main argument against opening R1725 as it stands.

**Targeted knockouts** (18 route genes: Kennedy, PSD, PS / PI / CDP-DAG synthase, PE methylation, PLD,
PLB, lyase, R2187). There is **no class change at any cutoff**, closed vs both open:
- **Still essential (ratio 0):** `YALI1D10930g` PS synthase, `YALI1E07347g` PE methyltransferase,
  `YALI1E15432g` phospholipid methyltransferase, `YALI1E17431g` CDS, `YALI1F27002g` PI synthase. All five are
  on the consensus list.
- **Ratio ≥ 0.9998:** every other gene, including `YALI1C15420g` itself and the PLD `YALI1E22633g`
  (consensus-listed, FN before and after: this change does not touch that FN).

Choline / ethanolamine are absent from Task A, so PE methylation stays the only net PC route, as B6 predicts.

**Mechanism control** (not a benchmark condition): choline R1041 and ethanolamine R1115 opened to uptake 1.
With both reactions open, growth is 0.3120834 and PE methylation becomes dispensable (`YALI1E07347g`,
`YALI1E15432g` → 1.0). Supplied choline makes net PC through R1726, as the *S. cerevisiae* biology describes.
With the closures as built, choline cannot be used (both genes stay at 0).

## 6. Options

| | Change | For | Against |
| --- | --- | --- | --- |
| **A (recommended)** | R1726 (0, 0) → (0, 1000); GPR unchanged; R1725 stays (0, 0), with the reason recorded in the curation file and checked by the guard | GPR matches the gene (B2, B4, B5); unblocks UM12 (census fix test: R1726 alone 11/11); choline becomes usable; no EGC; no knockout class change in the targeted set | Inherited closure of unknown intent; growth +0.016 % comes from a turnover cycle that stands in for the closed DGK R1720 (a separate gap); UM19 (8) stays blocked |
| B | A + R1725 opened with its GPR changed to the EPT1 ortholog `YALI1E31498g` (new model gene), optionally EPT1 as an OR member of R1726 | Unblocks UM19 too; gene matches the reaction by orthology | Two changes, one adding a gene with orthology-only identity; the PSD-null rescue artefact remains (needs a lyase-capacity decision first); AND/OR evidence rule (AGENTS §5) only partly met (B4 in vitro, B5 5 % in vivo) |
| C | Open both as they are | Census 19/19 | R1725's GPR contradicted by B4; PSD-null rescue artefact |
| D | Keep both closed; record the findings only | No change | Leaves a pathway that the yeast biology and the gene annotation support, and leaves supplied choline unusable |

## 7. Proposed change (Option A) and implementation

- **Data:** `model/curation/direction/r1726_cpt_open.json` (unique name; schema as `r2041_direction.json`).
  - `curation_id` R1726-CPT-OPEN-<date>, the authorization quote, the reviewed model `00e8e397…`.
  - R1726 `before_bounds` [0, 0], `after_bounds` [0, 1000]; its exact stoichiometry, species formulas,
    charges and compartments, and GPR `YALI1C15420g`.
  - A `kept_closed` record for R1725 with the reasons B3 / B4 / B7 and the PSD-null control.
  - The evidence rows of section 4, each with URL, tier, scope and access date; the reviewer.
- **Step:** `platform/scripts/gem_annotate/r1726_cpt_open.py`, `apply_r1726_cpt_open(model)`.
  - It fails closed if R1726 is missing, its stoichiometry, species, GPR or reversibility differ, its bounds
    are neither [0, 0] nor [0, 1000], or R1725 is not (0, 0).
  - It returns `applied` / `already_correct` and touches only R1726's bounds. Notes follow the
    `r2041_direction` precedent; any note change is reported in the diff as a deliberate change.
- **Flag:** `--r1726-cpt-open` in `cli.py` / `main.py`. It runs after `--cofactor-biomass-terms` (the last
  opt-in step) and does not require it (a bound change does not depend on `biomass_C`).
- **Test:** `platform/tests/test_r1726_cpt_open.py`, on an in-memory model.
  - It applies on the 00e8e397 state and is idempotent.
  - It refuses changed stoichiometry, a changed GPR, unexpected bounds, and an opened R1725.
  - Nothing else in the model changes.
- **Build:**
  `python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" "${BASE[@]}" --cofactor-biomass-terms --r1726-cpt-open --output-model ../model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r1726_cpt_<date>.xml`
  (`BASE` = the a1237821 argv as a zsh array).

## 8. Expected effect (inferred from Stage 1 unless marked)

- **Whole-model diff:** R1726 bounds only, plus notes if the step adds them.
- **WT, Task A:** 0.3120049 → 0.3120558 (observed in memory).
- **Blocked reactions, open condition:** UM12, 11 fewer (census fix test on a1237821; the re-census gives the
  exact number on 00e8e397). UM19 stays blocked. No newly blocked reaction expected.
- **Task A:** R1726, R212, R213, R1603, R1042, R1806, R1842, R1776 become active (observed in pFBA).
- **Knockouts:** no class change in the 18 route genes (observed). The full sweep may show small ratio shifts
  in lipid-particle lipase / acyl-CoA ligase genes, because the cycle replaces their detour; a class change
  there would be a stop-and-explain item.

## 9. Acceptance criteria

1. Guarded step, test and flag as in section 7. The full pytest suite is equal to the previous run test by
   test, plus the new file passing. Code-reviewer findings are fixed or reported.
2. A control build without the flag reproduces `00e8e397…` byte for byte. Seeds 0 and 1 with the flag give
   identical bytes.
3. The whole-model diff is limited to R1726 (bounds, and notes if any).
4. WT on Task A + PO1f is optimal and passes the check. The relaxed no-uptake ATP maximum is 0 (optimal) on
   the built file.
5. Two-model knockout sweep (Gurobi 1e-9, reset per knockout, knocked-flux check): every class change is
   traced to the R1726 route. New false positives are reported, never fixed. Scores at 1 / 5 / 10 / 15 %,
   with the user-defined negative class, mtDNA out of scope, and a development-reference label.
6. Re-census (open condition and Task A, census scripts, GLPK): the freed set is UM12 or explained, and no
   reaction is newly blocked unless explained by the solver fringe and checked at Gurobi 1e-9.
7. Independent read-only audit verdicts recorded; disagreements reported.
8. REPORT.md and a STATE.md entry, including the R1725 rationale and the follow-ups in section 11.

## 10. Budget and stop conditions

- **Budget:**
  - phase-1 re-check on the built file, about 20 LPs;
  - 3 builds;
  - 2 × 1,082 knockout LPs;
  - 4 census FVA scans (2 models × 2 conditions);
  - the pytest suite twice;
  - one code review and one audit.
- **Stop and report on:**
  - a guard rejection;
  - a diff outside R1726;
  - seed-dependent bytes;
  - a non-optimal WT;
  - a non-zero or non-optimal relaxed EGC;
  - a test regression;
  - a knockout class change not traceable to the R1726 route;
  - a newly blocked reaction that Gurobi 1e-9 confirms.

  Nothing is tuned or relaxed to pass.

## 11. Not in scope (findings and follow-ups only)

- R1725 / UM19 and the EPT1 ortholog `YALI1E31498g` (Option B). This needs a decision on sphingoid-base
  lyase capacity first: SPHPL is reversible (−1000, 1000) and supplies unbounded phosphoethanolamine.
- The closed DAG kinase R1720 (`YALI1F25444g`, census UM80). It is the reason the R1726 cycle raises growth.
- R2187 (GPR `YALI0E12442g`, a YALI0 id in a YALI1 model, 0.01 stoichiometry, product cytosolic PC):
  anomalous, not examined.
- Compartment `C_em` is named "endosomal membrane" in the model, but its reactions are named "ER
  membrane".
- R1725's ec-code 2.7.8.2 conflicts with its `PROTEIN_CLASS` 2.7.8.1.
- **PLD `YALI1E22633g`:** consensus-essential, but deleted in a viable W29-derived strain (B8). This is an
  experimental-list question; no label change.
- **Lipid-unlump line (colleague-owned, review only):** its candidates keep R1725 / R1726 closed, and its
  Stage 0 back-solved PE / PC formulas from them. If the lines merge, the unlumped Kennedy reactions need the
  same decision.

## 12. Authorization sentence (265 characters)

> Approved: Option A of the Kennedy TASK draft: open R1726 to (0,1000) via a guarded --r1726-cpt-open step after --cofactor-biomass-terms on 00e8e397; R1725 stays closed; test, builds, diff, Gurobi KO sweep, re-census, review, audit, REPORT, STATE. No commit or push.

For Option B or D, say so and the scope is redrafted. Option B adds a gene and needs the lyase decision first.

## 13. Re-run on e9d35857 (authorized 2026-10-27)

**Authorization.** The user, 2026-10-27, after reading the re-run draft (shared notes folder,
`kennedy_r1726_rerun_TASK_draft_20261027.md`; its section 6 gives the full approval sentence):

> Approve re-run

The stop of the first run is kept in [STOP.md](STOP.md). The re-run's outputs are in [rerun/](rerun/). The
draft follows unchanged, its headings demoted one level.



**Status (when drafted): draft for the user's decision. Approved as "Approve re-run"; carried out, see REPORT.md.**

- **Pre-check:** the only computation for this draft is
  `model/reports/kennedy_closures_20261027/scripts/rerun_precheck.py` → `rerun_precheck.json` (25 LPs,
  read-only).
- **Report location:** on approval this becomes section 13 of `model/reports/kennedy_closures_20261027/TASK.md`.
  That folder keeps its STOP.md as history and gets the REPORT.md and AUDIT.md of the re-run.

### 1. What changes since the stop

- **The leak is fixed.** The stop was caused by the R2135 acetyl-CoA → acyl-CoA mapping. It is fixed in
  candidate `e9d358574eeb4369ecd39f63e5a2ad6c113c4681232c6504cb1ef42c4fe89b1b` (= `00e8e397` +
  `--r2135-acyl-coa-pool`; REPORT and AUDIT in `model/reports/r2135_acyl_coa_pool_20261027/`; not committed).
- **The scientific decision is unchanged:** R1726 goes from (0, 0) to (0, 1000) with GPR `YALI1C15420g`, and
  R1725 stays closed, for the reasons of the first draft (CPT1-ortholog GPR, EPT1 ortholog absent, PSD-null
  rescue).
- **The step, data and test already exist** (uncommitted). Three things change:
  1. **Requirement:** `--r1726-cpt-open` requires `--r2135-acyl-coa-pool` (the R2135 reviewer's
     out-of-scope finding). Without it, opening R1726 brings the acyl leak back. The R1726 test gets the
     matching requirement case.
  2. **Curation file** `model/curation/direction/r1726_cpt_open.json`:
     - `reviewed_model` becomes `e9d35857…`;
     - `authorization` records both approvals;
     - the note `r1726_cpt_open_limits` is rewritten. Its clause "on the Task A medium its growth gain
       stands in for the closed DAG kinase R1720" is no longer true: on `e9d35857` R1726 adds no growth
       on Task A.
     - Notes are model text, so this lands in the diff as a deliberate change.
     - The guard is unchanged: it pins R1726 and R1725 only.
  3. **Build base:** the a1237821 argv + `--cofactor-biomass-terms --r2135-acyl-coa-pool
     --r1726-cpt-open`.

### 2. Evidence on the new base (verified unless marked)

| Check | Result | Source |
| --- | --- | --- |
| WT, Task A + PO1f | 0.3116824 with or without R1726 (no growth gain) | `r2135_acyl_coa_pool_20261027/built_checks.json` |
| R1726 at 99.9 % of the optimum | FVA [0, 0.0217]; 0 at the pFBA point | same |
| No-uptake ATP and biomass with R1726 open | 0 / 0 (optimal, maintenance relaxed) | same |
| Knockouts with R1726 open (1,082, checked) | **no class change** versus `e9d35857` at any cutoff; FAS / OLE1 / FAA genes lethal | `r2135_acyl_coa_pool_20261027/ko_sweep_r1726_stacked_summary.json` |
| Open condition (every exchange ±1000), UM12 reactions | all 11 [0, 0] as built → unblocked with R1726 open (R1726, R212, R213, R1603, R830, R1041, R1042, R1743, R1842, R1809, R354) | `rerun_precheck.json` |
| Choline + ethanolamine uptake 1, PE-methylation KOs (`YALI1E07347g`, `YALI1E15432g`) | lethal as built → 0.3117034 (= WT) with R1726 open: supplied choline makes PC, as the *S. cerevisiae* biology describes (PMID 7961735, P17898) | `rerun_precheck.json` |

**What the change is now worth:**
- **Network completeness:** UM12, 11 reactions, unblocked.
- **Use of supplied choline.**
- **No Task A growth or essentiality effect** at the 18 targeted genes or in the full stacked sweep.

The +5e-5 growth of the first run was the leak, not the pathway.

### 3. Implementation and acceptance

As in the first draft (sections 7–9), with these differences:
- **Control build:** the base argv without `--r1726-cpt-open` must reproduce `e9d35857…` byte for byte. Seeds
  0 and 1 with it must be identical, and `PYTHONHASHSEED` is recorded in the logs this time (an audit gap of
  the R2135 task).
- **Whole-model diff vs `e9d35857`:** R1726 bounds and notes only.
- **Built-file checks:**
  - WT;
  - no-uptake ATP and biomass;
  - FVA of R1726 and R1725;
  - the PSD1 + PSD2 double knockout (stays lethal);
  - the choline / ethanolamine control;
  - the four FAS-group knockouts (stay lethal).
- **Two-model knockout sweep** `e9d35857` vs the new candidate: no class change expected. Any change is
  stop-and-explain.
- **Re-census (GLPK):**
  - open condition: UM12's 11 reactions freed (pre-checked above);
  - Task A: the R1726 route may become active (FVA range above), with the exact number from the census;
  - no Gurobi-confirmed newly blocked reaction.
- **pytest, review, audit, REPORT, STATE**, as before.

### 4. Budget and stop conditions

- **Budget:**
  - 3 builds;
  - about 40 built-file LPs;
  - 2 × 1,082 knockout LPs;
  - 2 census FVA scans; the `e9d35857` census is reused;
  - pytest twice;
  - one review and one audit.
- **Stop and report on:**
  - a guard rejection;
  - a diff outside R1726;
  - seed-dependent bytes;
  - a non-optimal WT;
  - non-zero no-uptake ATP or biomass;
  - a test regression;
  - any knockout class change;
  - a Gurobi-confirmed newly blocked reaction.

### 5. Not in scope

- R1725 / EPT1;
- the DAG kinase R1720;
- the other 0.01-scaled reactions;
- other media.

### 6. Authorization sentence (267 characters)

> Approved: re-run Kennedy Option A on e9d35857 as in the re-run draft: --r1726-cpt-open requires --r2135-acyl-coa-pool, R1726 notes and reviewed model updated, R1725 stays closed; test, builds, diff, Gurobi KO sweep, re-census, review, audit, REPORT, STATE. No commit.
