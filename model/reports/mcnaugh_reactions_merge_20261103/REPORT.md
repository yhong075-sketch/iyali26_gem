# REPORT: 17 reviewed mcnaugh_pass reactions, lumped Kennedy steps closed (2026-11-03)

Scope, authorization, baseline and acceptance criteria are in [TASK.md](TASK.md). The audits are in [AUDIT.md](AUDIT.md).
The large outputs are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/mcnaugh_reactions_review_20261103/` (below:
`artifacts/`). The look-alike review is in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/mcnaugh_duplicate_biology_20261103/`.

## 1. What changed and why

Andrew McNaughton's lipid-unlump handoff (`mcnaugh_pass`) splits our lumped lipid-particle Kennedy pathway into
chain-specific steps and adds two biomass reactions. The user wants these merged so that the model follows the real
biology. Andrew owns the lipid-unlump line and approves.

- **Builder step:** `apply_mcnaugh_reactions` (`platform/scripts/gem_annotate/mcnaugh_reactions.py`), opt-in flag
  `--mcnaugh-reactions`. It runs before `--bigg-crosswalk-annotations`.
- **Curated data:** `model/curation/mcnaugh_pass/mcnaugh_reactions_20261103.json`, revision 2 (`a62050e3…`). It is
  written by `artifacts/scripts/make_mcnaugh_reactions.py` from our v2 model and Andrew's `mcnaugh_pass/iyali26.xml`
  (`e74a0f50…`). It holds the 17 reactions, the 4 closures, the 12 not added, the held COQ2 and the known artefact.
- **Added:** Andrew's equations translated term for term into our existing species, with his bounds and GPRs. No new
  species or genes.
  - G3PAT160/161/180/181/182_lp: glycerol-3-phosphate + acyl-CoA (C16:0, C16:1, C18:0, C18:1, C18:2) → 1-acyl-G3P
    (m577) + CoA;
  - AGPAT160–182_lp: 1-acyl-G3P + acyl-CoA → PA (m1705) + CoA;
  - DGAT160–182_lp: DAG (m1642) + acyl-CoA → TAG (m1640) + CoA;
  - BIOMASS_C_core and BIOMASS_C_full, objective 0.

  The 15 lipid steps join our group "Glycerolipid metabolism", which also holds the lumped R353 and R1771 they
  replace. Andrew's reactions carry no subsystem, and R1846 is in no group. The biomass reactions join no group, like
  biomass_C.
- **Closed** (bounds 0..0; the reactions stay):
  - R353, R1846, R1771: the lumped lipid-particle steps with the generic acyl-CoA m570;
  - R2165: cytosolic lumped TAG synthesis, removed in Andrew's design. It carries no flux in v2.
- **Candidate (v4)**, model SHA-256 `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2`. Two files
  hold it, byte-identical:
  - `model/candidates/E5_…_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_20261103.xml`, **the file of record**. Its
    `.build.json` records the final step code (`3e1b593b…`). Command in
    [logs/build_command_final.txt](logs/build_command_final.txt).
  - `…_mcnaugh_rxns_v2_20261103.xml`, built before the duplicate-guard fix with step code `b9473e58…`. Command in
    [logs/build_command_v4.txt](logs/build_command_v4.txt).

  Both are v2 (`9158e9f2…`, commit `048e18ec`) plus the flag, built from the uncommitted tree on `048e18ec`.
- **Superseded:** candidate v3 `…_mcnaugh_rxns_20261103.xml` (`ec7e3465…`). It has the same LP as v4, but the 15
  lipid steps are in no group and the biomass notes are short (audit F1, F3). Its inputs are kept in
  `artifacts/superseded_ec7e3465/`, and the v4 guard code in `artifacts/superseded_guard_b9473e58/`.

**Genes on the new reactions.** The identities below come from the reviewers, who opened the sources listed. I did not
re-open them for this report, so they are second-hand here. YALI0 IDs are from `w29_tag_protein_audit.csv`, where
every one is flagged `ok`.

| Reaction set | Gene | Name | Function | Evidence |
| --- | --- | --- | --- | --- |
| G3PAT | YALI1C00230g = YALI0C00209g | SCT1 (Tsakraklides 2018, [PMC5941336](https://pmc.ncbi.nlm.nih.gov/articles/PMC5941336/)) | glycerol-3-phosphate / DHAP acyltransferase | Y. lipolytica experimental |
| AGPAT | YALI1E22736g = YALI0E18964g | SLC1, supported only through KEGG K13509; none verified | 1-acyl-G3P acyltransferase | database/homology |
| AGPAT | YALI1F25951g = YALI0F19514g | none verified | ALE1-family (K13519) acyltransferase | database/homology; lipid-particle location unverified |
| DGAT | YALI1E38810g = YALI0E32769g | DGA1 (Beopoulos 2012, [PMC3275733](https://pmc.ncbi.nlm.nih.gov/articles/PMC3275733)) | DGAT2-family DAG acyltransferase, lipid bodies | Y. lipolytica experimental |
| DGAT | YALI1D10264g = YALI0D07986g | DGA2 (Beopoulos 2012) | DGAT1-family DAG acyltransferase; ER-like location (Gajdos 2016) | Y. lipolytica experimental |

A GPR in the model is an assignment, not experimental validation.

## 2. How the scope was set (evidence)

- **Andrew's 56 reactions** with no recorded counterpart in our model were reviewed by a multi-agent workflow: a
  reviewer, a skeptic told to refute, and a tie-breaker where they disagreed. Each reaction was simulated on Task A +
  PO1f with Gurobi 1e-9.
  - **26 look-alikes:** all are duplicates of existing reactions. None was added. The side findings (RNR GPR,
    folate shortcuts, annotation swaps and more) went to separate sessions.
  - **29 remaining** (`artifacts/review_output.json`, `review_summary.txt`):
    - 17 added;
    - 10 are duplicates once Andrew's species are mapped to ours. Re-checked after the audit:
      - 7 are term for term: 4ABUTtm = R975, ABUTt2r = R809, PAC2HYD = R2191 and DIPEP1–4_v =
        R2021/R2029/R2034/R2039.
      - 3 match up to one cytosolic H+ (m10), the protonation convention: 3HPCOAHYD = R2140, COMT2 = R2097 and
        COMT3 = R2098.
      - 4ABUTtm and ABUTt2r match only with Andrew's 4abut_c mapped to m133. His crosswalk maps it to m2040, a
        second cytosolic GABA species already in our model (m133 MNXM1144700, m2040 MNXM192). That duplicate
        species is pre-existing and belongs to another line.
    - 2 would reverse an authorized curation: ALCD2irm re-adds R132 (removed by R132-MITO-ADH-20261013; +37.9 % WT
      growth), and G3PD1irm re-adds R349 (removed by CoQ9 respiratory package item C4).
  - **COQ2 (cytosolic):** held by the CoQ9 gate.
- **R1771 must close in the same change.** With the new DGATs present and R1771 still reversible, the pair forms a
  loop and WT growth reaches 0.5127 (`artifacts/inputs/facts.json`, `together.growth_all_added`; observed).
- **R351 (DHAP acyltransferase, SCT1) stays open.** Andrew's design has no counterpart for it. Keeping it departs
  from his design, and it is on his list.

**User decisions that differ from the reviewers' recommendations:**
- **DGAT GPR:** the reviewers recommended DGA1 alone, because only DGA1 is localised to lipid bodies and DGA2 already
  carries the ER DGAT R2168. The user chose Andrew's "DGA1 or DGA2".
- **Acyl pools:** the skeptics called making every xPOOL_AC_*/xPOOL_FA_* forward-only "required". The user chose
  "Record it, change later" (§5).
- **BIOMASS_C_core:** the reviewer proposed our biomass_C without the TAG term. The user chose Andrew's reaction
  verbatim ("Both verbatim from Andrew").

## 3. Code logic

`apply_mcnaugh_reactions(model, spec=None)` works as follows:
1. **Scope check:** it requires schema 1, change ID `MCNAUGH-REACTIONS-20261103`, status `user_authorized`, exactly
   17 reactions and 4 closures, and `biomass_C` as the only objective.
2. **Group check:** every group a row names must exist exactly once.
3. **Pre-validation, before any edit:**
   - every species and GPR gene of each new reaction must already exist;
   - **no duplicates:** no reaction under another ID may have the row's stoichiometry, either as written or reversed.
     The map from stoichiometry to reaction IDs keeps every ID, so the result does not depend on reaction order;
   - a new reaction ID that already exists must match the curated row in name, stoichiometry, bounds, GPR,
     annotation, notes and group membership, so it counts as already applied;
   - each closure must keep its curated stoichiometry and GPR. Stoichiometry is compared, not the equation string,
     because the string's arrow changes with the bounds this step sets;
   - each closure's bounds must be the before or the after bounds.
4. **State check:** the per-row states must be all `before` or all `after`. A mix, or any stale bound, raises without
   edits.
5. **Edit:** from `before`, it adds the 17 reactions, then for each one sets:
   - its GPR;
   - its annotation (SBO:0000176, or SBO:0000629 for biomass);
   - its notes;
   - its group membership. The subsystem is the group's name, because cobra writes a subsystem only through groups.

   It then sets the 4 closures to 0..0. From `after`, it returns `already_correct`.

Note texts avoid characters that SBML escapes, so the notes survive a round trip unchanged.

## 4. Results (observed this session)

- **Whole-model diff, v3 against v2** ([verify_v3.json](verify_v3.json)): only the 17 added reactions and the bounds
  of the 4 closures differ. No other species, reaction, gene, group, compartment, objective, annotation or note
  changed.
- **Whole-model diff, v4 against v3** ([diff_v4_v3.json](diff_v4_v3.json)): 18 expected differences and none
  unexpected:
  - the 15 lipid steps gain subsystem "Glycerolipid metabolism";
  - the 2 biomass reactions get the new `composition_note`;
  - group `group33` gains the 15 lipid steps.

  Stoichiometry, bounds, objective, species, GPRs and genes are identical, so it is the same LP. The re-reviewer
  confirmed this with a raw text diff. The third build reproduced the v4 file byte for byte.
- **Balance of the new reactions:**
  - The 15 lipid steps are mass-unbalanced. This is the recorded exception: the generic species m577, m1705, m1642 and
    m1640 carry no balanced chain formula. In Andrew's species the imbalance is exactly the R group. In ours the
    numbers in [verify_v3.json](verify_v3.json) come from our species' formulas (audit F4):
    - AGPAT shows only the acyl group (for example −C16H30O), because m577 and m1705 both have no formula;
    - G3PAT shows the whole product missing (for example −C19H39O7P), because m577 has no formula;
    - DGAT shows O +56, a parse artefact: m1640's formula `C6H5O6*3` reads as O63.

    The `*` formulas (m1640, m1123, m1648) are pre-existing and belong to the formula/charges line.
  - The three C18:2 steps are also +4 in charge. This comes from m1735[C_lp] (linoleoyl-CoA, charge −4 on a neutral
    formula), a pre-existing species defect named in their notes. Andrew's `linoleoyl_neutral_correction` is recorded
    in his curation and is not applied.
  - The biomass balance cannot be computed: some species have no formula.
- **FVA (fraction 0, on v3; same LP in v4):** every new reaction can carry flux. The upper bounds range from 3.59
  (G3PAT180_lp, AGPAT180_lp) to 40.78 (DGAT160_lp).
- **Growth (Task A + PO1f, Gurobi 1e-9, checked solves):**

  | Objective | v2 | v3 | v4 |
  | --- | ---: | ---: | ---: |
  | biomass_C (the objective) | 0.31205354 | 0.31248856 (+0.139 %) | 0.31248856100769146, identical to v3 |
  | BIOMASS_C_core | n/a | 0.31335790 | same LP as v3 |
  | BIOMASS_C_full | n/a | 0.31259249 | same LP as v3 |

- **Energy check:** with every uptake closed (35), maximum growth is 0.0 on v3 and on v4
  ([wt_energy.json](wt_energy.json), [wt_energy_v4.json](wt_energy_v4.json)).
- **Certified knockout screen, biomass_C** (run on v3, then re-run on v4 on 2026-11-04, "do a new screen test";
  [ko_screen_biomass_C_summary.json](ko_screen_biomass_C_summary.json),
  [ko_compare_v3_vs_v2.json](ko_compare_v3_vs_v2.json); v4 run in the artifacts `screen_v4_20261104/`):
  - **v4 reproduces v3 (observed):**
    - the v4 TSV is byte-identical to the v3 TSV;
    - the summary differs only in the model SHA;
    - the log differs only in the model prefix and timings.
    - An independent re-score from the raw consensus CSV matches all 16 confusion matrices (rules A and B, two scopes,
      four cutoffs).
  - 1,082 sweep LPs, 1 ceiling and 40 certificate LPs per run;
  - no sweep point fails the check, and all 20 non-optimal knockouts are certified;
  - 2 genes change class, at every cutoff and under both rules:
    - YALI1C07638g = YALI0C05951g: consensus essential, now non-essential (TP → FN). It is the Δ9 desaturase; the
      name OLE1 is recalled, not verified this round. It catalyses R1434, R1435 and R1869. KO growth is 0.31248856,
      equal to WT.
    - YALI1D22124g = YALI0D17864g: consensus non-essential, model-essential in v2, non-essential in v3 (FP → TN). It
      is a fatty-acid–CoA ligase; the name FAA1 is recalled, not verified.
      - Its KO growth, 0.3120960274, equals v3 WT with forward-only pools (observed by this session and by the
        auditor).
      - So it probably blocks the same chain-interconversion route behind the +0.139 % artefact. Its flip is then the
        same composition escape as OLE1 (inferred; the route was not traced).
      - The net score change is two effects of one design, not an independent fix.
  - The genes on the new and closed reactions do not change class:
    - SCT1 stays model-essential, a consensus non-essential (FP);
    - SLC1 stays non-essential, a consensus essential (FN), through its OR with YALI1F25951g;
    - DGA1, DGA2, TGL3 and TGL4 stay non-essential.

  **Scores.** Rule B, mtDNA and untargeted genes out of scope. Model genes not on the consensus list count as
  non-essential, a user-defined negative class. This is a development reference, not independent validation.

  | Cutoff | v2 TP/FN/FP/TN | v2 MCC | v3/v4 TP/FN/FP/TN | v3/v4 MCC |
  | --- | --- | ---: | --- | ---: |
  | 1 % | 114/212/37/701 | 0.3957 | 113/213/36/702 | 0.3957 |
  | 5 % | 121/205/39/699 | 0.4105 | 120/206/38/700 | 0.4105 |
  | 10 % | 125/201/44/694 | 0.4084 | 124/202/43/695 | 0.4082 |
  | **15 %** | **142/184/44/694** | **0.4563** | **141/185/43/695** | **0.4562** |

  Positive coverage: 326 consensus positives are model genes. Rule A and the other scopes are in the comparison file.
- **Screen with BIOMASS_C_core as the objective** (for information, on v3; same model, medium and threshold;
  [ko_screen_core_summary.json](ko_screen_core_summary.json),
  [ko_compare_core_vs_biomass_C.json](ko_compare_core_vs_biomass_C.json)):
  - 15 %, rule B, untargeted out: 97/229/29/709, MCC 0.3685.
  - 58 genes are essential with biomass_C and non-essential with the core: 44 consensus positives and 14 negatives.
  - The cause is composition, not mechanism. The core reaction lacks 18 of biomass_C's terms (§5), so the genes that
    make them become non-essential.

## 5. Known limits and what they mean

- **OLE1 false negative (accepted by the user, "Accept, finish").** In v4, each chain-specific step makes the same
  generic product (m577, m1705, m1640). The optimum can therefore build every lipid from saturated chains, and the
  desaturase is no longer needed.
  - **Pools ruled out (observed):** making the acyl pools forward-only does not restore essentiality (KO = WT =
    0.312096).
  - **Generic-product design (inferred):** this cause is read from the stoichiometry. No positive test was run, for
    example one that re-imposes the chain composition.
  - In v2, the lumped steps drew on the generic acyl-CoA m570[C_lp]. It is made only by xPOOL_AC_LP, from a fixed
    chain mix whose largest shares are C18:1 0.472, C18:2 0.182, C16:0 0.174 and C16:1 0.079 (statically verified
    in v3). In v4 only R351 and R1766 still use m570.
  - Side observation, a candidate that was not checked: R1869 (oleoyl-CoA → linoleoyl-CoA, a Δ12 step) carries
    YALI1C07638g, the Δ9 desaturase. In Y. lipolytica the Δ12 step is usually assigned to FAD2 (recalled). This needs
    its own GPR review.
  - The biology says Ole1 is essential without unsaturated fatty acids in the medium. This is recalled textbook
    knowledge, uncited here, and the consensus list calls the gene essential. v4 therefore loses a correct call.
  - Restoring it needs a chain-composition constraint, chain-specific lipid species, or a pool-drawn composition.
    That is a design question for Andrew.
- **+0.139 % WT growth (recorded artefact, "Record it, change later").** Chain-specific acylation into generic
  species, combined with our reversible xPOOL_AC/FA pools, lets chains interconvert. With forward-only pools the gain
  would be +0.014 % (observed). The gain is under the 0.5 % gate.
- **R-group mass imbalance** on the 15 lipid steps: copied from Andrew's design and user-accepted.
- **The +4 charge on the C18:2 steps** stays until m1735[C_lp] is fixed by the charges task or by Andrew's correction.
- **BIOMASS_C_core and BIOMASS_C_full are verbatim from Andrew and not for screens.** Their notes list the
  differences from biomass_C, term by term:
  - Both lack 17 of its terms:
    - 13 cofactors: NAD, NADP, FAD, riboflavin, CoA, THF, TDP, ubiquinone-9, heme a, iron, lipoylprotein, lipoyl-E2,
      dolichyl phosphate;
    - the ceramide pool;
    - the m2062/m2064 apo-protein by-products;
    - the m1401 biomass product.
  - The core reaction also lacks the TAG term, so it lacks 18.
  - Both draw 20 charged tRNAs directly, in place of our `trna_biomass_residue_*` terms.
  - The full reaction keeps trehalose, which was removed from biomass_C on 2026-10-30. As the objective it would make
    TPS1 lethal again (reviewer run).
  - Both are drains with objective 0. They can carry flux in FVA or sampling under biomass_C.
- **libSBML notes error on R304:** every build since the v2 reproduction logs "Setting notes on sbase R_R304". It is
  pre-existing, unrelated to this change, and was not investigated.

## 6. Could the method be better, and how to test it

- **Composition escape** is the main design weakness of the generic-product Kennedy split. Two candidate fixes:
  - chain-specific lipid species with a fixed chain-composition drain;
  - one constraint per lipid class that ties chain fluxes to the measured acyl profile.

  Either can be tested in memory. Check that OLE1 returns to essential, that WT growth moves back toward 0.31205, and
  that the energy check stays 0.
- **Pool reversibility:** make xPOOL_AC/FA forward-only as a separate, authorized change. Test it with WT growth,
  the energy check and a full screen. Here it changed WT growth to 0.312096.
- **Guard coverage:** audit findings F1–F3, F5 and F6 are closed (AUDIT.md). Two gaps remain, neither affecting the
  current data:
  - **G1 (open, low, latent):** rows are not compared with each other, so a curated file in which one row copies or
    reverses another would apply once and then fail on the second run. The tests on the real file would catch it.
  - **Scaled duplicates:** a duplicate with the same metabolites but scaled coefficients is not detected. Both
    reviews found no existing reaction that shares a metabolite set with any of the 17.

## 7. Meeting list for Andrew

1. R351 kept open (no counterpart in his design).
2. DGAT GPR: DGA1 alone (reviewers) or DGA1 or DGA2 (kept, his rule).
3. His model re-adds two reactions we removed under authorized curations: ALCD2irm (R132) and G3PD1irm (R349).
4. The OLE1 loss from the generic-product design, and how he would restore the composition requirement.
5. Our reversible xPOOL_AC/FA pools and the +0.139 % artefact.
6. m1735[C_lp] charge −4: apply his `linoleoyl_neutral_correction`?
7. BIOMASS_C_core and BIOMASS_C_full differ from biomass_C in 18 and 17 terms (§5).
8. Crosswalk links missing for species that are ours: 4abut_c → m133, gly_v → m1863, 3hpcoa_c → m1990,
   pac_c → m2026, 34dhoxpeg_c → m1954 and the other COMT species.
9. From the look-alike review: his CYOOm is our removed non-pumping R2206; RNR OR GPRs; R2069.
10. Our two cytosolic GABA species: m133, which our reactions use, and m2040, the one his crosswalk links.
11. Cytosolic COQ2, held by the CoQ9 gate.
12. Which subsystem he wants for the 15 lipid steps (his reactions have none; we used "Glycerolipid metabolism").
13. Database links (KEGG, Rhea, MetaNetX, EC) for the 17 reactions, which carry only an SBO term (memote, §11).

## 8. Tests

`iyali26_gem/.venv` (Python 3.13.5, cobra 0.30.0), `python -m pytest --continue-on-collection-errors` from
`platform/`. In every run, 4 `test_gem_review_lit_*` modules cannot be collected (no `pydantic_ai` / `httpx2` in this
venv).

| Run | Passed | Failed | Errors | Skipped |
| --- | ---: | ---: | ---: | ---: |
| Before (HEAD `b57272d3`, before the annotation round) | 426 | 3 | 23 | 24 |
| v3 code (`048e18ec` plus the first step) | 450 | 3 | 23 | 24 |
| v4 code (step `b9473e58`) | 451 | 3 | 23 | 24 |
| Final code (step `3e1b593b`) | 451 | 3 | 23 | 24 |

Compared test by test:
- **Against the before run:** the only differences are new tests, all passing:
  - 17 from the committed annotation round (`test_bigg_crosswalk_annotations`, `test_export_bigg_ids`);
  - 8 in `test_mcnaugh_reactions`.
- **v4 against v3 code:** one test was added (rewritten metadata).
- **Final against v4 code:** no test changed outcome. Subtests went from 1,120 to 1,122: two were added to an
  existing test (reverse duplicate; earlier duplicate in an applied model).
- **Failing and erroring tests:** the same set in every run (see [logs/](logs/)).

**The new guard cases are not vacuous.** Run against the v4 guard code, the earlier-duplicate case returns
`already_correct` and the reverse case applies the change. Both fail there and pass with the final code.

## 9. Compute

- **Review stage:**
  - `artifacts/scripts/build_facts.py`: one WT solve, then per addable reaction one growth solve and one two-sided
    FVA, then the combined and lumped checks.
  - The reviewer and skeptic agents ran their own LPs (for example, one skeptic reports 20).

  These were not tallied into one count.
- **Verification:**
  - `verify_v3.py`: FVA over the 17 new reactions, plus 3 growth solves;
  - `wt_energy.py`: 2 LPs on v3 and 2 on v4;
  - two certified screens on v3, each 1,082 sweep + 1 ceiling + 40 certificate LPs;
  - one certified screen on v4 (biomass_C, 2026-11-04), with the same LP counts, 1 min 47 s.
- **OLE1 diagnosis:** in-memory solves on v3 (forward-only pools, WT and OLE1 KO). They were not logged, and the
  count is unknown.
- **Audits:** the first audit ran 9 LPs. The two fix reviews ran no LPs.
- **Builds:** three. v3 and v4 used the budgeted build and its spare. The third, which reproduced v4 with the final
  code, was approved by the user ("Fix and rebuild").
- **Memote ("do a memote test"):**
  - two snapshots run at the same time, v2 and v4, 1:12 and 1:13 wall time;
  - direct calls of the two memote functions that raised errors;
  - 9 LPs on the default-medium growth of newBiom, xBIOMASS and biomass_C.

## 10. Files

- **Repository:** the code, tests and curated data in §1, plus this folder:
  - [TASK.md](TASK.md), [AUDIT.md](AUDIT.md) and [before.json](before.json);
  - the v3 and v4 verification files;
  - [verification_inputs_sha256.txt](verification_inputs_sha256.txt);
  - [logs/](logs/): build commands and pytest summaries.
- **Artifacts:**
  - `inputs/`: facts and review groups;
  - `review_output.json` and `review_summary.txt`;
  - `scripts/`: facts, generator, verification, diff and objective-wrapper scripts;
  - `verification/`, `logs/`;
  - `superseded_ec7e3465/`: the v3 inputs and the fix diffs;
  - `superseded_guard_b9473e58/`: the v4 guard code and its diff;
  - `memote/`: the two memote snapshots (§11), the comparison and its README;
  - `screen_v4_20261104/`: the v4 screen (§4), its input SHAs, the comparisons with v3 and v2, the re-score and its
    README;
  - `SHA256SUMS`.
- **Publishing:** committed on branch `claude/lipid-unlump-collage-zip-merge-a59051` on the user's "Yes, commit"
  (2026-11-04); not pushed. Left uncommitted: the superseded v3 candidate `…_mcnaugh_rxns_20261103.*` (`ec7e3465`), the
  `…_mcnaugh_rxns_v2_20261103.*` copy (the same bytes as v4), the superseded pre-audit annotation candidate
  `…_bigg_annot_20261103.*`, and the auditor's `platform/lps_out.json`.

## 11. Memote snapshot (2026-11-04)

**Setup.**
- memote 0.17.0 on v4 (`e584151f`) and the v2 baseline (`9158e9f2`), with Gurobi, a 10 s solver timeout and no
  experimental data.
- Memote's "default medium" is the exchange bounds stored in the file, not Task A.
- A memote score is not evidence for or against a model change. This section is a quality check only.
- Run details, the comparison script and one script defect are in `artifacts/.../memote/README.md`. The defect: each
  FVA pool re-ran the script in its workers. The final files are the parent results; their totals match what each
  parent printed.

| Score | v2 | v4 | Change |
| --- | ---: | ---: | ---: |
| Total | 0.603724 | 0.603836 | +0.000112 |
| Consistency (weight 3) | 0.492448 | 0.492841 | +0.000393 |
| Metabolite annotation | 0.571033 | 0.571033 | 0 |
| Reaction annotation | 0.746291 | 0.745714 | −0.000577 |
| Gene annotation | 0.565310 | 0.565310 | 0 |
| SBO (weight 2) | 0.695204 | 0.695208 | +0.000004 |

34 of the 77 tests differ (`memote_compare.txt`). Every difference traces to the 17 additions or the 4 closures.

**A rule for reading the flux tests.** Adding reactions can only enlarge the feasible flux space, so a reaction or
metabolite that loses capacity in v4 lost it through the closures (R353, R1846, R1771, R2165). This holds for the
open-bound tests. It also holds for the unbounded-flux test, which runs FVA at 100 % of the biomass_C optimum: that
optimum is 0 under the file's bounds in both models.

**Scored consistency tests:**
- **Mass balance:** 543 → 558 imbalanced reactions, adding the 15 lipid steps (R-group design, §5).
- **Charge balance:** 461 → 464, adding AGPAT182_lp, DGAT182_lp and G3PAT182_lp (m1735[C_lp], §5). The metric still
  drops slightly, from 0.21542 to 0.21531, because the reaction count grows.
- **Unbounded flux in the default medium:** 113 → 105 reactions.
  - The 8 that left are R1771 (closed) and R1728, R1766, R1776, R1778, R1785, R1787, R1789. Those are PE:DAG
    acyltransferase, lysoPE acyltransferase, and DAG, TAG, PE and lysoPE transport between the ER membrane and the
    lipid particle.
  - The closures caused this (the rule above).
  - That these 7 ran a loop through R1771 is inferred.
- **Unchanged:** stoichiometric consistency fails in both models; disconnected metabolites are 0 in both.

**Unscored flux tests:**
- **Blocked reactions:** 677 → 680, adding R353, R1771 and R1846, which are closed by design. R2165 was already
  blocked in v2. None of the 17 is blocked.
- **Metabolites that cannot be consumed with open bounds:** 527 → 530, adding three decanoyl-CoA species:
  m2010[C_cy], m1737[C_mm] and m1947[C_cy].
  - Their only reactions are R2167 (m2010 ↔ m1737), R2266 (m1947 ↔ m1737) and R2165. R2165 is the lumped TAG
    synthesis, which drew 0.02 decanoyl-CoA and is now closed.
  - All three were already not producible in v2, so no flux is lost.
  - Andrew's 15 steps have no C10:0 chain.
- **Metabolites that cannot be produced:** 454 in both, the same set.
- **Stoichiometrically balanced cycles:** 545 → 535. The same 8 reactions left, plus xPOOL_AC_EM and xPOOL_AC_LP.
- **Energy-generating cycles:** the same 7 energy couples fail in both models, ATP (MNXM3) among them. The failure
  predates this change and also appears in the June memote run
  (`artifacts/diagnostics/data/memote_summary.csv`).
  - Memote's method differs from this session's energy check, which gives 0. The two are not reconciled here.
  - The reaction lists, taken from one LP solution, changed. These are alternative optima and were not compared
    further.
- **Two tests raised errors in the v4 run:** unconserved metabolites and minimal inconsistent stoichiometries. Both
  are unscored. Called directly, each in its own process:
  - Unconserved metabolites: 1,253 in both models, the same set.
  - Minimal inconsistent stoichiometries: 10 sets each.
  - On v4 the unconserved-metabolite LP solves or reports "unbounded" depending on the solver settings. This is a
    numerics effect, not a changed result.

**The two new biomass reactions:**
- **Growth under the file's bounds:** BIOMASS_C_core 1.3592 and BIOMASS_C_full 1.3559, both passing. biomass_C
  gives 0 in both models; this predates the change, and on Task A it grows 0.3125.
- **Old biomass reactions grow faster in v4 under the file's bounds:**
  - newBiom 1.24910 → 1.26449 (+1.23 %); xBIOMASS 0.93015 → 0.93530 (+0.55 %).
  - At the v4 optimum, newBiom runs mostly on the C16:0 steps (G3PAT160, AGPAT160, DGAT160) plus AGPAT182 (observed).
  - With the 15 lipid steps closed, growth is 0, because the lumped steps are closed too.
  - This is the same free chain choice as the +0.139 % artefact and the OLE1 loss (§5). The mechanism is inferred.
- **Blocked precursors:** exactly the 20 charged tRNAs, in both reactions and in both media.
  - Memote tests each precursor with its own demand. A charged tRNA can only be drained if its free tRNA comes back,
    and only a biomass reaction returns it (inferred).
  - Both reactions grow, so this is a test artefact of a tRNA-based biomass.
- **Missing "essential precursors":** 35 for each reaction:
  - the 20 free amino acids, which enter as charged tRNAs;
  - GTP, CTP, UTP, dATP, dCTP, dGTP, dTTP, FMN, PLP, NAD, NADP, FAD and TDP;
  - two MNX IDs (MNXM12, MNXM16) that memote did not map.

  biomass_C misses 31, the same set without NAD, NADP, FAD and TDP. Those four are among the 17 terms the new
  reactions lack (§5).
- **Not computed:** biomass molar mass gets no value for the two new reactions, as for biomass_C and xBIOMASS in both
  runs. The fast-growth and direct-metabolite tests give the same outcome as for all four old biomass reactions.
  Growth-associated maintenance is present.

**Annotation and structure:**
- **Reaction annotation:** the 17 reactions carry only an SBO term, as in the curated JSON. They therefore count as
  missing in all 9 reaction databases, which is the −0.00058.
- **Shared GPRs:** the chain copies share their step's GPR, as expected:
  - G3PAT joins the YALI1C00230g group (4 → 9 reactions);
  - the 5 DGAT copies share `YALI1D10264g or YALI1E38810g`;
  - the 5 AGPAT copies join R1843 under `YALI1E22736g or YALI1F25951g`.
- **Reactions without a GPR:** adds the two biomass pseudo-reactions.
- **Duplicates:** 13 in both, the same set.
- **Matrix:** rank 1,777 → 1,778; degrees of freedom 552 → 568.

**What this means.** The score barely moves. Every change has a known cause in this round's design. Two of them were
already listed as limits: the 15 lipid steps' R-group imbalance and the free chain choice. Memote raises one new
follow-up: database links for the 17 reactions, for the meeting with Andrew.
