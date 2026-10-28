# TASK: join the split dolichol route (census UM22 + UM167), close R2178, drain Dol-P (2026-10-28)

## 0. Authorization

The user, 2026-10-28, after reading the Stage 1 draft of this plan (shared notes folder,
`dolichol_chain_TASK_draft_20261028.md`; its section 9 gives the full approval sentence):

> Approve Option B

"Option B" is section 6 of the draft and the approval sentence of section 9:

> Approve Option B: close R2178, add the lumped gene-less dehydrodolichol-PP → dolichol step with the
> three-step-detour net stoichiometry, and add dolichyl phosphate m503[C_er] to biomass_C at 1e-4 mmol/gDW,
> through a guarded --dolichol-route step on ce9e1112 with the test, builds, diff, Gurobi 1e-9 knockout sweep,
> re-census, code review, audit, REPORT.md and STATE.md entry of section 7. No commit or push.

Sections 1–9 below are the draft unchanged apart from this section, the status line and the paths of this
folder. The task was started from the gap-fix line (parent session) as the next census target.

**Status: done 2026-10-28. Stopped at the knockout gate ([STOP.md](STOP.md)), accepted with the gate restated by the user (section 11), evidence corrected after the audit (section 12); final candidate `9d655c2b…`. Committed on branch `claude/dolichol-chain-gap-dfec4d` on the user's instruction, not pushed.** Stage 1 outputs (scripts, JSON, logs) are in
[stage1/](stage1/) and [scripts/](scripts/); `before.json` records the ce9e1112 state of every target.

## 1. Starting state (read 2026-10-28)

- **Branch / HEAD:** `claude/dolichol-chain-gap-dfec4d` at `b2a66dd` (= `main` after PR #26), clean.
- **Baseline model:** `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_20261027.xml`,
  SHA-256 `ce9e1112ee940531469457009234a2fd11de7f589c5948a7dec6b91118005294` (re-hashed). Built 2026-10-28
  with the argv in its `.build.json` (last flags `--cofactor-biomass-terms --r2135-acyl-coa-pool
  --r1726-cpt-open`). Newest candidate on `main`; **not** a designated reference.
- **Reference model:** none is designated (`model/candidates/README.md`).
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` (`IYALI26_RESEARCH_ROOT`).
- **Condition and solver:** Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` +
  `state/strain_profiles/po1f_sd_leu.json`; Gurobi 13.0.2 at FeasibilityTol / OptimalityTol 1e-9 through
  `model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py` (reset before each solve,
  knocked-flux check). The open-condition blocked scans use the census scripts (GLPK, 1e-7) in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_metabolites_ponce2013_20261026/scripts/`.
- **Census baseline (recalled from the Kennedy re-run census, 2026-10-27):** 709 blocked reactions in the
  open condition on ce9e1112; UM6 (19), UM22 (6) and UM167 (1) are all still blocked (re-verified below).

## 2. The three fragments as they are in ce9e1112 (verified: `map_ce9e1112.txt`, `provenance.txt`)

| Fragment | Reactions | Compartment | Equation(s) | Bounds | GPR | Provenance |
| --- | --- | --- | --- | --- | --- | --- |
| **UM6** | R217–R235 "cis-prenyltransferase step 01–19" | `C_lp` | R217: 3 H+ + FPP `m411` + IPP `m412` → PPi + GGPP `m414`; R218–R235: + IPP each, up to docosaprenyl-PP `m432` (C110, 22 units) | (0, 1000); R230 / R231 reversible | `YALI1C26036g or YALI1D21543g` | Yeast 7 / yeast-GEM `r_0281–r_0299` (same names, same GPR shape `YBR002C or YMR101C`, same odd reversibility of steps 14–15); identical in `iyli21.xml`, `iyali26.xml`, legacy, ce9e1112 (R217 gained 3 H+ from the proton gate) |
| **UM22** | R2179, R2180, R2178, R2177, R2182, R2181 | `C_cy` / `C_er` | R2179: GPP + IPP → FPP′ `m2020` + PPi (duplicate of R322, different FPP species, no formula); R2180: FPP′ + IPP + 2 H+ → dehydrodolichol-PP `m2019` + PPi (lumped, "adding 10–55 units", seed rxn15542); **R2178 "NO NAME": 6.5 O2 + `m2019` → Dol-PP `m2018` + 3 H2O + 5 CO2** (no GPR, no annotation, no name); R2177 / R2181: Dol-PP + H2O → Dol-P + Pi + 3 H+ (cy / er, BiGG `DOLDPP_Ler`, EC 3.6.1.43); R2182: Dol-PP cy ⇌ er | all (0, 1000); R2182 reversible | R2179 `YALI1D21022g or YALI1E06759g`; R2180 `YALI1C26036g or YALI1D21543g`; R2177 / R2181 `YALI1D26639g`; R2178, R2182 none | iYli21 additions (R2130–R2200 block, BiGG / HMR ids); identical in all four models; R2178 has no annotation in the raw `iyli21.xml` either |
| **UM167** | R287 "dolichol kinase" | `C_cy` | CTP + dolichol `m499` → CDP + Dol-P `m501` + H+ | (0, 1000) | `YALI1E22718g` | Yeast 7 / yeast-GEM `r_0360` (`YMR013C` SEC59); identical in all four models |

**Dolichyl phosphate is a carrier that nothing drains** (verified): `m501[C_cy]` is in R1101 (cy ⇌ er), R2177,
R287, R288; `m503[C_er]` is in R2102 (Dol-P-Man synthase), R289 (PMT, returns Dol-P), R372 / R374 / R375 (GPI
mannosyltransferases, return Dol-P), R2106 (ALG7, Dol-P → GlcNAc-PP-Dol; blocked, UM18), R2138 (ALG5 Dol-P-Glc;
blocked, UM18), R2181. At the Task A optimum Dol-P cycles through R2102 / R289 at 0.0735 mmol/gDW/h and
`biomass_C` has no dolichol term (its 61 species contain none).

**yeast-GEM 9.0.2 has the same gaps** (verified on the local copy, SHA-256 `9fd2c572…`): farnesyl diphosphate
[lp] has a transport `r_4604` but IPP [lp] is made only by the reversible steps 14–15; docosaprenyl diphosphate
[lp] (`s_0640`) is in `r_0299` only; dolichol [c] (`s_0642`, placeholder formula C20H36O) is in `r_0360` only (no
producer); dolichyl phosphate is in no pseudoreaction. So UM6 and UM167 are inherited Yeast 7 fragments and the
dehydrodolichol-PP block is an iYli21 addition that was never joined to them.

**Species formulas (verified):** `m499` dolichol has formula C15H28O (a farnesol formula; the builder's
MetaNetX name match was rejected for the C count); `m501` / `m503` Dol-P are C20H37O4P(C5H8)n; `m2018` /
`m2021` Dol-PP are C25H46O7P2; `m2019` dehydrodolichol-PP and `m2020` FPP′ have no formula. Mass balance of any
reaction on these species cannot be checked.

## 3. Biology and gene identity

| # | Claim | Source | Tier |
| --- | --- | --- | --- |
| B1 | `YALI1C26036g` = YALI0C18799g, KEGG yli:2909965, RefSeq "Dehydrodolichyl diphosphate synthase complex subunit RER2", KO K11778, EC 2.5.1.87, 267 aa, UniProt Q6CBH1. On the consensus essential list ("similar to RER2"). | [KEGG yli:2909965](https://www.kegg.jp/entry/yli:2909965); `s2_metabolic_genes.csv`; consensus list | verified (opened); orthology, not experiment |
| B2 | `YALI1D21543g` = YALI0D17556g, KEGG yli:2911304, RefSeq "Dehydrodolichyl diphosphate synthase 6", K11778, EC 2.5.1.87, 326 aa. Not on the consensus list. | [KEGG yli:2911304](https://www.kegg.jp/entry/yli:2911304) | verified; orthology |
| B3 | A **NUS1 ortholog exists and is not in the model**: `YALI1F26726g` = YALI0F20108g, KEGG yli:2908051, "DDS complex subunit NUS1", K19177, 301 aa; S2 `in_model = no`; **on the consensus essential list**. | [KEGG yli:2908051](https://www.kegg.jp/entry/yli:2908051); S2; consensus list | verified; orthology |
| B4 | `YALI1D26639g` = YALI0D21032g, KEGG yli:2911055, "Palmitoyl-protein thioesterase-dolichyl pyrophosphate phosphatase fusion 1", K07252 dolichyldiphosphatase EC 3.6.1.43, 216 aa. Not on the consensus list. | [KEGG yli:2911055](https://www.kegg.jp/entry/yli:2911055) | verified; orthology |
| B5 | `YALI1E22718g` = YALI0E18942g, KEGG yli:2911738, "Dolichol kinase", K00902, EC 2.7.1.108, 588 aa. **On the consensus list** ("similar to SEC59"). | [KEGG yli:2911738](https://www.kegg.jp/entry/yli:2911738) | verified; orthology |
| B6 | A **polyprenol/polyprenal reductase ortholog exists and is not in the model**: `YALI1E09357g` = YALI0E07788g, KEGG yli:2912457, "Polyprenol reductase", K12345 (SRD5A3 / DFG10), EC 1.3.1.22 / 1.3.1.94, 333 aa; S2 `in_model = no`; not on the consensus list. | [KEGG yli:2912457](https://www.kegg.jp/entry/yli:2912457); [KEGG EC 1.3.1.94](https://www.kegg.jp/entry/ec:1.3.1.94) (R13444: dolichal + NADP+ = polyprenal + NADPH + H+; lists sce YIL049W and yli 2912457) | verified; orthology |
| B7 | *S. cerevisiae* Rer2 (P35196) and Srt1 (Q03175) each form the dehydrodolichyl diphosphate synthase complex with Nus1 (Q12063); Rer2 is ER-membrane, Srt1 lipid droplet, Nus1 both. Sec59 (P20048) is the dolichol kinase (final step of Dol-P synthesis). Cax4/Cwh8 (P53223) is a **non-essential** dolichyl-diphosphatase that recycles Dol-PP. Dfg10 (P40526) reduces polyprenal to dolichal with NADP. | UniProt REST, reviewed entries (fields function, location, phenotype) | verified (entries read); curated, mostly *S. cerevisiae* experiment |
| B8 | Rer2 and Srt1 are two cis-prenyltransferases: Srt1 products are longer (about 21 units vs about 15 for Rer2), are "not sufficiently converted to dolichol and dolichyl phosphate", and Srt1 is "mainly localized to lipid particles"; RER2 is expressed in early log phase, SRT1 in stationary phase. | Sato et al. 2001, Genes Cells 6:495, [PMID 11442630](https://europepmc.org/abstract/MED/11442630), doi 10.1046/j.1365-2443.2001.00438.x (abstract read); Hoffmann et al. 2017, Genetics 207:1371, [PMC5714454](https://pmc.ncbi.nlm.nih.gov/articles/PMC5714454) (read: Rer2 also in lipid droplets; long Srt1 polyprenols not efficiently converted) | verified (abstract / article read) |
| B9 | Δrer2 is viable but "grew very slowly even at 23 °C", with about 2.7 % of wild-type cis-PT activity; Δsrt1 grows normally; the Δrer2 Δsrt1 double disruption is lethal. SGD: rer2Δ "inviable" in the S288C systematic set (Giaever 2002) but viable with a severe growth defect (Sato 1999); srt1Δ viable; nus1Δ inviable; sec59Δ inviable (Giaever 2002; Bernstein 1989); cax4Δ viable with decreased growth; dfg10Δ viable (fitness 0.94). | Sato et al. 1999, MCB 19:471, [PMC83905](https://pmc.ncbi.nlm.nih.gov/articles/PMC83905) (read); SGD phenotype records for [YBR002C](https://www.yeastgenome.org/locus/YBR002C), [YMR101C](https://www.yeastgenome.org/locus/YMR101C), [YDL193W](https://www.yeastgenome.org/locus/YDL193W), [YMR013C](https://www.yeastgenome.org/locus/YMR013C), [YGR036C](https://www.yeastgenome.org/locus/YGR036C), [YIL049W](https://www.yeastgenome.org/locus/YIL049W) (backend phenotype pages read) | verified (*S. cerevisiae*) |
| B10 | SRD5A3 (human) and its yeast ortholog Dfg10 are required to reduce the α-isoprene of polyprenol to make dolichol; depleted cells still contain some dolichol ("unexpected alternative pathway"). | Cantagrel et al. 2010, Cell 142:203, [PMID 20637498](https://europepmc.org/abstract/MED/20637498) (abstract read) | verified (abstract) |
| B11 | **2026 revision:** the three-step detour is conserved in yeast: Tda5 (DHRSX ortholog) oxidises polyprenol to polyprenal, Dfg10 reduces polyprenal to dolichal, Env9 reduces dolichal to dolichol (NADPH); tda5Δ lowers dolichol and accumulates polyprenol; DFG10 / SRD5A3 do not rescue tda5Δ. | Hanaoka et al. 2026, PNAS 123:e2613147123, [PMID 42201967](https://europepmc.org/abstract/MED/42201967), doi 10.1073/pnas.2613147123; Wilson et al. 2026, bioRxiv doi 10.64898/2026.04.22.720056 (abstracts read) | verified (abstracts); TDA5 / ENV9 orthologs in *Y. lipolytica* **not identified** here |
| B12 | The enzyme(s) that dephosphorylate polyprenyl diphosphate to polyprenol are not named in any source read; the reviews describe this step as uncharacterised. | B8, B10, B11 sources; [KEGG EC 3.6.1.43](https://www.kegg.jp/entry/ec:3.6.1.43) covers only Dol-PP → Dol-P (R01004) | inferred (absence in the sources read, not an exhaustive search) |
| B13 | Dolichol content: about 6 µg/g wet weight for yeast; *S. cerevisiae* 32 µg/g (late log) and 48 µg/g (early stationary), basis not stated in the sentence (the paragraph's other figure is wet weight), citing Surmacz & Swiezewska 2011 BBRC (unread). Dol-15 / Dol-16 dominate (Dol-14–18 family), with a second Dol-19–22 family. Sato 2001 reports roughly equal dolichol and Dol-P (search summary, not read). | Gryz et al. 2019, IJMS 20:3043, [PMC6628320](https://pmc.ncbi.nlm.nih.gov/articles/PMC6628320) (passage read) | verified for the 2019 passage; the primary measurement is unread; no dry-weight value found (two searches) |
| B14 | No *Y. lipolytica* paper on dolichol, polyprenol or cis-prenyltransferase was found (Europe PMC title and abstract searches, 0 hits). | this session | verified (searches run); absence, not evidence |
| B15 | The BiGG metabolite `dedoldp` (dehydrodolichol diphosphate) is in iRC1080, iCHOv1 and iLB1027_lipid, formulas C55H89O7P2 / C95H153O7P2; BiGG lists no reaction for it. | [BiGG dedoldp](http://bigg.ucsd.edu/universal/metabolites/dedoldp) (legacy site, read) | verified |

**What follows (inferred):**
- **UM22 + UM167 is the real route** (Rer2 / Nus1 → polyprenyl-PP → polyprenol → dolichol → Sec59 → Dol-P),
  split by the missing dephosphorylation / reduction steps, which the model replaces with the orphan R2178
  (an O2-consuming, CO2-releasing equation that matches no enzyme in B7–B12) and never connects to dolichol.
- **UM6 is the Srt1 lipid-droplet chain**, not a duplicate: its product (22-unit polyprenyl-PP) is a different
  molecule that *S. cerevisiae* does not efficiently convert to dolichol (B8). Its dead end is faithful to the
  biology as far as the model's species go; what is missing is a lipid-droplet polyprenol pool, for which no
  consumer or measurement exists in the model's sources. The census fix H5 (IPP + FPP import plus a sink on
  docosaprenyl-PP) unblocks it only through an artificial sink.
- **The GPR `YALI1C26036g or YALI1D21543g`** on R217–R235 and R2180 is the Yeast 7 `RER2 or SRT1` shape. B3
  says the Nus1 partner is missing from the model. By AGENTS.md §5 an OR member needs evidence that it
  catalyses the step on its own: Srt1 does (B9, residual activity, double null lethal), but its product is the
  long one. Not changed in this task.

## 4. In-memory LP results (Task A + PO1f, Gurobi 1e-9, checked; GLPK open-condition scan)

Scripts: `scripts/stage1_dolichol.py` (1,148 LPs → `stage1/stage1_dolichol.json`), `scripts/stage1_optionB.py`,
`scripts/stage1_optionB_net.py`, `scripts/stage1_dolichol_census.py` (→ `stage1/stage1_dolichol_census.json`).

**Variants.** `as_built`; `dolp_demand` (test demand on `m503[C_er]`); `dolp_biomass_4e-5` and `_1.6e-4`
(`m503[C_er]` in `biomass_C`); `dfg10_route` (biomass 4e-5 + R2178 closed + test step `m2019` + H2O + NADPH +
H+ → dolichol `m499` + PPi + NADP+); `h5_um6` (biomass 4e-5 + census H5); `optionB_*` (R2178 closed + the test
step + biomass at 4e-5 / 1e-4 / 1.6e-4); `optionB_net` (biomass 1e-4; test step with the three-step net
`m2019` + H2O + NAD+ + 2 NADPH + H+ → `m499` + PPi + NADH + 2 NADP+); `no_new_step` (biomass 1e-4 + R2178 closed,
nothing added).

| Variant | WT growth (status) | No-uptake max R_NTP1 / max biomass (maintenance relaxed) | Route reactions unblocked (Task A, FVA ≥ 0) | SEC59 `E22718g` KO | CWH8 `D26639g` KO |
| --- | ---: | --- | --- | ---: | ---: |
| as_built | 0.3116824 (optimal) | 0 / 0 (optimal) | none of UM6, UM22, UM167 | 1.00 | 1.00 |
| dolp_demand | 0.3116824 | 0 / — | UM22 all 6 (FVA max 4.6e-4) | 1.00 | 1.00 |
| dolp_biomass_4e-5 | 0.3116740 (−2.7e-5 rel.) | 0 / — | UM22 all 6; **R2178 forced to 1.25e-5 (min = max)** | 1.00 | **0.00** |
| dolp_biomass_1.6e-4 | 0.3116489 (−1.1e-4 rel.) | 0 / — | same; R2178 forced to 5.0e-5 | 1.00 | **0.00** |
| dfg10_route (4e-5) | 0.3116738 | 0 / — | R2179, R2180, R287 (+ test step) at 1.25e-5; R2177 / R2181 / R2182 stay 0 | **0.00** | 1.00 |
| h5_um6 (4e-5) | 0.3116740 | 0 / — | UM6 all 19 (FVA max 8.4e-5, 0 at the optimum) + UM22 | 1.00 | 0.00 |
| optionB_4e-5 / 1e-4 / 1.6e-4 | 0.3116738 / 0.3116608 / 0.3116479 | 0 / 0 (optimal, all three) | R2179, R2180, R287, test step (min = max = 1.25e-5 / 3.1e-5 / 5.0e-5) | 0.00 | 1.00 |
| optionB_net (1e-4) | 0.3116608 | 0 / 0; **also 0 / 0 with the test step made reversible** | same | 0.00 | 1.00 |
| no_new_step (1e-4) | **0 (optimal)** | 0 / 0 | — | — | — |

Other knockouts (all variants): `YALI1C26036g` and `YALI1D21543g` 1.00 (OR rule); the mevalonate / FPP
genes `YALI1E06759g`, `YALI1F06018g`, `YALI1F08363g`, `YALI1D21022g` and the Dol-P-Man synthase genes
`YALI1D08050g`, `YALI1B07538g`, `YALI1A19526g` are lethal before and after (unchanged); the ALG7 / ALG13 /
ALG14 / ALG5 / PMT genes stay 1.00 (their reactions stay blocked, UM18).

**Corrected energy check:** every exchange (0, 1000), the only positive lower bound (`xMAINTENANCE` 7.8625)
relaxed to 0, R_NTP1 maximised: **optimal at 0 in every variant** (status recorded in the JSON). Max biomass
with no uptake is also optimal at 0 for the Option B variants.

**Open-condition blocked scan (GLPK 1e-7, census method, every exchange open; `stage1_dolichol_census.json`):**

| Variant | Blocked | Freed vs as built | Newly blocked vs as built |
| --- | ---: | --- | --- |
| as_built | 709 | — | — |
| dolp_demand | 703 | R2177, R2178, R2179, R2180, R2181, R2182 (UM22) | none |
| dolp_biomass_4e-5 | 705 | the same 6 | R1343, R1344 |
| dfg10_route (Option B topology) | 708 | R2179, R2180, R287 | R1343, R1344 |
| h5_um6 | 686 | the 6 + all 19 of UM6 | R1343, R1344 |

R1343 / R1344 are the pair the cofactor task (REPORT.md 2026-10-27, "GLPK reports R1343, R1344 and R970 newly
blocked") classed as a GLPK 1e-7 fringe that Gurobi at 1e-9 does not reproduce; they appear here with every
biomass term, demand-only excluded. Stage 2 must re-check them with Gurobi. Under Option B the three UM22
reactions R2177 / R2181 / R2182 stay blocked on purpose: Dol-PP has no producer until the N-glycan chain
(UM18) is linked, and R2178 moves from "blocked" to "closed".

## 5. What the Stage 1 question asked

- **Where the biology says the route runs:** ER membrane (Rer2 / Nus1 complex, Sec59, the Dol-P-using
  glycosyltransferases), with Rer2 and especially Srt1 also on lipid droplets (B7, B8). The cytosolic ids of
  R2179 / R2180 / R287 are the model's convention for ER-membrane enzymes facing the cytosol (as for the
  mevalonate route); the Dol-P transport R1101 cy ⇌ er already exists.
- **Which fragment is real:** UM22 (R2179, R2180, R2177 / R2181 / R2182) plus UM167 (R287) together, once
  joined by the reduction step (B10, B11).
- **What the others are:** R2178 is an iYli21 placeholder with no gene, name, annotation or recognisable
  chemistry (B12 sources name no O2-dependent step); UM6 is the Yeast 7 Srt1 chain whose long product is a
  real but separate molecule with no modelled fate (B8). R2179 duplicates R322 with a formula-less FPP
  species `m2020` (not changed here; a merge is a separate cleanup).

## 6. Options

| | Change | For | Against |
| --- | --- | --- | --- |
| A | Dol-P `m503[C_er]` in `biomass_C` only | one coefficient; frees all 6 of UM22; matches the census's "cofactor pseudo-gap" remedy | makes the orphan R2178 obligatory (min = max); predicts CWH8 essential and SEC59 dispensable, the reverse of B9 and of the consensus list; leaves dolichol / Sec59 unconnected |
| **B (recommended)** | R2178 → (0, 0) with the reason recorded; one lumped irreversible gene-less step `m2019` + H2O + NAD+ + 2 NADPH + H+ → dolichol `m499` + PPi + NADH + 2 NADP+ (three-step detour net, B11), annotated with the candidate orthologs (DFG10 `YALI1E09357g`; TDA5 / ENV9 not identified); Dol-P `m503[C_er]` in `biomass_C` at 1e-4 mmol/gDW | the route runs as the biology describes (polyprenyl-PP → dolichol → Sec59 → Dol-P); SEC59 essential and CWH8 dispensable, as in *S. cerevisiae*; R2177 / R2181 keep their real role (Dol-PP recycling, active once the N-glycan chain UM18 is linked); no EGC even if the step were reversible; growth −5e-5 relative | frees only 3 reactions now (R2179, R2180, R287) and closes one; the lumped step has no gene (a single-gene GPR on a three-enzyme step would not meet AGENTS.md §5, and dfg10Δ is viable, B9–B10); formulas on `m499` / `m2019` are wrong or missing, so the step cannot be balance-checked; the coefficient rests on a wet-weight figure (B13) |
| C | B with three explicit reactions (polyprenol, polyprenal, dolichal species) and GPRs on the DFG10 ortholog plus TDA5 / ENV9 orthologs to be identified | gene-level predictions for the reductases | three new species without formulas, two ortholog searches, a NAD+/NADPH exchange across separate reversible-looking steps; larger than the smallest set |
| D | B plus H5 for UM6 (IPP / FPP into the lipid particle and a sink on docosaprenyl-PP) | frees 19 more | the sink is not biology; the Srt1 product has no modelled fate; not recommended |

Option A is not recommended: its 6 freed reactions include forcing flux through a reaction that describes
nothing. By AGENTS.md §2, freeing reactions is not a reason for a change.

## 7. Recommended implementation (Option B), following the curation-change skill

1. **TASK.md** in `model/reports/dolichol_route_20261028/` quoting the authorization; `before.json` with
   R2178, R287, R2179, R2180, R2177, R2181, R2182, `biomass_C` stoichiometry and the species `m2019`, `m499`,
   `m501`, `m503`, `m2018`, `m2021` (name, formula, charge, reactions).
2. **Curated data** `model/curation/dolichol/dolichol_route.json` (schema v1, `DOLICHOL-ROUTE-20261028`):
   the R2178 closure with its equation pinned; the new reaction's id (proposed `R_DOLSYN`), name, stoichiometry,
   bounds (0, 1000), group "N-Glycan biosynthesis", annotation (EC 1.3.1.94, KEGG R13444, candidate ortholog
   `YALI1E09357g` as a note, not a GPR); the biomass term (`m503[C_er]`, 1e-4, source B13 with its derivation
   and the wet-weight caveat); the evidence table of section 3 with tiers; the accepted before-state of
   `biomass_C` (ce9e1112) so the step fails closed if biomass_C differs.
3. **Guarded step** `platform/scripts/gem_annotate/dolichol_route.py` (`apply_dolichol_route`): fails closed
   unless R2178, R287 and the species match the recorded state, `biomass_C` is the ce9e1112 composition, the
   objective is `biomass_C`, and the new id is absent (or present with exactly this definition: idempotent).
   Flag `--dolichol-route`, requires `--r1726-cpt-open`, runs after it in `main.py`.
4. **Focused test** `platform/tests/test_dolichol_route.py` on an in-memory model: applies, re-applies
   (no-op), and every guard fails with its own message. `python -m pytest` before / after, compared test by test.
5. **Builds** (budget 3: control, seed 0, seed 1): the control without the flag must reproduce `ce9e1112` byte
   for byte; seeds 0 and 1 with the flag must be identical; `PYTHONHASHSEED` recorded in the logs.
6. **Whole-model diff** vs ce9e1112: only R2178 bounds + notes, the new reaction, `biomass_C` stoichiometry +
   notes. Anything else: stop.
7. **Built-file checks** (Gurobi 1e-9): WT; FVA of the route; no-uptake R_NTP1 and biomass (relaxed
   maintenance, status reported); the no-new-step control is not rebuilt (in-memory evidence above).
8. **Knockout sweep** 2 × 1,082 genes, reset per knockout, knocked-flux check; **stop condition:** any class
   change other than `YALI1E22718g` non-essential → essential, or any growth change above 1e-6 outside the
   route.
9. **Re-census** (GLPK, census scripts) on the open condition and Task A: expected 709 → 706 blocked in the
   open condition (R2179, R2180, R287 freed; R2178 counted as closed) apart from the R1343 / R1344 GLPK fringe,
   which Gurobi 1e-9 FVA must show unblocked; nothing else newly blocked.
10. **Code review** (fresh reviewer on the diff) and **independent audit** (TASK.md, diff, logs, deliverables);
    REPORT.md, AUDIT.md, `model/STATE.md` entry. No commit or push unless instructed.

**Acceptance criteria:** control byte-identical; seeds identical; diff confined to the four entries; WT within
1e-4 relative of 0.3116824; EGC and no-uptake biomass 0 and optimal; sweep as in step 8; census as in step 9.
**Compute budget:** 3 builds, about 2,200 LPs (sweep) + 2 census scans; about 60 Gurobi LPs for the checks.

## 8. Risks and limits

- **Chemistry:** the new step and the kinase R287 cannot be mass-balanced (polymeric or wrong formulas, section
  2). A formula set for the dolichol family (Dol-16: dolichol C80H132O, Dol-P C80H133O4P, Dol-PP C80H134O7P2,
  polyprenyl-PP C80H132O7P2) would touch about 20 reactions and is a separate chemistry task.
- **Gene level:** the lumped step is gene-less, so no knockout prediction for the DFG10 / TDA5 / ENV9
  orthologs. NUS1 (`YALI1F26726g`, consensus-essential) stays out of the model; an AND partner on R2180 /
  R217–R235 needs its own complex-evidence review (AGENTS.md §5). `YALI1C26036g` (RER2, consensus-essential)
  stays non-essential through the OR with `YALI1D21543g`; whether the Srt1-type product should feed the same
  dolichol step is exactly what B8 argues against, so a GPR change is not proposed.
- **Coefficient:** 1e-4 mmol/gDW is an order of magnitude from a wet-weight figure (B13: 32–48 µg/g, Dol-16
  1,110 g/mol, dry matter assumed 20–25 % of wet weight → 1.2e-4–2.2e-4; if the figure were dry weight,
  2.9e-5–4.3e-5). Growth differs by less than 1e-4 relative across 4e-5–1.6e-4, and essentiality does not
  change. The user may set another value.
- **Organism:** no *Y. lipolytica* dolichol data (B14); all evidence is orthology plus *S. cerevisiae*.
- **Not addressed:** UM6 (Srt1 chain), the UM18 N-glycan chain (which would give Dol-PP its producer and
  activate R2177 / R2181 / R2182 / ALG7), the R2179 / `m2020` duplicate, the `C5H8)n` formulas.

## 9. Approval sentence (to be quoted verbatim in TASK.md if given)

> Approve Option B: close R2178, add the lumped gene-less dehydrodolichol-PP → dolichol step with the
> three-step-detour net stoichiometry, and add dolichyl phosphate m503[C_er] to biomass_C at 1e-4 mmol/gDW,
> through a guarded --dolichol-route step on ce9e1112 with the test, builds, diff, Gurobi 1e-9 knockout sweep,
> re-census, code review, audit, REPORT.md and STATE.md entry of section 7. No commit or push.

## 10. Implementation findings (Stage 2, 2026-10-28, before any build)

Two details of the approved step were settled while writing the data file. Neither changes the decision.

1. **Proton count of R_DOLSYN: 2 H+, not the 1 H+ of the draft's textbook notation.**
   - The model writes NAD dehydrogenases without H+: NAD `m122` is C21H27 and NADH `m123` C21H29, both
     charge 0. Of the 55 cytosolic reactions with both, 23 carry no H+ and balance (e.g. R45, R130, and
     malate dehydrogenase R534, checked separately).
   - It writes NADPH reductases with one H+ consumed per NADPH: NADP `m41` is C21H29 and NADPH `m40`
     C21H30, both charge 0. Of the 63 cytosolic reactions with both, 34 do this, H-balanced with the charge
     off by 1 per NADPH (e.g. R16, R60, R134). Scan:
     [scripts/proton_convention_scan.py](scripts/proton_convention_scan.py) →
     `stage1/proton_convention_scan.json` (verified, this session).
   - The three-step net written in that convention is `m2019 + H2O + NAD + 2 NADPH + 2 H+ → dolichol m499 +
     PPi + NADH + 2 NADP`. With a Dol-16 formula set filled in memory (m2019 C80H132O7P2, m499 C80H132O), C,
     H, N, O and P balance and the charge is −2, the model's NADP charge convention (−1 per NADPH, as in R16 or the 16-NADPH R1393).
     The 1 H+ form leaves H +1 unbalanced (verified, `stage1/proton_variant_h1.json`, `_h2.json`).
   - Every LP result is identical for the two forms (Task A, Gurobi 1e-9): WT 0.3116608, R2179 / R2180 /
     R287 / R_DOLSYN FVA at 99.9 % [3.113e-5, 3.117e-5], no-uptake ATP and biomass 0 (optimal) also with
     the step reversible, SEC59 KO 0, CWH8 KO 1.0, RER2-type KOs 1.0 (verified, same files).
2. **Annotation of R_DOLSYN: SBO only.** EC 1.3.1.94 and KEGG R13444 cover only the Dfg10 step of the four
   lumped steps, so they are recorded in the notes and in the data file's `component_steps`, not as
   `ec-code` / `kegg.reaction` annotations, which would name the whole lumped step after one enzyme.

## 11. Gate decision after the stop (2026-10-28)

The run stopped at the knockout gate of section 7, step 8 ([STOP.md](STOP.md)). The user, after reading the
stop report:

> Accept, restate gate

This is the first option of STOP.md section 5: the candidate is accepted with the gate restated, for a change
that adds a biomass term, as "no class change other than `YALI1E22718g`, and no knockout gains absolute growth
above 1e-6". Both parts are met (STOP.md sections 2–3). The run continues with the independent audit,
REPORT.md, AUDIT.md and the STATE.md entry. No commit or push.

## 12. Erratum after the independent audit (2026-10-28)

Sections 1–9 are the Stage 1 draft as approved and are left as written. The audit ([AUDIT.md](AUDIT.md))
found two evidence errors in them, re-checked here against the sources:

1. **B11 attributed too much to Hanaoka et al. 2026.** Read again this session (abstract and full text
   PMC13229312, verified): Hanaoka et al. name TDA5 as the yeast ortholog of DHRSX (tda5Δ lowers dolichol,
   accumulates polyprenol, rescued by DHRSX but not DFG10 / SRD5A3). They name no cofactor, and they do not
   assign the dolichal → dolichol step in yeast; Env9 appears only as a tunicamycin-screen hit. The
   corrected attribution (verified, abstracts read):
   - **Human:** DHRSX does the first and third steps "as a NAD+-dependent dehydrogenase and as a
     NADPH-dependent reductase"; SRD5A3 does only the second (Wilson MP et al. 2024, Cell 187:3585,
     [PMID 38821050](https://pubmed.ncbi.nlm.nih.gov/38821050/)).
   - **Yeast:** Env9 catalyses the NADPH-dependent reduction of dolichal to dolichol, and Tda5, not Env9,
     converts polyprenol to polyprenal (Wilson MP et al. 2026, bioRxiv preprint, not peer reviewed,
     [doi 10.64898/2026.04.22.720056](https://doi.org/10.64898/2026.04.22.720056)). No cofactor is named for Tda5.
   - **So:** the 1 NAD + 2 NADPH of R_DOLSYN is sourced from human DHRSX and SRD5A3 (and KEGG R13444 for
     Dfg10). The NAD+ of the Tda5 step in yeast is inferred by orthology, not shown.
2. **B1–B6 locus mappings.** The KEGG pages show only YALI2 locus tags. The YALI1 ↔ YALI0 ↔ KEGG gene-id
   mapping comes from `$IYALI26_RESEARCH_ROOT/reference/locus_map/s2_metabolic_genes.csv` (rows for
   `YALI1_E22718g`, `YALI1_E09357g`, `YALI1_D26639g`, `YALI1_C26036g` re-read, verified). The KEGG claims
   themselves (KO, EC, names) hold.

**Model consequence.** The R_DOLSYN notes and the curation file's `component_steps` and `evidence` carried
the wrong attribution. They were corrected in the data file only, and the candidate was rebuilt to a new
file, `…_dolichol_route_v2_20261028.xml` (`9d655c2b…`, seeds 0 and 1 identical). Versus `a801aad0` the
whole-model diff is one entry, the R_DOLSYN notes, and the optimization problem is identical (all species,
stoichiometry, bounds, GPRs, objective and solver constraints; [v2_same_lp.json](v2_same_lp.json)). So
every LP result of the run applies to it unchanged. `a801aad0` is kept, superseded.
