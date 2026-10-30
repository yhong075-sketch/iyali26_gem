# R24 / R395 homocitrate steps: review on e1e7c868 (2026-10-29, no model change)

**Scope.** Investigation only ([TASK.md](TASK.md) §0). Baseline `e1e7c868…`, the aconitase-round candidate that is
still uncommitted on `claude/aconitase-aco2-split-1c492f` (not merged into `origin/main`). Nothing in the model,
curated data, code or tests was changed. Nothing was committed.

**Short answer.**
- **R24's gene.** The current gene, LYS4-like `YALI1E03387g`, is not supported for R24. The supported genes are
  aconitases: the ACO1-type `YALI1D11984g`, and possibly the ACO2-type `YALI1E17935g`. Whether the ACO2-type is
  active in *Y. lipolytica* is unknown.
- **R395** stays `YALI1E03387g`.
- **R24's annotation is wrong.** Its name and every cross-reference describe 2-methylcitrate dehydratase, and so do
  the cross-references of its product m76. The chemistry is homocitrate → cis-homoaconitate.
- **The 3 H+ are a protonation artefact,** with no effect on FBA.
- **Direction:** reversible is fine.
- **Knockouts (Task A medium):** only option D (ACO2-type alone) changes a class. `YALI1E17935g` goes TN → FP at 10 %
  and 15 %. Options B and C change nothing.

## 1. Evidence

Tiers: **verified** = I fetched or ran it this session and read the value, with the file in [evidence/](evidence/)
or [logs/](logs/) (SHA-256 in [evidence/SHA256SUMS.txt](evidence/SHA256SUMS.txt)). **inferred** = reasoned, no single
source. **recalled** = from memory or another report, not re-checked.

### 1a. Who catalyses homocitrate → cis-homoaconitate

| # | Claim | Tier |
| --- | --- | --- |
| E1 | Yeast and *A. fumigatus* homoaconitases (Lys4p, LysF) show no homocitrate activity ("n.d." in Table 1). They interconvert homoaconitate and homoisocitrate only (Lys4p: homoaconitate 0.29, homoisocitrate 0.5 U/mg). Fazius 2012, [PMC3556520](https://pmc.ncbi.nlm.nih.gov/articles/PMC3556520/), Table 1 and Results. | verified (Europe PMC full-text XML, [evidence/PMC3556520_fulltext.txt](evidence/PMC3556520_fulltext.txt)) |
| E2 | Homocitrate is converted only by aconitase plus homoaconitase together. In the coupled assay (AfAcoA + LysF + Lys12p) homocitrate gave NAD reduction, and AcoA + Lys12p alone gave none. With yeast Aco1p or Aco2p in place of AcoA, "homocitrate was converted to homoisocitrate in the presence of … Lys4p". The yeast part is "not shown". | verified (same file); the yeast part is the authors' statement of unshown data |
| E3 | Yeast Aco2p: 950-fold more active with homoaconitate than with aconitate (4.76 vs 0.005 U/mg). Aco1p: aconitate 23.7, homoaconitate 0.9 U/mg. The 240 nm assay "cannot distinguish" the hydration direction. That aconitases make homocitrate rather than homoisocitrate rests on E2. | verified (Table 1, Results) |
| E4 | *S. cerevisiae* genetics: *aco1* and *aco2* single mutants grow without lysine but more slowly (*aco2* more so). The *aco1 aco2* double mutant "was completely auxotrophic for lysine". So both aconitases contribute and either alone is enough to grow. | verified (Results) |
| E5 | Filamentous fungi: the TCA aconitase "seems to be responsible for homocitrate dehydration". The Aco2-like AcoB proteins were inactive in vitro and did not complement yeast mutants. | verified (Discussion) |
| E6 | *T. thermophilus*: HACN shows no homocitrate dehydratase activity, while aconitase dehydrates (R)-homocitrate. Jia 2006, [PMID 16524361](https://pubmed.ncbi.nlm.nih.gov/16524361/). | verified, abstract only ([evidence/pubmed_abstracts_selected.txt](evidence/pubmed_abstracts_selected.txt)) |
| E7 | `YALI1E03387g` = YALI0E02728g = [Q6C791](https://rest.uniprot.org/uniprotkb/Q6C791.json), Swiss-Prot reviewed: "Homoaconitase, mitochondrial", gene LYS4, EC 4.2.1.36. Its only catalytic activity is RHEA:15485 (homoisocitrate ⇌ cis-homoaconitate + H2O). Function and location are ECO:0000250 (by similarity); protein existence is "inferred from homology". 69.5 % identity to Sc Lys4p P49367. | verified (UniProt JSON; locus map; [logs/aco2_k610_position.log](logs/aco2_k610_position.log)) |
| E8 | `YALI1D11984g` = Q6C9P6, ACO1-type (PANTHER SF3, 77.7 % to Aco1p). `YALI1E17935g` = Q6C5V2, ACO2-type (PANTHER SF2 "homocitrate dehydratase, mitochondrial", 70.3 % to Aco2p). | verified this session (UniProt JSON); the identity percentages are recalled from the [aconitase round](../aconitase_aco2_split_20261029/evidence/pairwise_alignment.tsv) |
| E9 | Fazius' specificity residue: Aco2p K610R raised aconitate activity 115-fold (0.005 → 0.575 U/mg) and kept 56 % of homoaconitate activity (Table 2). At the aligned position, Aco1p, `YALI1D11984g` and Af AcoA have **R**. `YALI1E17935g` (K609) has **K**, but so do the **inactive** Af AcoB and An AcoB. The K therefore marks the Aco2 clade, not an active enzyme. | Table 2: verified. Alignment: verified (homology only; RefSeq XP_750430.2, version used by the paper unknown) |
| E10 | *Y. lipolytica* has no enzymology, knockout or localisation data for either aconitase on homocitrate. Screened: 4 PubMed queries, 35 titles ([evidence/pubmed_screen_titles.txt](evidence/pubmed_screen_titles.txt)). The *Y. lipolytica* aconitate-hydratase kinetics paper (Rakhmanova 2023, [PMID 37108831](https://pubmed.ncbi.nlm.nih.gov/37108831/)) tests citrate and isocitrate only (abstract). | verified for that search; absence beyond it is unknown |
| E11 | Whether *Y. lipolytica* follows the yeast pattern (both aconitases) or the respiring-fungus pattern (ACO1-type only) | **inferred; open.** Pointing to yeast: it carries an Aco2-clade paralog (E8, E9). Pointing to fungi: it respires, and the abstract links respiration to the single-aconitase pattern. Neither is evidence of activity. |

**Genes** (AGENTS.md §4):
- `YALI1E03387g`: UniProt-reviewed name LYS4, assigned by similarity (no *Y. lipolytica* experiment cited).
  Homoaconitase, cis-homoaconitate ⇌ homoisocitrate. Homology plus *S. cerevisiae* / *A. fumigatus* biochemistry.
  Consensus-essential.
- `YALI1D11984g`: no verified name. ACO1-type aconitate hydratase. Homology plus the yeast and fungal biochemistry
  above. Consensus-essential.
- `YALI1E17935g`: no verified name. ACO2-type aconitase-family protein; active in *S. cerevisiae*, inactive in
  *Aspergillus*. Homology only. Not on the list.
- No gene here has *Y. lipolytica* experimental validation; a GPR assignment is not one.

### 1b. Name, annotation, direction, protons

| # | Claim | Tier |
| --- | --- | --- |
| N1 | R24's cross-references all describe **2-methylcitrate dehydratase**: Rhea 57496 = "2-methylcitrate = 2-methyl-cis-aconitate + H2O"; MNXR145545 → rh:57496, EC 4.2.1.79; BiGG MCITD. Its `ec-code` 4.2.1.36 is the *homoisocitrate* step (KEGG EC 4.2.1.36 → R04371). Its note `PROTEIN_CLASS` 4.2.1.114 is "methanogen homoaconitase" (both steps in one enzyme). | verified ([rhea_selected.tsv](evidence/rhea_selected.tsv), [mnx_MNXR145545.html](evidence/mnx_MNXR145545.html), [kegg_ec4.2.1.36.txt](evidence/kegg_ec4.2.1.36.txt), [kegg_ec4.2.1.114.txt](evidence/kegg_ec4.2.1.114.txt)) |
| N2 | The correct identifiers for R24's chemistry are KEGG [R03444](https://www.kegg.jp/entry/R03444) (C01251 ⇌ C04002 + H2O; EC 4.2.1.114, 4.2.1.-) and [RHEA:26101](https://www.rhea-db.org/rhea/26101) ((2R)-homocitrate = cis-homoaconitate + H2O; no EC). | verified (KEGG REST, Rhea TSV) |
| N3 | **m76** is named "but-1-ene-1,2,4-tricarboxylic acid" (= cis-homoaconitate), but every cross-reference is **2-methyl-cis-aconitate**: KEGG C04225 "(Z)-But-2-ene-1,2,3-tricarboxylate", ChEBI 57872, MNXM1792 "2-methyl-cis-aconitate", BiGG 2mcacn, plus the model's own InChI with a methyl branch. cis-Homoaconitate is KEGG C04002 / ChEBI 58174. Both are C7H5O6, −3, so balance is unaffected. This is a local mis-annotation (the databases are consistent with each other), not a database error. | verified (KEGG, ChEBI API, MetaNetX; model dump) |
| N4 | Where the name came from: Yeast-GEM `r_0027` "2-methylcitrate dehydratase" is homocitrate ⇌ but-1-ene-1,2,4-tricarboxylate + H2O on YDR234W (LYS4). iMM904 `MCITDm` is the same chemistry with no gene. eciYali5 reaction 27 "homoaconitase" puts YALI0E02728g on homocitrate → homoaconitate. Yeast-GEM `r_0542` (homoisocitrate step) reads `YJL200C or YDR234W`, putting ACO2 on the step Fazius assigns to Lys4. R24's rule and name look inherited along this line. | verified that these external models contain these reactions ([logs/external_models_hits.txt](logs/external_models_hits.txt), SHAs in [logs/external_models_sha256.txt](logs/external_models_sha256.txt)); the inheritance itself is inferred |
| N5 | **3 H+.** The source `iyali26.xml` writes R24 as m75 ⇌ H2O + m76, with no formulas and every charge 0. In `e1e7c868`, m75 and m624 are neutral C7H10O7 while m76 is the trianion C7H5O6(−3), so R24 releases 3 H+ and R395 takes them back. R1170 lost and R398 gained one H+. ChEBI has homocitrate(3−) 58884, cis-homoaconitate(3−) 58174 and homoisocitrate(3−) 15404 all as C7H7O7 / C7H5O6 / C7H7O7 at −3, so the physiological R24 has **no** H+. | source vs candidate: verified ([logs/static.log](logs/static.log)). ChEBI: verified. That the build's charge-aware proton/water balancer (`metabolites.fix_proton_water_balance`) added the H+: **inferred** from its docstring; not traced in a build log. |
| N6 | The 3 H+ change no FBA result. m76 has only R24 and R395, so at steady state v24 = v395 and the H+ cancel. | statically verified (metabolite reaction sets in [logs/static.log](logs/static.log)) |
| N7 | **Direction.** R1170 (only producer of m75) and R398 (only consumer of m624) are irreversible, and m75, m76, m624 each have exactly two reactions. Any steady state therefore has v1170 = v24 = v395 = v398 ≥ 0, so R24 can carry only forward net flux; the reversible bounds are harmless. Biochemically the step is reversible and uphill ("the equilibrium of the reaction is far on the side of homocitrate", Fazius 2012, Results), pulled by the decarboxylation. KEGG and Rhea write it reversibly. | network: statically verified; biochemistry: verified (full text) |
| N8 | **Compartment.** Mitochondrion is consistent: UniProt puts Q6C791, Q6C9P6 and Q6C5V2 there (by similarity or automatic), and Fazius removed mitochondrial leaders before expression. Not reviewed further. | verified that UniProt states it; not experimental for *Y. lipolytica* |

## 2. Knockouts of each option (in-memory, verified)

[scripts/r24_variants.py](scripts/r24_variants.py), run by [scripts/run_variants.zsh](scripts/run_variants.zsh) on
a SHA-checked short copy of `e1e7c868`.
- Inputs: Task A medium (`6eca70fd…`) + PO1f (`35307853…`); Gurobi at 1e-9, reset per solve, targeted check.
- Each variant changes only R24's rule; R395 stays `YALI1E03387g`.
- Output in [variants.json](variants.json) and [logs/variants.log](logs/variants.log). 25 LPs, all optimal, all
  passed the check.
- WT 0.3116608216 in every variant.

| Option | R24 rule | `YALI1E03387g` KO/WT | `YALI1D11984g` | `YALI1E17935g` | double `D11984g`+`E17935g` |
| --- | --- | --- | --- | --- | --- |
| A status quo | `YALI1E03387g` | 0.0694 | 0.0244 | 1.0 | 0.0244 |
| B ACO1-type only | `YALI1D11984g` | 0.0694 | 0.0244 | 1.0 | 0.0244 |
| C both aconitases | `YALI1D11984g or YALI1E17935g` | 0.0694 | 0.0244 | 1.0 | 0.0244 |
| D ACO2-type only | `YALI1E17935g` | 0.0694 | 0.0244 | **0.0694** | 0.0244 |
| E keep LYS4 and add both | `YALI1E03387g or YALI1D11984g or YALI1E17935g` | 0.0694 | 0.0244 | 1.0 | 0.0244 |

**Classes** (essential when KO/WT < cutoff; 1 / 5 / 10 / 15 %):
- `YALI1E03387g` (listed): FN / FN / TP / TP in every option. R395 still closes.
- `YALI1D11984g` (listed): FN / TP / TP / TP in every option.
- `YALI1E17935g` (not listed): TN everywhere, except option D: TN / TN / **FP / FP**.

**Why nothing else moves.**
- The `YALI1E03387g` knockout closes R395 in every option, and growth is capped by the medium's lysine bound. Lysine
  uptake R1220 sits at its bound, 0.00595.
- The `YALI1D11984g` knockout already grows only 0.0076 (set by the medium's amino-acid bounds, inferred in the aconitase round §6). At that rate the
  medium's lysine covers demand (uptake 0.0021 < 0.00595), so also closing R24 under option B costs nothing.
- Only genes in R24's old or new rule can change class, because no other rule, bound or stoichiometry changes.
  Statically verified.

**Scores against `e1e7c868`.**
- Options A, B, C and E: no change at any cutoff.
- Option D: FP +1, TN −1 at 10 % and 15 %.
  - 10 %: 108 / 207 / 43 / 691. Specificity 0.9428 → 0.9414, precision 0.7200 → 0.7152, MCC 0.3740 → 0.3712.
  - 15 %: 125 / 190 / 43 / 691. Precision 0.7485 → 0.7440, MCC 0.4255 → 0.4228.
- Baseline counts come from the aconitase round's `ko_compare_summary.json` (`e0123ada…`), recomputed this session.
- Common to every score: negatives are the user-defined class (model genes not on the list); the 16 mtDNA genes are
  out of scope; 16 unresolved; 1,286 listed positives are out of scope. The list was used in development, so these
  are development scores, not validation.
- Match rate was not used to rank the options.

## 3. What the results mean

- **Mechanism.** On this medium the R24 rule cannot move `YALI1E03387g` or `YALI1D11984g`: the lysine supply in
  Task A and the TCA block dominate. The `YALI1E03387g` 1 % / 5 % FN is a lysine-bound effect, not a GPR effect.
  (Inferred from the uptake sitting at its bound; no bound was varied.)
- **Biology.** The literature points to an aconitase, not Lys4, for R24. Which aconitase is the open question. The
  only discriminating evidence would be *Y. lipolytica*-specific: an *aco2*-type null on lysine-free medium, or
  activity of recombinant Q6C5V2 with homoaconitate. Homology cannot settle it: E9 shows the Aco2-clade K is shared by
  active and inactive members.

## 4. Could the method be improved, and how to test

- A lysine-free medium would separate options B and C (the `YALI1D11984g` knockout under B would also lose lysine
  synthesis). That is a new condition and needs its own authorization; it is not benchmark evidence.
- A gene tree including more Saccharomycotina Aco2 homologues (active in *S. cerevisiae*) and Pezizomycotina AcoB
  (inactive) would show which side the *Y. lipolytica* paralog falls on. That is still homology.
- The proton/water balancer mixes a neutral homocitrate with a trianionic homoaconitate. A scan for other chains with
  mixed protonation of one carbon skeleton would find similar artefacts (see the user's memory note on protonation
  conventions). Not done.

## 5. Options for the user

Each is a separate `curation-change` round. The 3 H+ fix (option H) is independent of the GPR choice.

| Option | Change | For | Against | Acceptance criteria |
| --- | --- | --- | --- | --- |
| **A** keep | none | no new claim | R24 carries a gene whose only reviewed activity is the other step (E1, E7) | n/a |
| **B** ACO1-type only | R24 `YALI1E03387g` → `YALI1D11984g` | E2, E5, E6; Aco1p hydrates homoaconitate (E3); respiring-fungus pattern | drops the ACO2-type paralog that is active in *S. cerevisiae* (E4) | diff = R24 rule + notes + gene links only; control byte-identical; KO classes as §2 row B (no change); full pytest only new tests differ |
| **C** both aconitases (my recommendation) | R24 → `YALI1D11984g or YALI1E17935g` | the yeast pattern, where each enzyme alone supports lysine synthesis (E4); *Y. lipolytica* has both clades (E8) | `YALI1E17935g` activity is unproven (E9 and E5 say Aco2-clade members can be inactive); the OR member stands on transfer from *S. cerevisiae* only, to be labelled as such in notes | as B; KO classes as §2 row C (no change); notes state the tier "homology + *S. cerevisiae* biochemistry, inactive homologues exist in *Aspergillus*" |
| **D** ACO2-type only | R24 → `YALI1E17935g` | "Aco2p specifically and exclusively contributes to lysine biosynthesis" (abstract) | contradicts E4: Aco1p also contributes in yeast | as B; `YALI1E17935g` TN → FP at 10/15 % must be reported, not tuned away |
| **E** keep LYS4, add both | R24 → `YALI1E03387g or …` | none | keeps an OR member with no evidence for the step (E1) | not recommended |
| **H** chemistry (any option) | m75, m624 (and nuclear m628) → C7H7O7, −3; remove the 3 H+ from R24 and R395; rebalance R1170, R396, R398 | matches ChEBI (N5) | charge and formula edits ripple to every reaction using these metabolites; FBA-neutral for this chain (N6) | every touched reaction mass- and charge-balanced; model-wide balance re-check; WT and all KO growths identical |
| **N** annotation (any option) | R24 name → homocitrate dehydratase (or "homoaconitase, aconitase-catalysed step"); R24 refs → KEGG R03444, RHEA:26101, EC 4.2.1.-; drop BiGG MCITD, Rhea 57496, MNXR145545, `ec-code` 4.2.1.36; review `PROTEIN_CLASS` 4.2.1.114. m76 refs → KEGG C04002, ChEBI 58174 (MetaNetX ID to look up) | N1–N3 | annotation only; no FBA effect | whole-model diff shows annotation and name fields only |

Why C over B: on this medium they give the same knockouts, so the choice rests on evidence alone. C follows the only
organism with genetics for both clades (E4). B is the more conservative choice if an unproven OR member is
unacceptable under AGENTS.md §5. Either is defensible; the user decides.

## 6. Noticed in passing (not reviewed)

- Nuclear homocitrate m628 has a single reaction (R396, nuclear homocitrate synthase), so R396 is a dead end
  (verified by m628's reaction set).
- R95 "aconitase" (2-methylcitrate → m217) carries the citrate aconitase's cross-references (KEGG R01325, BiGG
  ACONTa). Its product m217's InChI (`c1-14-7…`) is an O-methyl ether, not 2-methylisocitrate. Seen in the dump only.
- R238 `PROTEIN_CLASS` 4.2.1.4 (inherited; flagged by the aconitase round).

## 7. Limits

- One medium, one profile. No lysine-free test.
- The *Y. lipolytica* aconitase activities on homocitrate and homoaconitate are unknown (E10, E11).
- Fazius' yeast coupled assays were "not shown"; full text was read as Europe PMC XML, not the PDF.
- The assignment of the H+ addition to a specific build step is inferred (N5).
- The baseline is an uncommitted candidate on another branch. If the aconitase round changes before it merges, these
  results must be re-checked against the merged SHA.
