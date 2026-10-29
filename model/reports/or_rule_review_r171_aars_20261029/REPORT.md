# OR-rule review: R171 ATP synthase and the cytosolic aminoacyl-tRNA synthetase rules (2026-10-29)

Scope and authorization: [TASK.md](TASK.md). Investigation only: rules were edited in memory and restored; no model
XML, GPR, curated data, medium, benchmark or platform code changed; nothing committed or pushed. Per-gene table:
[or_rule_review_genes.csv](or_rule_review_genes.csv); raw results: [variant_knockouts.json](variant_knockouts.json),
[solver_calls.jsonl](solver_calls.jsonl); scripts: [scripts/](scripts/); fetched evidence: [evidence/](evidence/).

## 1. Summary

- **Both rule families are inherited, not introduced here** (verified). The five rule strings in the candidate are
  byte-identical to those in `model/source/iyli21.xml` (`6974b758…`).
- **R171 (ATP synthase).** The rule is `(14-subunit AND) or YALI1B28025g or YALI1B15733g`. The two single-gene OR
  branches are subunit g (ATP20) and subunit k (ATP19), 138 and 70 residues (UniProt, verified). An OR branch asserts
  that the gene alone performs the reaction; a single accessory subunit cannot synthesise ATP, so the rule is wrong
  by construction. Mechanism (verified): every single knockout of the 16 genes leaves R171 open (KO/WT 1.0000);
  closing R171 gives KO/WT 0.1097, below the 15 % cutoff and above 10 %.
- **What the fix would do.** Removing the two OR branches turns the 9 consensus-essential core subunits from FN to
  TP at 15 % and nothing at 1 / 5 / 10 % (the 0.1097 energy ceiling is the model's fermentative escape, out of scope
  here). How subunits e (`YALI1E38131g`), j (`YALI0D17490g`), g and k appear in the rule decides the false positives:
  0 to 4 of them. The *Y. lipolytica* screens call e, g and k non-essential; j was not screened. In the
  *Y. lipolytica* cryo-EM work (Hahn 2016) e, g and k are dimer-specific and the monomer lacking them hydrolyses ATP
  at the dimer's rate; in *S. cerevisiae* (Arnold 1998) they are not required for an enzymatically active ATP
  synthase; the *S. cerevisiae* subunit-j null (Vaillier 1999) is viable with half the ATPase activity and a low
  ADP/O ratio. None of this is a knockout result in *Y. lipolytica*.
- **Cytosolic synthetase rules.** `R759` (Thr; `R760` in the mitochondrion carries the same rule), `R723` (Ser) and
  `R673` (Pro) each read `cytosolic-type gene or mitochondrial-type paralog`. All six genes are consensus-essential and
  all six are FN. Closing any of the three cytosolic reactions gives growth 0 (verified). Assigning each cytosolic
  reaction to its cytosolic gene turns 3 genes FN → TP at every cutoff with no new FP; the three mitochondrial-type
  paralogs stay FN because `R760` is blocked and no mitochondrial Ser or Pro charging reaction exists (category d
  of the 2026-10-27 study, not a rule problem).
- **Localisation evidence** (section 7.2): the Thr and Ser pairs are supported by UniProt names and locations, by
  length and UniRef50 clustering against the yeast orthologs, and by the N-termini; the Pro pair rests on the
  genome annotation's ancestral-locus (synteny) assignment and an N-terminal heuristic. No *Y. lipolytica*
  experimental localisation was found for any of the six. A formal presequence predictor was not run.
- **Recommendation** (candidates, not changes; section 8): (1) the synthetase split for Thr and Ser now, Pro after a
  predictor run; (2) R171 with the OR branches removed and a decision on e and j taken on the catalytic evidence,
  never on the FP count. Separate follow-ups: the `YALI0D17490g` identifier (a CLIB122 id inside a YALI1 model; the
  S2 map gives `YALI1_D21476g`), and the 0.1097 energy ceiling.

## 2. Run and inputs (verified)

| Item | Value |
| --- | --- |
| Model | `candidates/E5_…_dolichol_route_v2_r2171_1to1_20261028.xml`, SHA-256 `965ab947f2aadc2c69443e7a794358c502acfdf0892e78582278ff07a37429d7` (PR #31). Not a designated reference. `5f54eb22…` (PR #30) has the same gene set and identical `R171`, `R759`, `R760`, `R723`, `R673` (checked); not run. |
| Medium / profile | Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile `po1f_sd_leu.json` (`35307853…`); `PO1f_plasmid_LEU2` excluded at runtime |
| List / scope rule | `consensus_essential_genes.csv` (`1e887f5a…`), 1,612 positives, 326 in the model; mtDNA rule `16b3e9b4…` (16 genes out of scope) |
| Solver | Gurobi 13.0.1 (worktree venv, Python 3.13.5, cobra 0.30.0), FeasibilityTol / OptimalityTol 1e-9, Presolve 0, 1 thread, `problem.reset()` before every solve; a point is scorable only if optimal, \|S v\| ≤ 1e-6, bound violation ≤ 1e-6 and knocked flux ≤ 1e-9 |
| WT | 0.3116608216 h⁻¹ (residual 5e-14), equal to the baseline of the R2171 and UM7 sweeps |
| LPs | 1,286 in the complete run (run 3, 62.7 s): WT, 6 reaction-level knockouts, 9 variants × 22 genes, 1,081-gene sweep. Across the three attempts recorded in TASK.md §2: 0 + 344 + 1,286 = 1,630 |
| Checks | Every variant knockout point scorable (0 failing); the full sweep's failing points and unresolved knockouts are listed in section 6 |
| Negatives | User-defined class (model genes not on the list); the list was used during development, so all scores are a development reference |

## 3. The rules as inherited

| Reaction | Compartment | Rule in `965ab947` (= iYli21, verified) |
| --- | --- | --- |
| `R171` ATP synthase, EC 7.1.2.2 | mitochondrion (4 H⁺ cytosol → matrix per ATP) | `(YALI1M00260g and YALI1M00267g and YALI1M00318g and YALI1F04433g and YALI1F04110g and YALI1D27826g and YALI1D15734g and YALI1D14639g and YALI1B08885g and YALI1F07384g and YALI1B05364g and YALI0D17490g and YALI1F26967g and YALI1E38131g) or YALI1B28025g or YALI1B15733g` |
| `R759` threonyl-tRNA synthetase | cytosol | `YALI1A09121g or YALI1F14882g` |
| `R760` threonyl-tRNA synthetase | mitochondrion | `YALI1A09121g or YALI1F14882g` |
| `R723` seryl-tRNA synthetase | cytosol | `YALI1B07957g or YALI1F03908g` |
| `R673` prolyl-tRNA synthetase | cytosol | `YALI1D19110g or YALI1E05929g` |

No other ATP synthase reaction exists (EC 7.1.2.2 / 3.6.3.14 search: only `R171`; `R_NTP1` is the balanced NTP
hydrolysis). The mitochondrial synthetase reactions `R760`, and the Ser / Pro mitochondrial steps (absent), are not
touched by this review.

**Convention in Yeast-GEM 9.0.2** (`reference/external_models/yeast_gem/Yeast-GEM.xml`, read this session; a
convention, never a reason, AGENTS §5): `r_0226` ATP synthase is a 17-gene AND that includes e (YDR322C-A), g
(YPR020W), k (YOL077W-A) and j (YML081C-A); `r_1042` cytosolic Thr = THS1 alone, `r_1043` mitochondrial Thr = MST1
alone; `r_0995` cytosolic Ser = `SES1 or DIA4` (the same defect as `R723`); `r_0941` cytosolic Pro = YHR020W alone,
`r_4187` mitochondrial Pro = AIM10 alone.

## 4. The genes: identity and experimental calls (verified)

Identities come from UniProt (CLIB122 `YALI0` entries, mapped to the model's `YALI1` ids through
`reference/locus_map/s2_metabolic_genes.csv`; TSV in [evidence/](evidence/README.md)). Screen calls are the Cas9 and
Cas12a rows of `cas9_cas12a_fitness.csv` (fitness score FS, q value); the consensus list is "essential in at least 2 of
3 screens", and the third screen is not in that normalized table. Gene names are homolog names; none is a verified
*Y. lipolytica* name.

| Gene | Protein (UniProt) | aa | Consensus | Cas9 | Cas12a |
| --- | --- | ---: | --- | --- | --- |
| `YALI1B28025g` (YALI0_B21527g) | ATP synthase subunit g, ATP20 ([B5FVB8](https://www.uniprot.org/uniprotkb/B5FVB8/entry)) | 138 | no | non-essential (FS −2.80, q 0.085) | non-essential (−1.51, q 0.32) |
| `YALI1B15733g` (YALI0_B11913g) | ATP synthase subunit k, ATP19 ([B5FVB3](https://www.uniprot.org/uniprotkb/B5FVB3/entry)) | 70 | no | non-essential (−0.27, q 1) | non-essential (−0.63, q 0.91) |
| `YALI1E38131g` (YALI0_E32164g) | ATP synthase subunit e, TIM11/ATP21 ([B5FVG3](https://www.uniprot.org/uniprotkb/B5FVG3/entry)) | 90 | no | non-essential (−1.24, q 1) | not in table |
| `YALI0D17490g` (= YALI1_D21476g) | ATP synthase subunit j, ATP18 ([Q6C8S0](https://www.uniprot.org/uniprotkb/Q6C8S0/entry)) | 62 | no (id never matches) | not in table | not in table |
| 9 core subunits (α ATP1, β ATP2, γ ATP3, δ ATP16, OSCP ATP5, b ATP4, d ATP7, h ATP14, f ATP17) | see CSV | | all yes | all essential | 5 essential, 4 non-essential |
| 3 mtDNA (`YALI1M00260g`, `M00267g`, `M00318g`) | ATP6/8/9 set | | out of scope | not in table | not in table |
| `YALI1F14882g` (YALI0_F11231g) | Threonine–tRNA ligase, cytoplasmic ([Q6C234](https://www.uniprot.org/uniprotkb/Q6C234/entry)) | 755 | yes | essential (−3.09, q 0.021) | essential (−2.35, q 6e-4) |
| `YALI1A09121g` (YALI0_A09427g) | threonine–tRNA ligase, "Mitochondrion matrix" (ARBA) ([Q6CHF9](https://www.uniprot.org/uniprotkb/Q6CHF9/entry)) | 461 | yes | essential (−4.73) | essential (−3.67) |
| `YALI1F03908g` (YALI0_F02629g) | Serine–tRNA ligase, cytoplasmic ([Q6C350](https://www.uniprot.org/uniprotkb/Q6C350/entry)) | 453 | yes | essential (−3.07, q 0.023) | essential (−2.97) |
| `YALI1B07957g` (YALI0_B05918g) | serine–tRNA ligase ([Q6CFL2](https://www.uniprot.org/uniprotkb/Q6CFL2/entry)) | 470 | yes | essential (−4.49) | non-essential (−0.78, q 1) |
| `YALI1E05929g` (YALI0_E05027g) | proline–tRNA ligase ([Q6C6Z6](https://www.uniprot.org/uniprotkb/Q6C6Z6/entry)) | 517 | yes | essential (−4.11) | essential (−2.58) |
| `YALI1D19110g` (YALI0_D15620g) | proline–tRNA ligase ([Q6C8Z3](https://www.uniprot.org/uniprotkb/Q6C8Z3/entry)) | 445 | yes | essential (−4.94) | essential (−3.76) |

`YALI0D17490g` is a CLIB122 identifier inside a model whose genes and benchmark list use `YALI1` ids. The S2 locus map
pairs it with `YALI1_D21476g` (NCBI gene 2911301, KEGG map 00190), marked `in_model: no`, and `YALI1D21476g` is absent
from the consensus list and from both screen tables. Under the benchmark's label rule the model gene is a negative;
it was never screened.

## 5. Mechanism (verified, this run)

- **Baseline single knockouts** of all 22 genes: KO/WT 1.0000, and none closes a reaction (`closed_reactions` empty
  in the JSON): for every AND member the `or g or k` branches keep `R171` open; for every synthetase gene the OR
  partner keeps the reaction open.
- **Reaction-level knockouts** (original rules): `R171` closed → KO/WT **0.10970**; `R759` → **0.0000**; `R760` →
  1.0000 (the mitochondrial threonyl step carries no flux at the optimum: dead-end mitochondrial tRNA, 2026-10-27
  study); `R759`+`R760` → 0.0000; `R723` → **0.0000**; `R673` → **0.0000**.
- So the rules, not the network, hide 9 ATP synthase and 3 cytosolic synthetase knockouts. The three
  mitochondrial-type paralogs are hidden by network disconnection and would stay FN under any rule.
- The 0.1097 without ATP synthase is the substrate-level ATP the capped glucose (2.414) can supply after NGAM; it
  is the same energy ceiling as the complex III/IV infeasibility of 2026-10-21. Biologically *Y. lipolytica* cannot
  live on fermentation (established, uncited here), so a correct rule would still under-predict the severity at the
  10 % cutoff and below.

## 6. Rule variants scored in memory (verified)

Only the 22 genes of the five rules can change class under a rule edit (a knockout changes bounds only on reactions
whose rule contains the gene; provable from the knockout semantics). Transitions are relative to the in-run
baseline; the "after" matrices add them to the recorded baseline of `965ab947`
(`reports/r2171_pi_plc_1to1_20261028/ko_compare_summary.json`: 1 % 94 / 221 / 35 / 699, 5 % 100 / 215 / 37 / 697,
10 % 104 / 211 / 42 / 692, 15 % 112 / 203 / 42 / 692; 16 unresolved at each cutoff).

| Variant | Rule change | 15 %: transitions | 15 %: TP / FN / FP / TN, recall, MCC | 1 / 5 / 10 % |
| --- | --- | --- | --- | --- |
| baseline | none | — | 112 / 203 / 42 / 692, 0.356, 0.386 | as recorded |
| R171_A | drop `or g or k` (14-AND kept) | 9 FN→TP, 2 TN→FP (e, j) | 121 / 194 / 44 / 690, 0.384, 0.408 | no change |
| R171_B | 16-AND (g, k as AND members) | 9 FN→TP, 4 TN→FP (e, j, g, k) | 121 / 194 / 46 / 688, 0.384, 0.403 | no change |
| R171_C | 13-AND (also drop e) | 9 FN→TP, 1 TN→FP (j) | 121 / 194 / 43 / 691, 0.384, 0.411 | no change |
| R171_D | 12-AND (also drop e and j) | 9 FN→TP | 121 / 194 / 42 / 692, 0.384, 0.414 | no change |
| aaRS_split | `R759`←F14882, `R760`←A09121, `R723`←F03908, `R673`←E05929 | 3 FN→TP | 115 / 200 / 42 / 692, 0.365, 0.396 | 3 FN→TP at each |
| aaRS_split_alt | as above but `R673`←D19110 | 3 FN→TP (D19110 instead of E05929) | same counts | same counts |
| combined A + split | | 12 FN→TP, 2 TN→FP | 124 / 191 / 44 / 690, 0.394, 0.417 | 3 FN→TP at each |
| combined D + split | | 12 FN→TP | 124 / 191 / 42 / 692, 0.394, 0.423 | 3 FN→TP at each |

Every R171 variant gives KO/WT 0.1097 for each AND member, so the gain exists only at 15 %. The synthetase split gives
0.0000 for the cytosolic gene at every cutoff. The alternative Pro assignment changes which gene flips, not the count.

**Confirming full sweep** (all 1,081 genes under `combined_A_split`, 1,081 LPs, no failing point, max knocked flux 0):
1 % 97 / 218 / 35 / 699, 5 % 103 / 212 / 37 / 697, 10 % 107 / 208 / 42 / 692, 15 % 124 / 191 / 44 / 690 (recall 0.394,
precision 0.738, MCC 0.417), 16 unresolved at every cutoff (the same infeasible complex III / IV knockouts), 16 mtDNA
genes out of scope. This equals the recorded baseline plus the predicted transitions at every cutoff, and the number
of model-essential genes outside the 22 target genes at 15 % is 154, the baseline's TP + FP. The variant knockouts
and reaction-level values of run 3 are identical to those of the stopped run 2.

## 7. Evidence for and against each candidate edit

### 7.1 R171

Sources opened this session: Hahn et al. 2016, *Mol Cell*, PMID [27373333](https://europepmc.org/article/MED/27373333)
(full text through the Europe PMC REST service; *Y. lipolytica* dimer structure); Liu et al. 2015, *Biochem J*, PMID
[25759169](https://europepmc.org/article/MED/25759169) (full text; purified fungal ATP synthases including
*Y. lipolytica*); Arnold et al. 1998, *EMBO J*, PMID [9857174](https://europepmc.org/article/MED/9857174) (abstract
record); Vaillier et al. 1999, *J Biol Chem*, PMID [9867878](https://europepmc.org/article/MED/9867878) (abstract
record). The fetch tool returned paraphrases plus short verbatim quotes; the quotes below are the verbatim parts.

| Subunit (gene) | For removal from the AND / against OR | Against removal | Level |
| --- | --- | --- | --- |
| g, ATP20 (`YALI1B28025g`) | Screens: non-essential in both (Cas9 q 0.085). Hahn 2016, *Y. lipolytica*: e, g, k "known as dimer specific"; the DDM monomer lacking e, g, k hydrolyses ATP at about the dimer's rate (§2.1); e and g "do not participate directly in dimer formation" but induce the ~100° membrane bend. Arnold 1998, *S. cerevisiae*: "these three proteins are not required for the formation of enzymatically active ATP synthase"; "Su e/Tim11 and Su g are essential for the formation of the dimeric state." | Dimer loss affects cristae and respiratory growth in *S. cerevisiae* (fitness, not lethality; recalled). Cas9 FS −2.80 suggests a fitness cost. Hahn 2016 reports hydrolysis, not synthesis. | verified (sources above); the "not a catalytic requirement" conclusion is inferred from them |
| k, ATP19 (`YALI1B15733g`) | Screens: non-essential in both. Arnold 1998 as above. Liu 2015: in the purified *Y. lipolytica* enzyme "Subunits g, k and l were not found" (Table 2 footnote); roles of g, j, k "remain to be elucidated". Hahn 2016: k only tentatively placed. | none found | verified / inferred as above |
| e, TIM11 (`YALI1E38131g`), now an AND member | Screens: Cas9 non-essential (q 1); no Cas12a row. Same literature as g. | Liu 2015 lists e among the *Y. lipolytica* subunits (Table 2) though the DDM preparation "lacked subunit e". | as above |
| j, ATP18 (`YALI0D17490g`), now an AND member | Hahn 2016 calls subunit i "the yeast-specific, non-essential subunit i" (§2.5.5). Vaillier 1999, *S. cerevisiae* null: "low growth on nonfermentable medium", "2-fold less specific ATPase activity", "low ADP/O ratio". Not screened in *Y. lipolytica*. | The same null data show it is needed for full coupling and respiratory growth; UniProt (Q6C8S0) places it in the *Y. lipolytica* enzyme. | verified sources; the keep-or-drop call is a judgement |
| 9 core subunits | all consensus-essential; AND membership consistent with the screens | — | verified |

The removal of the two OR branches needs no biology: a 70- or 138-residue accessory subunit does not catalyse ATP
synthesis on its own (AGENTS §5: an OR member needs evidence that it catalyses the step alone; there is none). The
biology decides only whether g, k, e and j become AND members or leave the rule. No *Y. lipolytica* knockout of any
of the four exists in the sources opened; the screens are the only organism-specific phenotype data.

### 7.2 Cytosolic synthetase pairs

| Pair | Cytosolic candidate: evidence | Mitochondrial-type candidate: evidence | Level |
| --- | --- | --- | --- |
| Thr, `R759` | `YALI1F14882g` (Q6C234): UniProt name "cytoplasmic", location Cytoplasm (ARBA); 755 aa vs THS1 734; annotation "similar to THS1"; N-terminus acidic (net −1, 7 D/E in 30). | `YALI1A09121g` (Q6CHF9): UniProt location "Mitochondrion matrix" (ARBA, automatic); 461 aa vs MST1 462; annotation "similar to MST1 … Mitochondrial threonyl tRNA synthase" (its cited accession Q6H9P8 is deleted in UniProt); N-terminus `MLRTALRRAAKPSRVPRRTLTTS…` (8 R/K, 0 D/E). | verified (UniProt), heuristic (N-terminus) |
| Ser, `R723` | `YALI1F03908g` (Q6C350): UniProt "cytoplasmic"; member of [UniRef50_P07284](https://www.uniprot.org/uniref/UniRef50_P07284) with SES1 (verified); 453 aa vs SES1 462; annotation "highly similar to SES1"; N-terminus acidic (net −2). | `YALI1B07957g` (Q6CFL2): own 2-member UniRef50 cluster; 470 aa vs DIA4 446; annotation "similar to DIA4" (DIA4 = mitochondrial SerRS, UniProt P38705, verified); N-terminus `MLRRITRAVPRVNARLFSSTSVS…` (6 R/K, 2 D/E at positions 27–28). | verified / heuristic |
| Pro, `R673` | `YALI1E05929g` (Q6C6Z6): annotation "similar to YHR020W … ancestral locus Anc_1.352" (YHR020W = the cytoplasmic ProRS, UniProt P38708, 688 aa; no mitochondrial annotation); 517 aa; UniRef50 cluster A0A1Y2FP29 (27 eukaryotes), not YHR020W's; N-terminus `MPPKAAKAKNDNMLIGITVNKEDDFSG…` (4 D/E in 30, net 0). | `YALI1D19110g` (Q6C8Z3): annotation "similar to C. albicans ProRS cytoplasmic, similar to YER087W … ancestral locus Anc_7.368" (YER087W = AIM10, mitochondrial ProRS, UniProt P39965, verified); 445 aa vs AIM10 576; UniRef50 singleton; N-terminus `MLPRLSRLTTRRLKDVTPSTTASD…` (5 R/K, 2 D/E, 9 S/T). | annotation + synteny (recalled from the annotation strings) + heuristic; weakest pair |

The N-terminal numbers come from `scripts/nterm_heuristic.py` and are a rough indicator, not a prediction.
TargetP 2.0 / DeepLoc 2.1 (the services used in `reports/glypro_localization_predictors_20260924`) need a web form
submission and were not run. Europe PMC searches for *Y. lipolytica* synthetase localisation (99 and 16 hits; two
queries, top 20 and 15 read) returned nothing on these enzymes.

Both Pro assignments give the same scores (section 6). The reason to get it right is the model, not the score: the
wrong assignment would make the mitochondrial enzyme the cytosolic one.

## 8. Candidate edits, acceptance criteria, and what is not a reason

Each candidate must go through `curation-change` with the user's scoped authorization. Implementation path as for
`R2193` (`curation/gpr/r2193_iron_gpr.json`, `platform/scripts/gem_annotate/r2193_iron_gpr.py`): one curated JSON
per reaction (before / after rule, evidence), a guarded builder step behind an opt-in flag, a focused test, a build
to a **new** candidate on top of the newest candidate, a whole-model diff, an independent audit, REPORT and STATE.

| # | Candidate | Expected result (verified in memory) | Acceptance criteria |
| ---: | --- | --- | --- |
| 1 | Synthetase split, Thr and Ser: `R759` ← `YALI1F14882g`, `R760` ← `YALI1A09121g`, `R723` ← `YALI1F03908g`; `YALI1B07957g` kept in the model unlinked (precedent: `YALI1F07747g`, 2026-10-22) | +2 TP at every cutoff, no other class change | Whole-model diff = the three rules and notes; gene list unchanged; WT unchanged; full sweep = baseline + exactly those transitions; 16 unresolved unchanged; no failing point |
| 1b | Pro: `R673` ← `YALI1E05929g`, `YALI1D19110g` unlinked | +1 TP at every cutoff | Same, after a presequence prediction (TargetP 2.0 / DeepLoc 2.1) agrees with the synteny assignment; if it disagrees, stop and report |
| 2 | `R171`: remove `or YALI1B28025g or YALI1B15733g`; decide e and j: option A keeps both as AND members (+2 FP), C drops e (+1 FP), D drops e and j (0 FP) | +9 TP at 15 % only | Whole-model diff = `R171` rule and notes; the option is stated with its catalytic ground (section 7.1), the FP count is reported as a consequence; full sweep = baseline + the stated transitions; 16 unresolved unchanged |

Not reasons for any choice: the recall or MCC gain; Yeast-GEM's 17-gene AND; the number of false positives an option
avoids. Separate lines, not part of these edits: remapping `YALI0D17490g` to `YALI1D21476g` (annotation); the
fermentative energy ceiling (0.1097 without ATP synthase; NGAM / glucose cap; 2026-10-21 entry); mitochondrial Ser
and Pro charging reactions and the dead-end mitochondrial tRNAs (category d).

## 9. Verification levels and limits

- **Verified this session:** rule strings and their iYli21 origin; identities, lengths, UniProt locations and
  UniRef50 memberships (REST, TSV saved with SHA); screen calls; reaction-level and single-gene knockouts; the
  variant transitions; the literature statements quoted (through the Europe PMC REST records and full texts, with
  the fetch tool's paraphrase around the quotes).
- **Inferred:** that dimer-specific subunits are not catalytic requirements (from the sources, not from a
  *Y. lipolytica* knockout); the Pro localisation; that an edit would flip exactly the listed genes on a built
  candidate (checked in memory, and by the full sweep for one combination).
- **Recalled, not re-checked:** that *S. cerevisiae* dimer mutants have respiratory fitness defects; that
  *Y. lipolytica* cannot grow by fermentation.
- **Not done:** a presequence predictor; any medium other than Task A; the third screen's per-gene calls; an
  independent audit; pytest (no platform code changed). The two mitochondrial Pro / Ser reactions do not exist, so the
  mitochondrial paralogs' essentiality cannot be modelled by a rule change.
- The LP spend across three attempts is in TASK.md §2; each attempt stayed under the per-run budget.

## 10. Files and commands

- `scripts/variant_knockouts.py` — run 3 produced `variant_knockouts.json` and `solver_calls.jsonl`, and the full
  sweep TSV in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/or_rule_review_r171_aars_20261029/`.
- `scripts/build_gene_table.py` — `or_rule_review_genes.csv`.
- `scripts/nterm_heuristic.py` — the N-terminal table of section 7.2.
- `evidence/` — UniProt TSV (`303f78b1…`), per-screen rows, README with the queries.

```bash
OR_RULE_SCRATCH=<scratch dir> <worktree venv python> -B model/reports/or_rule_review_r171_aars_20261029/scripts/variant_knockouts.py
python3 model/reports/or_rule_review_r171_aars_20261029/scripts/build_gene_table.py
python3 model/reports/or_rule_review_r171_aars_20261029/scripts/nterm_heuristic.py
```
