# R2121 malonyl-CoA bypass: Stage 1 (in memory, no build), 2026-11-10

## 0. Authorization

The user, 2026-11-10, in the reaction-activation-options session:

1. "A2 malonyl-CoA first", choosing A2 from the ranked options for zero-flux reactions on e584151f.
2. "Approve Stage 1", approving the Stage 1 plan as written in that session (reproduced in sections 4–6).

Scope: computation and documentation only. No model XML, curated data, builder code or test changes; no build;
no commit or push. Stage 2 (the curation change itself) needs a separate approval.

## 1. Starting state (checked 2026-11-10)

- Branch `claude/reaction-activation-options-1ca120`, HEAD `154eb64ca58f827863bdf061e5b4d03b4db05f4b`, clean.
- Model executed: `model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_20261103.xml`,
  SHA-256 `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2`, read through a short-path copy with
  the same SHA (cobra rejects paths of 260 characters or more on macOS).
- Reference model: none designated (`model/candidates/README.md`).
- Condition: Task A medium `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) plus the PO1f
  profile `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` (`35307853…`), applied in memory.
- Baseline screen: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/mcnaugh_reactions_review_20261103/screen_v4_20261104/ko_screen_v4.tsv`
  (certified, Gurobi 1e-9; WT 0.3124885610076883; at 15 %, rule B: 141 / 185 / 43 / 695).
- Target state: `before.json` in this folder.

## 2. Question

In WT pFBA on e584151f, all cytosolic malonyl-CoA (0.186 mmol/gDW/h) is made by R2121, and R2121 also supplies
60 % of cytosolic NADPH (0.186 of 0.31). The route is Asp → β-alanine (R1936) → malonate semialdehyde (R2120) →
malonyl-CoA (R2121). The cytosolic acetyl-CoA carboxylase R88 carries 0. Stage 1 asks what the model does when
R2121 is closed, before any curation change is proposed.

## 3. Evidence for the proposed closure

| Claim | Source | Tier |
| --- | --- | --- |
| R2121 is KEGG R00740, EC 1.2.1.75, malonyl-CoA reductase (malonate semialdehyde-forming) of the 3-hydroxypropanoate CO2-fixation cycles of green non-sulfur bacteria and thermoacidophilic archaea | [R00740](https://www.kegg.jp/entry/R00740), [EC 1.2.1.75](https://www.kegg.jp/entry/1.2.1.75) | verified (KEGG REST, 2026-11-10) |
| YALI1C02550g is yli:2909116, K00140 (EC 1.2.1.18 / 1.2.1.27); that ortholog group catalyses R00705, semialdehyde + CoA + NAD⁺ → acetyl-CoA + CO2 + NADH, which is the model's R2119 (same gene, irreversible, WT flux 0) | [yli:2909116](https://www.kegg.jp/entry/yli:2909116), [R00705](https://www.kegg.jp/entry/R00705); model annotation `ec-code 1.2.1.27` | verified |
| KEGG assigns no K15017 (malonyl-CoA reductase) gene to Y. lipolytica (yli) or S. cerevisiae (sce) | [K15017](https://www.kegg.jp/entry/K15017) | verified (annotation-level absence, not an assay) |
| Yeast-GEM carries R00705 (`r_4573`, irreversible) and not R00740; neither Yeast-GEM nor eciYali5-GEM carries R00740 (corrected 2026-11-12 after the audit: "only" was wrong, Yeast-GEM also has the NADP form R00706, `r_4574`) | local copies in `$IYALI26_RESEARCH_ROOT/reference/external_models/` | verified (scan by KEGG id and name) |
| R2121 is inherited unchanged (reversible, YALI1C02550g) from iYli21 and the immutable source | `model/source/iyli21.xml`, `model/source/iyali26.xml` | verified |
| ACC1 is essential in S. cerevisiae | [Hasslacher 1993](https://pubmed.ncbi.nlm.nih.gov/8098706/) | verified (abstract) |
| Y. lipolytica ACC raises malonyl-CoA when expressed in S. cerevisiae | [Pereira 2022](https://pubmed.ncbi.nlm.nih.gov/36284710/) | verified (abstract) |
| Flagged before as a malonyl-CoA compensation route, undecided | `model/reports/model_static_screen_20260911/REPORT.md` | recalled (read 2026-11-10) |

No evidence for R2121 activity in yeast was found. Against the closure: nothing found; absence from KEGG and from the
yeast models is annotation-level evidence.

## 4. Genes

No name below is experimentally verified for Y. lipolytica; names are RefSeq/KEGG annotations. All six pass the W29
tag-protein audit (`ok` in `model/curation/annotation/w29_tag_protein_audit.csv`).

| Gene | Name (database) | Function | Reactions | v4 screen |
| --- | --- | --- | --- | --- |
| YALI1C02550g | methylmalonate-semialdehyde dehydrogenase (acylating), K00140 | (methyl)malonate semialdehyde → acyl-CoA + CO2 | R2112, R2119, R2121 | TN |
| YALI1E21915g | 4-aminobutyrate aminotransferase, K13524 | GABA / β-alanine transamination | R50, R2066, R2120 | TN |
| YALI1F11827g | glutamate decarboxylase, K01580 | Glu → GABA | R328, R1936 | TN |
| YALI1C23902g | "L-aspartate decarboxylase dtxS4", K01580 | Asp → β-alanine | R1936 | TN |
| YALI1C15991g | acetyl-CoA carboxylase, K11262 | acetyl-CoA → malonyl-CoA | R87, R88, R1393 | TP |
| YALI1E35955g | biotin–protein ligase, K01942 | biotinylation of carboxylases | R175, R88 (AND) | FN |

## 5. Stage 1 computation

Gurobi 13, FeasibilityTol / OptimalityTol 1e-9, the checked solve of
`model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py` (`69c42690…`), and the certified sweep of
`model/reports/ceramide_pool_biomass_20261101/scripts/ko_screen_certified.py` (`313d5d01…`). Variant V1 = e584151f
with R2121 bounds set to (0, 0) in memory after the condition is applied.

1. Baseline: reproduce WT growth (0.3124885610076883 within 1e-9).
2. V1: WT growth and pFBA; the producers of cytosolic malonyl-CoA, acetyl-CoA and NADPH.
3. FVA at 99.9 % of the optimum, baseline and V1, on 14 named reactions only: R88, R2121, R2119, R2120, R1936,
   R1575, R2003, R191, R90, R325, R639, R550, R175, R87 (56 LPs).
4. Energy check, baseline and V1: every exchange closed to uptake, the `xMAINTENANCE` lower bound relaxed to 0,
   maximise `R_NTP1`.
5. Certified single-gene sweep on V1 (1,081 knockouts), scored at 1, 5, 10 and 15 % under rules A and B with the
   mtDNA and untargeted scopes, and compared gene by gene with the v4 screen.
6. Controls on V1: knockouts of YALI1E35955g, YALI1C15991g and YALI1C02550g, and the YALI1E35955g knockout with the
   R88 rule set to YALI1C15991g alone (shows whether a YALI1E35955g call comes only from the AND rule).

Budget: at most 1,400 LPs; no builds.

## 6. Stop conditions

Stop and report, without continuing to the sweep or proposing Stage 2, if:
- any solve of steps 1–4 is non-optimal, or the baseline does not reproduce within 1e-9;
- the energy check is above 0;
- V1 WT growth drops more than 5 %;
- another route takes over the cytosolic malonyl-CoA supply (R88 FVA minimum 0 in V1);
- in the sweep, any TP becomes FN, or any knockout grows more than 1e-6 above WT.

## 7. Outputs

- This TASK.md and `before.json` (repo).
- Scripts, solves, sweep TSV, scores and a manifest in
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2121_malonyl_coa_20261110/`.
- Results are reported in the session; REPORT.md and a STATE.md entry follow at the end of the round.

Scoring note: model genes not on the consensus list count as non-essential (user-defined negative class); mtDNA and
untargeted genes are out of scope; the list was used in development, so scores are a development reference.

## 8. Stage 2: close R2121 through the builder (authorized 2026-11-12)

**Authorization.** The user, 2026-11-12, after reading the Stage 1 results and the Stage 2 proposal (REPORT.md §6):
"Approve Stage 2". The proposal approved: close R2121 to (0, 0), keeping the reaction and its rule, with a note
giving the KEGG evidence, through curated data and the builder (skill `curation-change`), with an Escher
before/after figure, against the acceptance criteria below. The two open decisions of REPORT.md §6 were not
answered: the `YALI1E35955g` call stays reported as an effect of the R88 AND rule, and R1936 and the chitin
coefficient stay unaddressed here. No commit or push.

**Baseline.** e584151f (path in section 1, SHA-256 `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2`),
built by the argv in its `.build.json` (HEAD `048e18ec` then; the V-ATPase A1 step merged since is opt-in and not
used). HEAD now `154eb64c`.

**Targets.**
- Changed: R2121 bounds [−1000, 1000] → [0, 0], plus notes. Reaction, stoichiometry, GPR, annotation and existing
  notes kept.
- Pinned, not changed: R2119, the retained reaction of YALI1C02550g (R00705 chemistry).
- Nothing else may differ.

**Implementation.**
- Data: `model/curation/direction/r2121_malonyl_coa_close.json` (evidence fields: source URLs, tiers, reviewer).
- Step: `platform/scripts/gem_annotate/r2121_malonyl_coa_close.py`, `apply_r2121_malonyl_coa_close`; guarded,
  idempotent, fails closed when R2121 or R2119 differs from the reviewed e584151f state.
- Flag: `--r2121-malonyl-coa-close`, opt-in, run as the last model step (after `--bigg-crosswalk-annotations`; it
  writes bounds and notes only).
- Test: `platform/tests/test_r2121_malonyl_coa_close.py` on an in-memory e584151f.

**Expected effect.** The Stage 1 V1 problem: WT 0.3120526953810365; at 15 % (rule B, mtDNA + untargeted out of
scope) 143 / 183 / 43 / 694.

**Acceptance criteria.**
1. A control build without the flag reproduces e584151f byte for byte.
2. The whole-model diff of the flag build against e584151f shows only R2121's bounds and the new notes.
3. Seed 0 and seed 1 builds are byte-identical.
4. WT growth on the built candidate equals 0.3120526953810365 within 1e-9; no-uptake ATP maximum 0 (optimal).
5. The certified sweep on the built candidate reproduces `ko_screen_v1.tsv` class by class at every cutoff, with KO
   growth within 1e-9.
6. pytest: the new tests pass and no existing test changes outcome.
7. Code review and independent audit leave no unresolved blocking finding.

**Budget.** 3 builds (control, candidate, seed-1 rebuild) plus 1 spare; at most 1,400 LPs (sweep about 1,123, WT
and energy checks, pFBA for the figure).

**Stop conditions.** Any acceptance criterion fails, or the diff shows anything outside the targets: stop and
report, without retrying or relaxing the criterion.
