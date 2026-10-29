# R171 ATP synthase rule: OR branches g and k removed, subunits e and j out of the AND (Option D), candidate ca43ae0c (2026-10-29)

Scope and authorization: [TASK.md](TASK.md) (user, 2026-10-29: "Option D (Recommended)" for the plan and the
baseline `d8d1bac7`). Evidence review: [../or_rule_review_r171_aars_20261029/REPORT.md](../or_rule_review_r171_aars_20261029/REPORT.md)
§4, §5, §7.1. Audit: [AUDIT.md](AUDIT.md). Figure: [escher/](escher/README.md). Nothing committed or pushed.

## 1. What changed and why

R171 (ATP synthase, EC 7.1.2.2, 4 H⁺ cytosol → matrix per ATP) read, as inherited from iYli21,
`(14-subunit AND) or YALI1B28025g or YALI1B15733g`. The two OR branches are subunits g (ATP20-like, 138 aa) and k
(ATP19-like, 70 aa). An OR branch says the gene alone catalyses the reaction; there is no evidence that a single
accessory subunit synthesises ATP (AGENTS §5), so every single knockout of the 16 genes left R171 open.

The user chose Option D: R171 is now the 12-gene AND of ATP6/8/9 (mtDNA, `YALI1M00260g`, `M00267g`, `M00318g`) and the
nine nuclear core subunits (α `F04433g`, β `B05364g`, γ `F04110g`, δ `D27826g`, OSCP `D15734g`, b `F26967g`, d `B08885g`,
h `F07384g`, f `D14639g`). g, k, e (`YALI1E38131g`) and j (`YALI0D17490g`) leave the rule and stay in the model with no
reaction. Catalytic ground (from the review, sources opened there and not re-opened here; *recalled* in this round):

| Subunit | Ground for leaving the Boolean rule | Level |
| --- | --- | --- |
| g, k | No evidence either catalyses ATP synthesis alone. Dimer-specific: the *Y. lipolytica* monomer lacking e, g, k hydrolyses ATP at about the dimer's rate (Hahn 2016, [PMID 27373333](https://europepmc.org/article/MED/27373333)); g and k "not found" in the purified *Y. lipolytica* enzyme (Liu 2015, [PMID 25759169](https://europepmc.org/article/MED/25759169)). | recalled (review §7.1, verified there); "not catalytic" inferred |
| e | Dimer-specific (Hahn 2016); *S. cerevisiae* e, g, k "not required for the formation of enzymatically active ATP synthase" (Arnold 1998, [PMID 9857174](https://europepmc.org/article/MED/9857174)). | as above |
| j | *S. cerevisiae* atp18Δ viable with about half the specific ATPase activity and a low ADP/O (Vaillier 1999, [PMID 9867878](https://europepmc.org/article/MED/9867878), abstract): loss impairs, does not abolish, activity. Not screened in *Y. lipolytica*. | recalled; the keep-or-drop call is a judgement on one *S. cerevisiae* abstract |

**Evidence against removal** (review §7.1, recalled): Hahn 2016 measured ATP *hydrolysis*, not synthesis, by the
monomer; g's Cas9 fitness score −2.80 (q 0.085) suggests a fitness cost; Liu 2015 lists e among the *Y. lipolytica*
subunits; UniProt (Q6C8S0) places j in the *Y. lipolytica* enzyme, and the *S. cerevisiae* j null shows j is needed
for full coupling and respiratory growth. Dimer loss in *S. cerevisiae* affects cristae and respiratory fitness
(recalled in the review, not lethality). None of this shows that losing e, g, k or j abolishes the activity.

The AND is read strictly: a gene is a member if its loss is taken to close the reaction. Under that reading A and C
keep subunits whose loss is not known to abolish activity. AGENTS §5 asks only for "evidence of a complex or
co-requirement" for an AND member, and e and j are complex members, so the strict reading is a choice, not the only
admissible rule. **The reading is applied unevenly:** for the retained nuclear members (in particular b, d, f, h and
OSCP) the review has no enzymological evidence that their loss abolishes activity; their membership rests mainly on
the consensus screens, which are the benchmark list. So the +9 TP at 15 % is not independent of the curation
evidence (beyond the general "development reference" caveat in §5). Not reasons for the choice: the recall or MCC
gain, the FP count of each option (A +2, C +1, D 0 in the review), Yeast-GEM's 17-gene AND.

**Gene names.** None of the 16 gene names is a verified *Y. lipolytica* name; "ATP1-like", "subunit g" etc. are
homolog names. The subunit assignment of the nine nuclear core genes comes from the model's own UniProt
cross-references (`YALI1F04433g` Q6C326 ATP1 α, `B05364g` Q6CFT7 ATP2 β, `F04110g` Q6C338 ATP3 γ, `D27826g` Q6C877
ATP16 δ, `D15734g` Q6C9B1 ATP5/OSCP, `F26967g` Q6C105 ATP4 b, `B08885g` Q6CFH9 ATP7 d, `F07384g` Q6C2V6 ATP14 h,
`D14639g` Q6C9E6 ATP17 f), checked against UniProt REST by the independent audit this round; the four removed
genes' identities are from the review (§4, UniProt verified there). These are database identities, not experimental
validation.

## 2. Implementation

| Part | File |
| --- | --- |
| Decision and evidence (data) | [model/curation/gpr/r171_atp_synthase_gpr.json](../../curation/gpr/r171_atp_synthase_gpr.json), written by [scripts/write_curation.py](scripts/write_curation.py) from [before.json](before.json) |
| Guarded step | `platform/scripts/gem_annotate/r171_atp_synthase_gpr.py`, `apply_r171_atp_synthase_gpr` |
| Flag | `--r171-atp-synthase-gpr` (`cli.py`), called in `main.py` after the aaRS split step; no flag dependency |
| Test | `platform/tests/test_r171_atp_synthase_gpr.py` (6 tests, in-memory model, no LPs) |

The step follows `aars_thr_ser_split.py`. Before any edit it checks the spec (id, status, before / after rule, the 16
genes, the four unlinked genes and their reaction lists) and R171's name, bounds, stoichiometry, species (formula,
charge, compartment), annotation, rule and notes, plus each gene's annotation and reaction list. A stale or foreign
state raises `ValueError` with no edit; the curated state returns `already_correct`; notes are restored on cobra
context exit. The test covers 13 stale-input cases, 6 scope cases, a whole-model state comparison and an SBML
roundtrip, the knockout logic before and after (every single knockout leaves R171 open before; after, each of the 12
closes it and none of the 4 does), context exit and the opt-in flag.

One fix during implementation (verified): the first data file had apostrophes in two note strings; libSBML writes
them as `&apos;`, which the roundtrip test caught. The note text was reworded without apostrophes and regenerated.

## 3. Builds (verified)

From `platform/` at HEAD `ff7fc62` plus the uncommitted change ([logs/build_head.txt](logs/build_head.txt)),
`/opt/miniconda3/bin/python`; commands in [logs/build_commands.txt](logs/build_commands.txt).

| Build | Seed | Output SHA-256 |
| --- | --- | --- |
| Control: `d8d1bac7` argv (flag off) | 0 | `d8d1bac78138f59ba6beefcbddae7f93d0a20e245fea6f5dd9c79760ee535751`, byte-identical to the baseline |
| Candidate: + `--r171-atp-synthase-gpr` | 0 | `ca43ae0c12e01c13c7b2ca4792b745d00d4a28ed950f8e401aece37b2a0843d0`, [candidates/E5_…_aars_thr_ser_split_r171_atp_synthase_gpr_20261029.xml](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_20261029.xml) |
| Candidate | 1 | `ca43ae0c…`, identical |

Both candidate logs report "R171 ATP synthase GPR: applied"; the control log has no such line. All three logs carry
the pre-existing libSBML `R_R304` notes error (seen in every recent build; not investigated). The control's and the
seed-1 build's `.build.json` are in `logs/`; their XMLs were scratch files.

## 4. Whole-model diff (verified)

`model_diff.py` (unchanged copy, `8a464ff8…`) vs `d8d1bac7` ([model_diff_vs_d8d1bac7.json](model_diff_vs_d8d1bac7.json)):
6 differences, all intended: R171 rule, R171 notes, and the reaction lists of `YALI0D17490g`, `YALI1B15733g`,
`YALI1B28025g`, `YALI1E38131g` (`R171` → none). Gene list unchanged (1,081 genes in both, sweep summary). No
stoichiometry, bound, formula or charge changed, so mass and charge balance of R171 are as before.

## 5. Knockouts (verified)

Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) + PO1f profile `po1f_sd_leu.json` (`35307853…`),
list `consensus_essential_genes.csv` (`1e887f5a…`; [logs/sweep_inputs_sha256.txt](logs/sweep_inputs_sha256.txt));
Gurobi 13.0.2 at 1e-9, reset per knockout, targeted check (`ko_compare_two.py`, unchanged copy, `bb3118a6…`, with
`checked_solve.py` from the 2026-10-22 task); models read from SHA-checked short copies
([logs/sweep_short_copies_sha256.txt](logs/sweep_short_copies_sha256.txt)); 2,164 LPs;
[ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv).

- A first attempt stopped before any LP: cobra rejected the 260+-character absolute path
  ([logs/ko_compare_attempt1_longpath.log](logs/ko_compare_attempt1_longpath.log)).
- WT 0.3116608216 on both (residual 2e-14).
- Growth changes above 1e-6: exactly the 12 AND members, 0.31166 → 0.034190 (KO/WT 0.10970); no other gene.
- Class changes: at 15 % only. The 9 nuclear core subunits FN → TP; the 3 mtDNA genes also flip but are out of scope
  (benchmark §Labels). No change at 1, 5 or 10 %. g, k, e, j: KO/WT 1.0 before and after (TN).
- No failing point, no fresh re-solve, the same 20 infeasible knockouts on both (16 nuclear unresolved + 4 mtDNA);
  max knocked flux 0.

| Cutoff | d8d1bac7 TP / FN / FP / TN | ca43ae0c TP / FN / FP / TN | Recall | Specificity | Precision | MCC |
| --- | --- | --- | --- | --- | --- | --- |
| 1 % | 96 / 219 / 35 / 699 | 96 / 219 / 35 / 699 | 0.305 → 0.305 | 0.952 → 0.952 | 0.733 → 0.733 | 0.356 → 0.356 |
| 5 % | 102 / 213 / 37 / 697 | 102 / 213 / 37 / 697 | 0.324 → 0.324 | 0.950 → 0.950 | 0.734 → 0.734 | 0.370 → 0.370 |
| 10 % | 106 / 209 / 42 / 692 | 106 / 209 / 42 / 692 | 0.337 → 0.337 | 0.943 → 0.943 | 0.716 → 0.716 | 0.368 → 0.368 |
| **15 %** | 114 / 201 / 42 / 692 | **123 / 192 / 42 / 692** | 0.362 → **0.390** | 0.943 → 0.943 | 0.731 → **0.745** | 0.393 → **0.420** |

16 unresolved at each cutoff. Positive coverage: 326 of 1,612 listed genes are model genes, 315 scored, 11 among the
unresolved; 1,286 listed positives are out of scope. Negatives are the user-defined class (model genes not on the
list); mtDNA genes are out of scope; the list was used during development, so these scores are a development
reference, not independent validation.

The result equals the review's in-memory variant `R171_D` (on `965ab947`) added to this baseline: +9 TP at 15 % and
nothing else.

## 6. Mechanism (verified, figure run)

With the alpha subunit knocked out ([escher/](escher/README.md), pFBA, 12 LPs):
- Before: R171 stays open (bounds 0–1000) and every mapped flux equals WT; growth 0.3117.
- After: R171 is (0, 0); growth 0.0342. The chain still oxidises NADH and succinate, faster than in WT (complex III
  R305 9.0 → 18.3, complex IV R304 4.5 → 9.2), and the proton leak R1384 (0 → 155.7) dissipates the gradient that no
  longer drives ATP synthase. ATP comes from substrate-level steps (R642 3.34 → 4.67, R694 3.08 → 4.65, R741 1.16 →
  4.55). This is the review's 0.1097 energy ceiling; that the capped glucose (2.414) cannot pay NGAM and growth
  without oxidative phosphorylation is recalled from review §5, not computed in this round.
- So the 15 % cutoff catches the 9 subunits and the 10 % cutoff does not. *Y. lipolytica* is regarded as unable to
  grow by fermentation (recalled, uncited), so the model still under-predicts the severity; that ceiling is a
  separate line (2026-10-21 entry), not part of this change.

## 7. pytest

Full suite, the workflow-automation venv (Python 3.13.5, cobra 0.30.0, pytest 9.0.2; this worktree has none; it
imports this worktree's `scripts` and `tools`, checked), before vs after
([logs/pytest_before.txt](logs/pytest_before.txt), [logs/pytest_after.txt](logs/pytest_after.txt),
[logs/pytest_compare.txt](logs/pytest_compare.txt)): 324 → 330 passed, the same 19 failed and 20 errors, subtests
858 → 909. Test by test, the only differences are the 6 new tests, all passing; no existing test changed outcome. The
before run started before any code was written and did not collect the new test (checked). The 19 failures and 20
errors are pre-existing and equal to the counts recorded by the aaRS split round.

## 8. Acceptance (TASK.md §4)

| # | Criterion | Result |
| --- | --- | --- |
| 1 | Control reproduces `d8d1bac7` byte for byte | met |
| 2 | Seeds 0 and 1 identical | met |
| 3 | Diff = R171 rule and notes + four gene links | met (6 differences) |
| 4 | Sweep = baseline + the §3 transitions; WT, unresolved unchanged; no failing point | met |
| 5 | Only the new tests differ in pytest | met (§7) |

## 9. Limits, what was not tested, next gates

- The catalytic evidence was read in the review, not re-opened here. No *Y. lipolytica* knockout of e, g, k or j
  exists in the sources; the j decision rests on one *S. cerevisiae* abstract.
- **Model note overstates slightly (audit finding, not fixed in this build):** `r171_atp_synthase_gpr_evidence` in
  R171's notes says the *Y. lipolytica* monomer lacking e, g, k "is catalytically active"; the measured activity was
  ATP hydrolysis, so "hydrolytically active" is accurate (the `…_evidence_tier` note does say "hydrolysis rather than
  synthesis"). Fixing it is a data change and a rebuild to a new file; left for the next revision or the user's call.
  **Fixed later the same day** in the curated data, with the user's approval:
  [../r171_r673_note_revision_20261029/](../r171_r673_note_revision_20261029/REPORT.md), candidate `25dabf62`. That
  supersedes `ca43ae0c` for note wording; this round's builds and results are unchanged.
- The retained AND members rest mainly on the screens (§1), so the gain is not independent of the curation data.
- **Could the method be improved?** The Boolean AND cannot express "loss halves activity" (j) or "loss costs
  fitness" (g); both are forced to all-or-nothing. A test would need a graded model (e.g. a capacity bound scaled per
  knockout) and *Y. lipolytica* knockout phenotypes of e, g, k, j on a non-fermentable carbon source; neither exists
  here.
- `YALI0D17490g` is a CLIB122 identifier inside the YALI1 model (S2 map: `YALI1_D21476g`); it is now unlinked. The
  identifier question is a separate annotation line.
- One medium (Task A + PO1f). The gain exists only at 15 % because of the fermentative energy ceiling.
- A pFBA observation outside this change: WT carries about 18 through R659 / r0210 / R_CYST_II in both models; not
  traced (escher/README.md).
- Provable from the knockout semantics: only the 16 R171 genes can change class under this edit. Statically
  verified: the diff. Observed in runs: the sweep and the figure. Untested: other media, the HTML viewer.
- Next gate: round 2, the prolyl pair R673, needs the user's approval of a TargetP 2.0 / DeepLoc 2.1 submission.

## 10. Independent audit

A fresh read-only subagent re-derived the claims with its own scripts and 14 LPs ([AUDIT.md](AUDIT.md)). Supported:
SHAs, control and seed builds, the whole-model diff (its cobra diff counts R171's gene set as a 7th field; raw XML
diff 2 hunks, both in R171), the rule and unlinked genes, the guarded step, the test (6 passed, 51 subtests), the
sweep spot check and the recomputed matrices, the figure numbers. Partially supported: the catalytic reasoning (the
report omitted the evidence against removal and did not say the retained members rest on the screens; the model note
says "catalytically active" for a hydrolysis measurement) and completeness (pytest, AUDIT.md, STATE, verified-name
statement, method note were missing at audit time). Fixes: §1 now gives the evidence against, the uneven application
of the strict reading and the gene-name statement; §6 marks the glucose-cap sentence as recalled; §9 records the
note wording, the screen dependency and a method note; §7 has the full pytest comparison. The note wording itself is
not changed (would need a rebuild). Cosmetic, not changed: the copied `ko_compare_two.py` docstring describes an
older run.
