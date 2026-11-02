# Ceramide pool in `biomass_C`: report (2026-11-01)

Scope, authorization and acceptance criteria: [TASK.md](TASK.md). Independent audit: [AUDIT.md](AUDIT.md): 7 claims
supported, C3 (the guarded step) supported with minor gaps (§7). Exploratory evidence (FVA, TEST screens, bypass
check, literature, build copies): `$IYALI26_RESEARCH_ROOT/artifacts/tasks/ceramide_pool_biomass_20261101/`.

## 1. What changed and why

**Before.** `biomass_C` had no sphingolipid term. The legacy pseudoreaction R1371 that once drained "complex
sphingolipid" is closed. So nothing in the model needed the sphingoid backbone. An FVA of 97 LPs on `c1f70ac4` on
2026-11-01 found serine palmitoyltransferase R722, 3-ketosphinganine reductase R37 and every SUR2 and SCS7 reaction
at min = max = 0 (verified). Their genes were consensus-essential false negatives.

**Change.** Candidate `2bc263f23603c27131d0b07143300a484dbd5a135ef7151665574400ab64a062`
([file](../../candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_20261101.xml))
is `c1f70ac4` plus `--ceramide-pool-biomass`:

- **Pool metabolite:** a new pseudo-metabolite `sphingolipid_ceramide_pool[C_er]`, with no formula and charge 0.
- **Feed reactions:** 10 gene-less pseudo-reactions `xPOOL_CER_<mid>`, each `<mid>[C_er] → pool` with bounds
  (0, 1000), for m386, m388, m390, m391, m392, m393, m397, m398, m399 and m400. These are ceramide-1, -2, -2′, -3
  and -4, at C24 and C26.
- **Biomass term:** `biomass_C` consumes the pool at **5.4e-4 mmol/gDW**, with 7 `ceramide_pool_biomass_*` notes.
  `biomass_C` is not re-normalised.

**Why a pool instead of one species.** A pool can be filled by any mix of ceramides, so it requires a sphingolipid
backbone without requiring any one modification. A *Y. lipolytica* sur2Δ strain is viable but grows slowly, and its
sphingolipids shift to the dihydro form (Shin 2025 abstract and the same group's patent, verified). A drain on phyto
or α-OH ceramide alone would therefore have forced SUR2 or SCS7 to be essential against that evidence.

## 2. Evidence, by tier

| Claim | Source | Tier |
| --- | --- | --- |
| Total long-chain bases are 0.7 % of polar lipids; polar lipids are 52.3 % of lipid; lipid is 4.9 % of DCW; 19-phytosphingosine is 90.8 % of LCB | [Rupcic 1996](https://europepmc.org/article/MED/8997745), *C. lipolytica* on methanol | verified (abstract, opened by a review agent and by the auditor) |
| 0.007 × 0.523 × 0.049 = 1.794e-4 g/gDW; ÷ 331.54 g/mol (C19H41NO3) = **5.41e-4 mmol/gDW** | arithmetic | inferred (re-done by the auditor) |
| One long-chain base per sphingolipid, so moles of LCB = moles of sphingolipid | established chemistry | inferred / uncited |
| Sphingolipids are "a minority" of W29 membrane lipids; membrane lipids are 12.98 mg/g | [Ianutsevich 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6950712/) | verified (full text) |
| Fig. 4A: SL is about 1.3–2 % of membrane lipids at pH 5.5, 29 °C, so about 0.15–0.25 mg/gDW, or 1.5–3.7e-4 mmol/gDW as intact lipid | same | figure reading (I viewed the figure); no number in the text |
| Ylsur2Δ is viable with retarded growth | [Shin 2025](https://europepmc.org/article/MED/39603335) | verified (abstract); full text closed access, unread |
| lcb1, lcb2, tsc10 and aur1 nulls are inviable | *S. cerevisiae*, SGD | homolog phenotype, checked by a review agent; no *Y. lipolytica* deletion data found |
| Neither yeast-GEM 9.1.1 nor iYali has a sphingolipid biomass term | model files | verified by a review agent |

The term is **provisional**, like the r_4598 cofactor terms. The value comes from methanol-grown cells measured by
TLC/GC in 1996 (abstract only), not from PO1f on the Task A medium, and the Rupcic figure includes free bases. The
two *Y. lipolytica* sources agree at the mass level (about 0.15–0.25 mg/gDW). The molar value spans about 1.5–5.4e-4
depending on whether LCB or intact-lipid mass was weighed.

## 3. Results (verified this run)

| Check | Result |
| --- | --- |
| pytest after vs before | No new failures: the same 19 failed and 20 errors; the 5 new tests pass |
| Control build without the flag | `c1f70ac4…7540`, byte for byte |
| Seeds 0 and 1 | Both `2bc263f2…a062` |
| Whole-model diff vs `c1f70ac4` | Exactly 13 entries: 1 metabolite, 10 reactions, the `biomass_C` stoichiometry (+pool, −5.4e-4) and its notes ([model_diff_vs_c1f70ac4.json](model_diff_vs_c1f70ac4.json)) |
| WT (Task A + PO1f, Gurobi 1e-9) | 0.3122964 → **0.3120535** (−0.078 %), optimal, checks passed ([wt_energy.json](wt_energy.json)) |
| Energy check (no uptake, `xMAINTENANCE` lower bound 0, maximise `xMAINTENANCE`) | 0.0 on both models |
| Certified knockout screen | 1,082 knockouts, 40 certificates, 0 fresh re-solves ([ko_screen_ceramide_pool.tsv](ko_screen_ceramide_pool.tsv), [summary](ko_screen_ceramide_pool_summary.json), [ko_compare.json](ko_compare.json)) |
| LPs used | 1,127 of the 1,300 cap |

**Class changes** (at every cutoff, all from KO/WT 1.0 to growth 0, optimal; no other gene changed):

| Gene | List label (homology only; no *Y. lipolytica* name verified) | Reaction(s) | Experimental | Result |
| --- | --- | --- | --- | --- |
| YALI1C04265g | LCB1-like, serine palmitoyltransferase | R722 (AND) | essential | FN → TP |
| YALI1F20637g | LCB2-like | R722 (AND) | essential | FN → TP |
| YALI1B22883g | TSC10-like, 3-ketosphinganine reductase | R37 | essential | FN → TP |
| YALI1B26350g | ELO/SUR4-like elongase | R1406–R1412 | essential | FN → TP |
| YALI1F15781g | PHS1-like dehydratase | R1420–R1426 | essential | FN → TP |
| YALI1A05116g | TSC13-like enoyl reductase | R1427–R1433 | essential | FN → TP |
| YALI1A06453g | IFA38/KAR ortholog (reciprocal best BLAST, 56 % identity) | R1413–R1419 | not listed | TN → **FP** |

**Scores** (rule B, mtDNA + untargeted genes out of scope; user-defined negatives; development reference, not
independent validation):

| Cutoff | Before (TP/FN/FP/TN, MCC) | After |
| --- | --- | --- |
| 1 % | 108/218/36/702, 0.381 | 114/212/37/701, 0.396 |
| 5 % | 115/211/38/700, 0.396 | 121/205/39/699, 0.411 |
| 10 % | 119/207/43/695, 0.394 | 125/201/44/694, 0.408 |
| **15 %** | **136/190/43/695, 0.442** | **142/184/44/694, 0.456** |

## 4. Mechanism (inferred from the FVA and the screen)

The pool needs one ceramide per unit. The only ER sphingoid-base entry is R722 → R37 (serine + palmitoyl-CoA →
sphinganine). Ceramide also needs C24/C26-CoA, which only the ER elongation cycle makes:
- the ELO condensation;
- the KAR reduction;
- the PHS1 dehydration;
- the TSC13 reduction.

Each of these steps has a single gene, so every one becomes essential.

Other routes keep three genes non-essential:
- Ceramide synthase YALI1F09542g is bypassed by the cytosolic copy R2127 (`YALI1F09542g or YALI1D05038g`, verified
  by 9 in-memory LPs on 2026-11-01). Cytosolic ceramide-1 m1980 reaches the ER pool feed m386 through the gene-less,
  reversible transport R2128. m1980 has no formula, so R2127 and R2128 are both elementally unbalanced as written
  (verified static, 2026-11-02).
- In the test variant's WT optimum the ER dihydroceramide comes from R267, the reversible dihydroceramidase
  (`YALI1D25664g`), carrying positive flux in its written condensation direction (sphinganine + lignoceric acid →
  ceramide-1 (C24) + H2O + H⁺, the reverse of ceramide hydrolysis), not from the ceramide synthase; R268 is its C26
  copy (`bypass_check.json` in the research-root folder). Both are charge-imbalanced by +1 (verified static,
  2026-11-02).
- SUR2 and SCS7 are optional because the dihydro ceramide-1 can fill the pool.

## 5. The new false positive YALI1A06453g

- **Identity:** the IFA38/KAR ortholog. It is the reciprocal best hit to *S. cerevisiae* IFA38 (P38286), and UniProt
  Q6CHP1 calls it a very-long-chain 3-oxoacyl-CoA reductase by homology (verified by a review agent).
- **In *S. cerevisiae*:** the null is viable but very slow-growing, because Ayr1 supplies the leftover activity. The
  double mutant is inviable ([Han 2002](https://europepmc.org/article/MED/12087109), abstract, verified by a review
  agent).
- **In *Y. lipolytica*:**
  - No backup reductase has been identified, and there is no clear AYR1 ortholog.
  - The CRISPR calls split: Cas9 nonessential (FS −1.40, q = 1.0), Cas12a essential (FS −2.16, q = 0.004). I checked
    this in the local fitness file.
  - The gene is not on the consensus list.
- **Conclusion (inferred):** the single-gene rule is supported by orthology. The false positive likely reflects an
  unknown backup activity, or a slow-growth phenotype that the screen calls non-essential. There's no evidence for an
  OR partner, so it stays as a known limitation.

## 6. What is not tested or not settled

- One medium (Task A) and one strain profile (PO1f). The coefficient is not a PO1f measurement.
- The full text of Shin 2025 is unread. No *Y. lipolytica* deletion data exist for LCB1, LCB2, TSC10 or the
  elongation genes; their essentiality rests on the screen plus *S. cerevisiae* orthology.
- The consensus essential labels for SUR2 and SCS7, against Ylsur2Δ, which is viable with retarded growth.
- The R2127 ceramide-synthase OR rule (`YALI1F09542g or YALI1D05038g`); both genes are consensus positives.
- The elongation-cycle chemistry (verified static):
  - R1414 (5 H⁺) and R1427 are charge-imbalanced by −5; R1407 and R1420 by +4.
  - R1421 and R1428 are elementally imbalanced: trans-hexadec-2-enoyl-CoA m1449 has no formula.
  - Owner unknown, possibly the lipid line. They carry the forced flux in this candidate but were not changed.
- R37 (3-ketosphinganine reductase, now forced by the pool) is unbalanced by charge −2 and H −1 as written; the
  stored NADP⁺ formula is C21H29 while its name says C21H28 (verified static, 2026-11-02). This is the same
  protonation-convention question as the elongation cycle; not changed.
- The routes into the pool that need no ER ceramide synthase (R2127 → R2128; R267/R268 in the condensation direction) are all unbalanced
  as written (§4).
- The colleague's `lipid_unlump/missing_lipids.csv` suggests ceramide-2, -3 and -4 formulas one O lower than this
  model has. This is the colleague's line, review only. The pool lumps all 10 species, so it isn't affected now, but
  the step pins the current names and formulas and will fail closed if they change.
- The complex-sphingolipid copies (census UM1) remain disconnected. They went to the disconnected-pathway session.

## 7. Audit follow-up

- **Gaps 1–3** (note values not pinned, spec pins come from data, order enforced only by the composition pin): kept.
  They follow the trehalose step's convention, and the composition pin catches a misplaced step at build time.
- **Gap 4** (no notes on the pool and the 10 pseudo-reactions): not changed in this round, because it would need a
  rebuild outside the approved scope. It is a candidate follow-up.
- **Gap 6:** this report and the STATE.md entry close it.

## 8. Next gates

1. Kept (2026-11-02); to land on main after PR #43.
2. Optional follow-ups, each its own round:
   - the elongation-cycle chemistry;
   - provenance notes on the pseudo-reactions;
   - the R2127 OR rule;
   - a *Y. lipolytica* sphingolipid measurement in PO1f.

## 9. Files and commands

- **Scripts:** [scripts/](scripts/) holds `before_state.py`, `run_builds.zsh` (with `args_c1f70ac4.json`),
  `run_pytest_pair.zsh`, `model_diff.py`, `wt_energy.py`, `ko_screen_certified.py` (copied unchanged from the
  baseline round) and `compare_screens.py`.
- **Logs:** [logs/](logs/) holds the build commands and results, the pytest before/after logs, the screen log and the
  input SHAs.
- **Build** (from `platform/`): `PYTHONHASHSEED=0 /opt/miniconda3/bin/python -m scripts.gem_annotate <c1f70ac4 argv>
  --ceramide-pool-biomass --output-model ../model/candidates/<name>.xml`. The exact text is in
  [logs/build_commands.txt](logs/build_commands.txt).
- **Exploratory work this round built on** (2026-11-01, scratch only): the FVA, the in-memory bypass check, and the
  scratch-SBML test screen. Summaries are in the session's shared notes, `blocked_wt_fn_20261101.md`.
