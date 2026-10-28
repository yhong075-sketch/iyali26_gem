# TASK: R2135 acyl-CoA pool remap and R2145 1:1 rescale (2026-10-27)

## 0. Authorization

The user, 2026-10-27, after reading the Stage 1 draft of this plan (shared notes folder,
`r2135_acyl_coa_pool_TASK_draft_20261027.md`; its section 12 gives the full approval sentence):

> Approve Option A

"Option A" is section 6 of the draft:
- R2135 becomes `m567[C_em] <=> m1984[C_er]`, with acetyl-CoA `m68[C_cy]` removed.
- R2145 becomes PA + H2O → DAG + Pi, with all coefficients 1.

It is applied by a guarded `--r2135-acyl-coa-pool` step after `--cofactor-biomass-terms` on 00e8e397, with
the test, builds, diff, Gurobi knockout sweep, re-census, review, audit, REPORT and STATE entry of sections
7–10. No commit or push. Sections 1–12 below are the draft unchanged apart from this section, the status line
and the paths of this folder.

The Kennedy R1726 task stays paused (`../kennedy_closures_20261027/STOP.md`). Its code stays uncommitted in
the same worktree and is not built into this candidate.

**Status: authorized 2026-10-27; Stage 2 in progress.** Stage 1 outputs:
- [r2135_fix_variants.json](r2135_fix_variants.json) and the scripts in [scripts/](scripts/);
- `../kennedy_closures_20261027/acyl_leak_diagnosis.json` and `../kennedy_closures_20261027/r2135_control.json`.

## 1. Starting state (read 2026-10-27)

- **Branch / HEAD:** `claude/kennedy-closures-864c91` at `31e7aa4`. It has the Kennedy task's uncommitted
  files, which are preserved and unrelated to this change.
- **`main`:** `e26ac81`; PR #24 is merged, checked with `gh` today.
- **Baseline:** `00e8e39775b64cfeb79cfafff65ef82eaf4760eff1b64fcddc3ce93b40919592`, the cofactor
  biomass-terms candidate, which is on `main`. Not a designated reference.
- **Also tested:** the unplaced R1726 build `95316a85…` (Kennedy task, artifacts folder).
- **Condition and solver:** Task A medium + PO1f; Gurobi 1e-9 via `checked_solve.py`, every reported point
  checked.

## 2. The block and its provenance

- **Inherited unchanged.** R2130–R2200 are in `model/source/iyali26.xml` and `model/source/iyli21.xml`
  with the same equations, bounds and GPRs (verified, cobra comparison). iyli21 is the published model of
  Guo et al. 2022, [doi:10.1016/j.csbj.2022.05.018](https://doi.org/10.1016/j.csbj.2022.05.018). The local
  file was recorded as byte-identical to the publisher's supplementary model by the 2026-09-09 R794 / R795
  review (`model/reports/r794_r795_gpr_review_20260909/provenance.md`; recalled, not re-checked).
- **Rationale:** none recorded; R2135 has no notes.
- **Species:** they carry BiGG suffixes `_TB` / `_LM` (`12dgr_TB`, `tag_TB`, `pa_LM`) and have no formulas.

| Reaction | Equation | Bounds | GPR |
| --- | --- | --- | --- |
| R2135 "AC-CoA transport" | acetyl-CoA `m68[C_cy]` ⇌ Acyl-coA `m1984[C_er]` | ±1000 | none |
| R2134 acyl-CoA:sphingosine N-acyltransferase (KEGG R01496) | `m1984` + ceramide-1 `m386[C_er]` → dihydroceramide + CoA | 0 / 1000 | `YALI1F09542g or YALI1D05038g` |
| R2168 acyl-CoA:DAG O-acyltransferase (KEGG R02251) | `m1984` + DAG `m2011[C_er]` → TAG `m2014[C_er]` + CoA | 0 / 1000 | `YALI1D10264g` |
| R2144 PA transport | PA `m378[C_em]` ⇌ PA `m1993[C_cy]` | ±1000 | none |
| R2145 "diacylglycerol pyrophosphate phosphatase, yeast-specific" | **0.01** PA `m1993` + H2O → **0.01** DAG `m1994[C_cy]` + Pi | 0 / 1000 | `YALI1D35593g` |
| R2169 DAG transport | DAG `m2011[C_er]` ⇌ DAG `m1994[C_cy]` | ±1000 | none |
| R2170 TAG transport | TAG `m1399[C_cy]` ⇌ TAG `m2014[C_er]` | ±1000 | none |

## 3. Identity evidence

| # | Claim | Source | Tier |
| --- | --- | --- | --- |
| I1 | `m1984[C_er]` is annotated MetaNetX MNXM1136895, named "acyl-CoA" (generic). | [MetaNetX MNXM1136895](https://www.metanetx.org/chem_info/MNXM1136895) | verified (page opened) |
| I2 | The ER-membrane pooled acyl-CoA `m567[C_em]` (produced by `xPOOL_AC_EM` from the seven chain-specific acyl-CoAs) carries the same MNXM1136895 and Reactome R-ALL-2534017. So do the lipid-particle `m570` and mitochondrial-membrane `m790`. | model annotations (cobra read) | verified (model) |
| I3 | `m68[C_cy]` is acetyl-CoA, C23H38N7O17P3S, BiGG `accoa`, = KEGG C00024. | model; [KEGG C00024](https://www.kegg.jp/entry/C00024) | verified |
| I4 | The reference reactions of R2134 and R2168 use generic acyl-CoA C00040 (C22H35N7O17P3SR), not acetyl-CoA. | [KEGG R01496](https://www.kegg.jp/entry/R01496), [R02251](https://www.kegg.jp/entry/R02251), [C00040](https://www.kegg.jp/entry/C00040) | verified |
| I5 | R2145 is ModelSEED rxn09551, which itself has the 0.01 coefficients on its "yeast-specific" lumped PA / DAG. | [ModelSEED rxn09551](https://modelseed.org/solr/reactions/select?wt=json&q=id:rxn09551) (Solr query) | verified |
| I6 | `YALI1D35593g` = YALI0D27016g. KEGG K15728 phosphatidate phosphatase LPIN (EC 3.1.3.4), RefSeq "Nuclear elongation and deformation protein 1"; named PAH1 in a Y. lipolytica study. | [KEGG yli:2911239](https://www.kegg.jp/entry/yli:2911239); Szczepańska 2022 [PMC9505966](https://pmc.ncbi.nlm.nih.gov/articles/PMC9505966/) | verified; orthology plus one paper's naming, not an enzyme assay |
| I7 | `YALI1D10264g` = YALI0D07986g. KEGG K11155 diacylglycerol O-acyltransferase 1 (EC 2.3.1.20), RefSeq "Diacylglycerol O-acyltransferase 1". | [KEGG yli:2910950](https://www.kegg.jp/entry/yli:2910950) | verified (orthology). That Y. lipolytica papers call it DGA2 is recalled, not opened. |

The KEGG entries were read through `rest.kegg.jp/get/…`. The `kegg.jp/entry` links are the same entries'
web pages, built from the id pattern.

**Inferred.**
- "AC-CoA" in R2135's name most likely meant acyl-CoA. Its only product feeds two acyl-CoA reactions, and its
  partner `m1984` is a generic acyl-CoA. Mapping it to acetyl-CoA is a species-identity error (cause (b) of
  CLAUDE.md §6a: the local entry points to a different entity).
- Acetyl-TAGs do exist in some plants (established / uncited). Nothing here suggests them in Y. lipolytica,
  and the TAG species `m2014` is the generic TAG.
- R2145's 0.01 coefficients are consistent within ModelSEED's lumped "yeast-specific" scheme. They are not
  consistent with this model's 1:1 links R2144 (PA) and R2169 (DAG) to per-molecule species. Through those
  links, every unit of PA dephosphorylated releases 100× its phosphate.

## 4. What the wild type uses R2135 for (verified)

On the baseline at the optimum (pFBA) and at 99.9 % (FVA):
- **TAG route:** TAG for biomass is made through PA `m378` → R2144 → R2145 → DAG `m1994` → R2169 → **R2168
  with acetyl-CoA from R2135** → TAG `m2014` → R2170 → `m1399` → lipid particle.
- **Forced flux:** R2135 and R2168 carry [1.9e-5, 2.4e-3], R2145 [1.9e-3, 0.24]. All three minima are above 0.
- **Why it is chosen:** it is cheaper than acylating with long-chain acyl-CoA from FAS.
- **Phosphate created by R2145:** at the optimum R2145 runs at 0.0593. It turns 5.9e-4 PA into DAG but
  releases 0.0593 Pi from water. Phosphate uptake is 0.0259 as built and 0.0846 when R2145 is rescaled, so
  about 70 % of the cell's phosphate comes from nowhere. Phosphate is not growth-limiting on Task A, so
  growth does not change.
- **Ceramide synthase R2134 carries no flux** (FVA [0, 0]); its ceramide-1 substrate route is blocked
  elsewhere.

## 5. In-memory fix variants (verified; `r2135_fix_variants.json`; about 400 LPs)

Each variant was applied to the baseline and to the R1726 candidate:
- **remap:** R2135 becomes `m567[C_em]` ⇌ `m1984[C_er]`.
- **close:** R2135 at (0, 0).
- **rescale:** R2145 becomes PA + H2O → DAG + Pi.

| Variant | WT baseline | WT R1726 candidate | Pi uptake (baseline) | FAS / OLE1 / FAA KO in the R1726 candidate | No-uptake ATP / biomass |
| --- | ---: | ---: | ---: | --- | --- |
| as built | 0.3120049 | 0.3120558 | 0.0259 | **1.0 (leak)** | 0 / 0 (optimal) |
| remap | 0.3116824 | 0.3116824 | 0.0846 | 0 (lethal) | 0 / 0 |
| close | 0.3116824 | 0.3116824 | 0.0846 | 0 | 0 / 0 |
| rescale | 0.3120049 | 0.3120558 | 0.0846 | **1.0 (leak)** | 0 / 0 |
| remap + rescale | 0.3116824 | 0.3116824 | 0.0846 | 0 | 0 / 0 |
| close + rescale | 0.3116824 | 0.3116824 | 0.0846 | 0 | 0 / 0 |

**Knockouts (18 targeted genes).** These were the four Kennedy class-change genes; DGAT `YALI1D10264g` and
`YALI1E38810g`; PDAT `YALI1E20049g`; PAH1 `YALI1D35593g`; the ceramide synthases; IPC synthase; the CPT and
PLD genes; PS / CDS / PEM / OPI3 / PIS. On the baseline, no gene changes class in any variant. The four
FAS-group genes, PS / CDS / PEM / OPI3 / PIS stay at 0, and the rest at ≥ 0.999.

**Readings:**
- The remap removes the leak and costs 0.10 % of WT growth.
- Remap and close give identical optima. With the remap, R2168 / R2135 remain usable (FVA up to 6e-3) but are
  not chosen, because acylating with real acyl-CoA costs FAS output.
- The rescale removes the phosphate source but not the acyl leak; it does not change growth.
- **With R2135 fixed, R1726 carries no flux at the optimum.** Its +5e-5 growth gain in the Kennedy task came
  entirely from feeding this leak.

## 6. Options

| | Change | For | Against |
| --- | --- | --- | --- |
| **A (recommended)** | R2135: `m68[C_cy]` → `m567[C_em]`; R2145: coefficients 0.01 → 1 on PA and DAG | Fixes the species error, not just the symptom; keeps the ER DGAT (R2168, K11155 gene) and ceramide-synthase (R2134) routes connected to real acyl-CoA; removes the phosphate source; no class change in the targeted set | Two edits in one inherited, formula-less block; the remap target (ER-membrane pool) is a judgement, not a measured transport; −0.10 % WT growth |
| B | Remap only | Smallest fix of the acyl leak | Leaves R2145's phosphate source, which a later change could again make growth-relevant |
| C | Close R2135 (± rescale) | Simplest; same optimum as A | Blocks R2168 / R2134, which have genes, with no evidence that they do not exist (AGENTS §3: deletion needs non-existence evidence) |
| D | No change | — | The leak stays in the baseline and blocks the Kennedy R1726 change |

## 7. Proposed implementation (Option A)

- **Data:** `model/curation/chemistry/r2135_acyl_coa_pool.json` (unique name).
  - The authorization; the reviewed model `00e8e397…`.
  - Per reaction: the exact before and after stoichiometry; bounds (unchanged); GPR (unchanged); species ids
    with formula, charge and compartment, and the MetaNetX ids of I1–I2.
  - The evidence rows I1–I7 with tier and URL; notes for R2135 and R2145.
- **Step:** `platform/scripts/gem_annotate/r2135_acyl_coa_pool.py`, `apply_r2135_acyl_coa_pool(model)`.
  - It fails closed unless both reactions match the before-state exactly, or both already match the
    after-state with the notes.
  - It touches only R2135's and R2145's stoichiometry and notes.
  - It returns `applied` / `already_correct`.
- **Flag:** `--r2135-acyl-coa-pool`. It requires `--cofactor-biomass-terms`, runs after it, and runs before
  `--r1726-cpt-open`, so the two stack.
- **Test:** `platform/tests/test_r2135_acyl_coa_pool.py`, on an in-memory model, with each guard checked by
  its message.
- **Balance:** not checkable. `m567`, `m1984`, `m1993` and `m1994` have no formulas. The rescaled R2145 is
  written in the per-molecule form of KEGG / EC 3.1.3.4 (PA + H2O → DAG + Pi). This is a stated limit, not a
  pass.
- **Build:** the same `BASE` argv as the Kennedy task (a1237821 flags + `--cofactor-biomass-terms`) plus
  `--r2135-acyl-coa-pool`, output
  `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_<date>.xml`.

## 8. Expected effect (from section 5 unless marked)

- **Whole-model diff:** R2135 stoichiometry (m68 out, m567 in) and notes; R2145 stoichiometry and notes;
  derived reaction lists of m68 and m567.
- **WT:** 0.3120049 → 0.3116824 (−0.10 %). Phosphate uptake 0.0259 → 0.0846 at the pFBA point.
- **Knockouts:** no class change in the targeted set. The full sweep may show small ratio shifts in TAG / PA
  genes (inferred); any class change is a stop-and-explain item.
- **Census (inferred):** R2135 / R2168 / R2145 stay unblocked in the open condition. No newly blocked
  reaction is expected; the re-census checks this.
- **Kennedy follow-up:** with this candidate as base, re-running the Kennedy Stage 2 (`--r1726-cpt-open` on
  top) should restore FAS / OLE1 / FAA essentiality, verified in memory. R1726 then adds no growth on
  Task A.

## 9. Acceptance criteria

1. Guarded step, test and flag as in section 7. pytest is equal to the previous run test by test, plus the new
   file. Code-reviewer findings are fixed or reported.
2. A control build without the flag reproduces `00e8e397…` byte for byte. Seeds 0 and 1 are identical.
3. The whole-model diff is limited to the section 8 entries.
4. WT is optimal and passes the check. No-uptake max ATP and biomass are 0 (optimal, maintenance relaxed).
5. **Two-model knockout sweep** (Gurobi 1e-9, reset per knockout, knocked-flux check), every class change
   explained:
   - repeated with `--r1726-cpt-open` stacked in memory, to show that the four Kennedy genes stay essential;
   - scores at 1 / 5 / 10 / 15 %, with the user-defined negative class, mtDNA out of scope, and a
     development-reference label.
6. Re-census (GLPK, open condition and Task A): no newly blocked reaction unless explained and checked at
   Gurobi 1e-9.
7. Independent audit; REPORT.md; STATE.md entry. The STATE entry also records the paused Kennedy task and its
   stop note.

## 10. Budget and stop conditions

- **Budget:**
  - about 20 LPs re-checked on the built file;
  - 3 builds;
  - 2 × 1,082 knockout LPs, plus 1,082 for the stacked R1726 check;
  - 2 census FVA scans; the census of `00e8e397` is reused from the cofactor task;
  - the pytest suite twice;
  - one review and one audit.
- **Stop and report on:**
  - a guard rejection;
  - a diff outside section 8;
  - seed-dependent bytes;
  - a non-optimal WT;
  - non-zero no-uptake ATP or biomass;
  - a test regression;
  - an unexplained class change;
  - a Gurobi-confirmed newly blocked reaction;
  - the FAS-group genes still non-essential in the stacked R1726 check.

  Nothing is tuned to pass.

## 11. Not in scope (findings only)

- The other 0.01-scaled reactions of the block:
  - IPC synthases R2143 / R2148 / R2151 / R2154 / R2157 / R2161 (0.01 PI + 1 ceramide → 0.01 DAG + 0.01 IPC);
  - MIPC synthases R2146 …;
  - R2171 (PI + 100 H2O → 100 H+ + DAG + 100 inositol-P);
  - R2187 (0.01 PDME + SAM → 0.01 PC; GPR is a YALI0 id).

  The IPC synthases have FVA [0, 0] at 99.9 % on Task A, and R2171 has a maximum of 2e-6. The MIPC synthases
  and R2187 were not checked.
- Two TAG representations linked 1:1 by R2170 (`m1399` / `m2014`), and two cytosolic DAGs (`m1670` /
  `m1994`) linked only through the lipid-particle route.
- The ceramide synthase R2134 is blocked by its ceramide-1 substrate.
- The DAG kinase closure R1720.
- The Kennedy R1726 change: paused, with its stop note; to be re-run after this task.
- Formula-less lumped lipid species: the lipid-unlump line (colleague-owned; review only).

## 12. Authorization sentence (277 characters)

> Approved: R2135 Option A as drafted: remap R2135 to m567[C_em] <=> m1984[C_er] and rescale R2145 to 1:1 via a guarded --r2135-acyl-coa-pool step after --cofactor-biomass-terms on 00e8e397; test, builds, diff, Gurobi KO sweep, re-census, review, audit, REPORT, STATE. No commit.

For Option B, C or D, say so and the scope is redrafted.
