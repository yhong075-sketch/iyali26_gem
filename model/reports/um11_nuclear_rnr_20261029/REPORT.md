# UM11: the nuclear ribonucleotide-reductase copies are kept blocked (Option B, 2026-10-29)

Scope, authorization and the Stage 1 analysis are in [TASK.md](TASK.md). The audit is in [AUDIT.md](AUDIT.md).

**Figure:** [escher/](escher/README.md) (PNG, SVG and an offline viewer). It shows the model as built next to the
Option A what-if on Task A + PO1f, with the two thioredoxin roots marked. The main panel is
[um11_nuclear_rnr_asbuilt_vs_whatif.png](escher/um11_nuclear_rnr_asbuilt_vs_whatif.png).

Evidence tags:
- **verified**: checked this session, with the output kept;
- **inferred**: reasoned;
- **recalled**: not opened this session.

## 1. Decision and what changed

- **Decision.** The user chose Option B ("Option B", 2026-10-29, in this session). The 11 UM11 reactions stay
  as they are:
  - the nuclear RNR copies R707, R709, R711 and R1914;
  - the nuclear UMP kinase R782;
  - the six gene-less nucleotide diffusions R1033, R1076, R1077, R1093, R1136 and R1353.
- **No model, curated-data or code change.** No model file was written.
- **What this task produced:**
  - this folder (TASK.md, Stage 1 outputs and scripts, this report, the audit, the figure);
  - the UM11 paragraph of
    [gap_metabolites_ponce2013_20261026/REPORT.md](../gap_metabolites_ponce2013_20261026/REPORT.md) section 4,
    which now records the decision;
  - a `model/STATE.md` entry.
- **The baseline was built on the user's instruction earlier in the session.** It is
  `14f73bda4b0772019e82ea18a27e198c5ba175e67fef205a51dbfb090809e757`, which is `9d655c2b` plus the UM7 and
  R2171 flags. It is documented separately in
  [combined_um7_r2171_20261029/REPORT.md](../combined_um7_r2171_20261029/REPORT.md).
- Committed in `423e194` on branch `claude/wizardly-wing-a0d6c0` on the user's instruction ("Commit the working tree
  changes", 2026-10-29). Pushed with a pull request to `main` on the user's instruction relayed by the census session
  ("Yes, tell the UM11 session to commit its work on its branch, push it and open a PR to main", 2026-10-29).

## 2. Inputs (verified)

- **Model:** `14f73bda…`, read from SHA-checked short-path copies.
- **Condition:**
  - Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) with the PO1f profile `po1f_sd_leu.json`
    (`35307853…`). Uracil R1354 is open to uptake 1000, and R612 is (0, 0).
  - "Open" means the raw model with every exchange at (−1000, 1000).
- **Solver:** Gurobi 13.0.2 at 1e-9 with the real HOME, through `checked_solve.py` (reset per solve, stale-point
  check).
- **Experimental list:** `consensus_essential_genes.csv` (`1e887f5a…`), 1,612 genes:
  - 326 have a model counterpart;
  - 1,286 are out of scope (no model gene);
  - none is an mtDNA gene.

## 3. The logic: why UM11 is blocked and what would unblock it

Read from the model and checked by LP (**verified**, `stage1/inventory_5f54eb22.txt`, `stage1/stage1_um11.log`).

**The two roots.** Both are nuclear thioredoxin:
- Reduced TRX1 `m931[C_nu]` is consumed by the four RNR copies and has no producer.
- TRX1-disulphide `m933[C_nu]` is produced by them and has no consumer.
- The nucleus has no thioredoxin reductase, no NADPH/NADP species and no thioredoxin transport. The model has no
  glutaredoxin species anywhere.

**Why the rest follows:**
- Every RNR copy is blocked.
- The CDP, GDP and UMP diffusions and the nuclear UMP kinase R782 are blocked, because their nuclear species have
  no other partner.
- The dADP, dCDP and dGDP diffusions are blocked, because their only nuclear producer is a blocked RNR copy.

**Elementary fixes, tested in memory on `14f73bda`:**
- A reduced-TRX link alone frees 0 of 11, and a disulphide link alone frees 0 of 11.
- Both links together (Option A) free 11 of 11.
- A nuclear TRR copy with NADPH/NADP links (Option A′) also frees 11 of 11.
- This matches the census's "chain disconnected at several roots".

**Gene-level bypass:**
- Each nuclear copy carries the same GPR string as its cytosolic counterpart: R707/R709/R711 ↔ R706/R708/R710,
  R1914 ↔ R1913, R782 ↔ R781.
- So no opening can let a nuclear copy bypass a cytosolic knockout.
- The full knockout sweeps confirm this (section 4).

## 4. Results

| What | Result | Level |
| --- | --- | --- |
| UM11 blocked on `14f73bda` (reaction FVA, open and Task A) | 11/11 under both | verified |
| Whole-model blocked, as built | 695 open, 828 Task A | verified |
| Option A, whole-model | 684 open, 817 Task A; exactly the 11 UM11 reactions freed, none newly blocked | verified |
| Option A′, whole-model | 817 Task A, 684 open; the same 11 freed, none newly blocked (scanned by the audit, [AUDIT.md](AUDIT.md) row 5d; not by Stage 1) | verified (audit) |
| WT, Task A + PO1f | 0.3116608216 as built; A, A+ (A plus the R1914 water fix), A′ and C1 within 3e-15 | verified |
| Knockout sweeps of A, A+ and A′ vs the baseline (2 × 1,081 each, all four cutoffs) | no class change, no knockout change > 1e-9 (max 1.6e-13), no stale point, no fresh re-solve; 20 infeasible knockouts, the same in every run | verified |
| Energy-generating cycles (no uptake, `xMAINTENANCE` and every positive lower bound relaxed) | R_NTP1 and biomass maxima 0, status optimal, in every variant | verified |
| pFBA | the nuclear copies carry 0 at the optimum, also when they are open | verified |
| Figure (Task A + PO1f) | as built vs the what-if A, 9,476 LPs; WT reproduced; viewer drawn in jsdom and headless Chrome | verified |

**Benchmark on the baseline** (`14f73bda`, Task A + PO1f, Gurobi 1e-9). Option B leaves the model unchanged, so
these are also the scores after the decision.
- Negatives are the user-defined class (model genes not on the list).
- 16 knockouts are unresolved (infeasible); 11 of them are positives.
- All 16 mtDNA genes are excluded as out of scope; 4 of them are among the 20 infeasible knockouts.
- 1,286 out-of-scope positives have no model gene.
- The list was used in development, so this is a development reference, not independent validation.

| Cutoff | TP | FN | FP | TN | Recall | Specificity | Precision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 % | 94 | 221 | 35 | 699 | 0.298 | 0.952 | 0.729 |
| 5 % | 100 | 215 | 37 | 697 | 0.317 | 0.950 | 0.730 |
| 10 % | 104 | 211 | 42 | 692 | 0.330 | 0.943 | 0.712 |
| **15 %** | **112** | **203** | **42** | **692** | **0.356** | **0.943** | **0.727** |

The 15 % counts equal those recorded for the two parent candidates (STATE.md, 2026-10-28).

## 5. Biology, kept apart from the mechanism

Details and links are in [TASK.md](TASK.md) section 4 and [stage1/evidence.md](stage1/evidence.md). The evidence
was gathered this session by a subagent from PubMed, UniProt and SGD fetches, which are saved in the research root.
From those saved fetches, the main session read the Yao 2003, Lee 2006, Camier 2007 and Pontarin 2008 abstracts
and part of Leeds 1985. The audit opened Yao, Lee, Camier, SGD TRR1 and about ten more sources itself.

- ***Y. lipolytica*:** no data on RNR localization or biochemistry. The gene identities rest on reciprocal best
  BLAST hits only:
  - YALI1F39330g: Rnr1-like large subunit, 67.5 % identity.
  - YALI1D04261g: a second large-subunit paralog with no S. cerevisiae reciprocal best hit.
  - YALI1B05072g: Rnr2-like small subunit, with no separate Rnr4 ortholog in W29.
  - YALI1F02309g: Trx2-like thioredoxin.
  - YALI1D35735g: Trr1, the only TRR in W29.
  - YALI1F12969g: Ura6.

  None of these names is experimentally verified in *Y. lipolytica*.
- ***S. cerevisiae*:** Rnr1 is cytoplasmic. Rnr2/Rnr4 are held in the nucleus by Wtm1 and Dif1 and move to the
  cytoplasm in S phase, after DNA damage or under low iron, to form the active holoenzyme (**verified**,
  abstracts). That the nuclear pool is itself inactive is the authors' model (**inferred**). No source found
  places RNR catalysis in the nucleus.
- **Mammals:** RNR is cytosolic (Pontarin 2008), and no reductase was found in nuclear preparations (Leeds
  1985). Pontarin suggests, but does not measure, that dNTPs diffuse into the nucleus. One human report recruits
  RNR to DNA-damage sites (Niida 2010).
- **Thioredoxin:**
  - It is the major in vivo RNR reductant (Camier 2007; Koc 2006).
  - Its nuclear location in yeast rests on Trx1 high-throughput GFP (SGD HDA), Trx2 computational (SGD IEA),
    and UniProt's citation of Huh 2003 (not opened).
  - It is small enough to cross the pore; this is inferred from Timney 2016.
  - No nuclear thioredoxin reductase is annotated in yeast.
  - The W29 RefSeq thioredoxin (`XP_504863.2`) has a 25-residue Arg/Lys-rich N-terminal extension that resembles
    a mitochondrial presequence (**inferred**; no predictor run). Which start is used in vivo is unknown.
- **Model history** (static scan of the local SBML files):
  - iYali4, iYli21, iYali5, Yeast-GEM 9.0.2 and iMM904 all carry nuclear RNR copies.
  - None of them has a nuclear source of reduced thioredoxin, so the copies are blocked in all of them. That
    follows from the stoichiometry; no LP was run on those models.
  - Yeast-GEM links only the disulphide (`r_4612`, a gap-fill with "Confidence Level: 0").

**What this means.** The mechanism gives no reason to choose: opening the module only adds capacity. The biology
gives a reason not to open it, because that would make the model claim a nuclear RNR activity no source
supports. Keeping the copies blocked also matches every ancestor model.

## 6. Could the method be improved, and how to test it

- **The census counts compartment copies with no source as gaps.** UM11 (and UM2 / UM7) are inherited
  compartment duplicates. A census flag such as "copy of an active reaction in another compartment, same GPR" would
  separate them from gaps that limit a function.
  - **Test:** for every blocked reaction, look for an active reaction with identical stoichiometry after
    stripping compartments and an identical GPR. Then check by hand that UM11's 4 RNR copies and R782 are
    flagged, and that R1104 (active) is not.
- **No-uptake loops.** FVA reports R789 and R1104 at ±1000 with no uptake, through a net-zero kinase /
  diffusion cycle. Capacity maps would be cleaner with loopless FVA (CycleFreeFlux or ll-FVA).
  - **Test:** run cobra's `loopless=True` FVA on the 25 mapped reactions. The loop should vanish while the
    cytosolic RNR capacity stays.
- **The knockout comparison** reused the Stage 1 scripts unchanged except for the in-memory variant.
  - **Test:** running `ko_compare_variant.py` with variant `as_built` should reproduce the baseline column
    exactly. That was not run.

## 7. What is proven, verified, observed or untested

- **Provable from stoichiometry:**
  - With no nuclear source of reduced TRX1, the four RNR copies must carry zero flux at steady state.
  - Closing reactions that are blocked in WT (C1) cannot change any knockout optimum.
- **Statically verified:**
  - the identical GPR strings of the nuclear and cytosolic copies;
  - the reaction inventory;
  - the copies in the five ancestor models;
  - the reciprocal-best-hit tables.
- **Observed in runs this session:**
  - all FVA, WT, pFBA and energy-check results;
  - the three knockout sweeps and the RNR controls;
  - the figure's LPs and the viewer checks.
- **Untested:**
  - any condition other than Task A + PO1f and "open";
  - nuclear RNR activity in *Y. lipolytica* (no data exist);
  - whether D04261g is an active large subunit;
  - the tooltips in an interactive browser.
- **Not reproduced by the audit:**
  - the seed-1 rebuild;
  - the WT values of A+ and C1;
  - the headless-Chrome viewer check.

  These are observed only in this session's runs.

## 8. Follow-ups found (not changed; each affects another line)

1. **RNR GPRs.** R706/R708/R710 (and the nuclear copies) use an OR of large and small subunits. R1913/R1914 use
   `F39330g or (B05072g and F02309g)`. The candidate form is (large) **and** (small).
2. **R702–R705.** These ribonucleoside-triphosphate reductases carry only the thioredoxin gene. They come from
   Yeast-GEM r_0970–0973.
3. **RNR false negatives.** Items 1 and 2 make the consensus-essential RNR genes YALI1F39330g and YALI1B05072g
   false negatives. On the baseline, every single RNR knockout grows fully. Only all three RNR genes knocked out
   with R702–R705 closed gives 0 (`stage1/rnr_controls.log`, verified).
4. **R1914's water** is written in the cytosol (inherited from iYali4).
5. **Thioredoxin formulas.** TRX1-disulphide has no formula, and TRX1 is `C6H9NO2S2R2`, so no thioredoxin
   reaction's mass balance can be checked.
6. **Duplicate thioredoxin pair.** `m2033` / `m2034` (used by R2198 and R2204) duplicate TRX1 / TRX1-disulphide.
7. **Nuclear dUMP-kinase loop.** R789, R1104, R1105 and R788 form a net-zero no-uptake cycle at ±1000.

## 9. Limits and accounting

- Only Task A + PO1f and "open" were tested.
- Stage 1:
  - 2 builds (the combined baseline at seeds 0 and 1);
  - about 35,600 LPs (stage1 1,307; whole-model FVA about 27,800; three sweeps 3 × 2,164; controls 20).
- Stage 2: 0 builds and 9,476 LPs for the figure, within the about 11,300 declared in TASK.md section 0.
- Large files (PubMed / UniProt / SGD fetches, `stage1_um11.json`, `blocked_sets.json`, the full knockout
  tables, the Escher work folder) are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um11_nuclear_rnr_20261029/`.
- The literature evidence was gathered by one subagent. The main session spot-checked five abstracts, the RBH
  table and the Yeast-GEM scan line, but not every citation. Huh 2003 (the GFP atlas) was not opened by anyone.
