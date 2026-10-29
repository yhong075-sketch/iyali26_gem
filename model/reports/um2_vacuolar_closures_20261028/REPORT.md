# UM2: the inherited vacuolar closures are kept (decision A0, 2026-10-28)

## 1. Result

- **Decision:** the user chose **Option A0** on 2026-10-28, after the Stage 1 draft
  ([TASK.md](TASK.md) section 0):
  - all 17 UM2 closures stay (0, 0);
  - the 19 dipeptide reactions and the open Ile import R882 stay as they are;
  - the decision and its reasons are recorded here, in `model/STATE.md` and on the UM2 entry of the census
    report.
- **No model change:** no build, no curated data and no builder code. Every candidate file is unchanged.
- **Model the conclusions apply to:**
  - Stage 1 ran on `9d655c2b…` (dolichol route v2).
  - UM7 was merged during the task (PR #30). The newest candidate on `main` is now `5f54eb22…`
    (`…_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml`). It differs from `9d655c2b` by R_ADPRtn and
    R_PPItn only.
  - All 69 vacuolar and UM2 reactions and every vacuolar species are identical in the two files (checked
    2026-10-28). The as-built, C′ and "open all 17" checks were re-run on `5f54eb22` and reproduce
    (`stage1/stage1_um2_5f54eb22.log`).
  - Not a designated reference.
- **Condition:** Task A medium + PO1f, Gurobi 13 at 1e-9, real HOME. The 15 % KO/WT cutoff is primary;
  1, 5 and 10 % are also reported.
- **Figure:** [escher/README.md](escher/README.md). It shows as built vs what-if C′ (not implemented), the
  status map, the loop of the census test, and an offline viewer checked in headless Chrome.

## 2. Why the closures stay

Evidence tiers follow TASK.md: **verified** means read on 2026-10-28; **inferred** means reasoning only.

**History (TASK.md section 2).**
- The 17 reactions are open in the *S. cerevisiae* Yeast 7.11 template and closed in the published iYali4
  (2016). Every later Y. lipolytica model in the line kept them closed: iYali4_corr, iYli21, iYali26 and
  iYali5. The Chalmers rebuild iYali v4.1.2 and Yeast-GEM have them open (read this session).
- No record says who closed them or why. The 2026-09-11 provenance task found no author statement, and its
  independent audit left only the negative search open.
- A plausible motive, **inferred**: in Yeast-GEM, the reversible vacuolar H+ diffusion and the H+-coupled
  vacuolar amino-acid transports run internal loops with every uptake closed (verified, run this session,
  `stage1/yeastgem_loops.log`). Closing the set removes them.

**Mechanism (TASK.md sections 4 and 6; computed on the model).**
- Every vacuolar amino acid has only an import and an export of itself. The polyamines have only an import.
  Vacuolar H+ is used only for transport.
- So at steady state, with no storage, pH or membrane-potential term, opening any subset gives one of two
  things:
  - an ATP → pump → H+ → shuttle/leak futile cycle. On Task A, R795 capacity is [0, 55.5] at 0 % of the
    optimum and pFBA flux is 0;
  - with R1163 in its written direction (cytosol → vacuole), an unbounded no-uptake loop: R1163 at 1000,
    ten amino-acid transports at 500.
- Under every option tested, WT growth is unchanged to within 3e-15. No single-gene knockout changes by
  more than 1e-6, and there is no class change at any of the four cutoffs. No knockout gains growth. The
  energy check is 0 everywhere (status optimal, `xMAINTENANCE` relaxed).
- 15 % scores, unchanged: TP 112, FN 203, FP 42, TN 692, 16 unresolved. Negatives are user-defined (model
  genes not on the consensus list; mtDNA genes excluded). These are a development reference, not
  independent validation.

**Biology (TASK.md section 5; homology only, no W29 name verified).**
- The vacuolar V-ATPase and the AVT-family transporters are real. W29 has reciprocal best hits for AVT1
  (`YALI1F37658g`), AVT3 (`YALI1E13514g`), AVT6 (`YALI1D14670g`) and ATG22 (`YALI1C04798g`).
- *Y. lipolytica* vacuoles pump protons (Kulakovskaya 1993), and the V-ATPase core genes are essential in
  the PO1f screens.
- So the closure does **not** claim that these reactions are absent from the organism. It records that this
  steady-state model has no function for them.
- The 11 V-ATPase false negatives come from a missing acidification requirement, not from the closure:
  opening R795 does not change them.

**Why not B, B′ or C′.**
- These options are biologically defensible and harmless, but they add no predictive content.
- Their only measurable effect is a lower census blocked count (Task A at best 37 → 25). That count is not a
  reason to change the model (`AGENTS.md`).
- They would also make FVA show non-zero ranges whose only mode is futile.
- The full option analysis and acceptance criteria are kept in TASK.md sections 6–7, in case the question
  is reopened. For example, if the model gains a vacuolar storage or pH term, C′ is the loop-free version.

## 3. Defects found in the closed or dead reactions (not changed, follow-ups)

1. **R925 spermidine transport is unbalanced:** −3 H and −3 charge. Its stoichiometry comes from the
   2026-09-10 metadata selection; Yeast-GEM `r_1261` is 1:1 and balanced.
2. **R870 is named "L-aspartate transport" but carries asparagine.** Inherited; the same name sits on
   Yeast-GEM `r_1188`.
3. **R1163 is written in the uphill direction** for an acidic vacuole. Yeast vacuole pH is 6.2 (Preston
   1989) and the passive leak runs vacuole → cytosol (Tarsio 2011).
4. **The dipeptide sets cannot work.** They were added by hand in iYli21 (Supplementary Table 1, verified),
   apparently for W29's Biolog PM1 carbon-source growth on all four dipeptides (Fig. 1, verified).
   - Their exchanges allow secretion only.
   - The hydrolases are gene-less, sit in the vacuole and are named "cytosol".
   - The dipeptides have no formula.
   - The yeast route is PTR2 uptake plus cytosolic hydrolysis. `YALI1D24836g` is the PTR2 reciprocal best
     hit.
5. **R2055 / R2056 (vacuolar His / Lys import) are open but blocked:** there is no exit from the vacuole.
   - The VBA1 reciprocal best hit `YALI1B04391g` (essential in the Cas12a screen) is not on their GPR.
   - The reactions are written as uniport, while the yeast mechanism is H+ antiport (Sato 1984).

## 4. What is verified, reproduced, untested

- **Verified (2026-10-28):**
  - the closure history in the three external model files read this session, and the Yeast-GEM loop;
  - all Stage 1 LPs: 5,811 + 447 (C′) + 1,341 (on `5f54eb22`) + 6,492 (knockout sweeps), all WT points
    passing the targeted check, 0 fresh re-solves;
  - the figure's 14,388 LPs;
  - the R925 imbalance;
  - the two literature spot checks (Russnak 2001, Preston 1989).
- **Reproduced:**
  - the census UM2 numbers on `9d655c2b` and `5f54eb22` (37/37 blocked in both conditions, C′ 3/37 open
    and 25/37 Task A);
  - the baseline 15 % scores, which match the dolichol report.
- **From the evidence subagent, not re-read by me:**
  - the remaining abstracts and UniProt entries;
  - the BLAST reciprocal best hits;
  - the PO1f screen calls.

  Two of my PMID links were wrong in the draft and were corrected against PubMed; the subagent's
  `stage1/evidence.md` is kept as written.
- **Untested:**
  - any vacuolar storage, pH or nitrogen-starvation condition;
  - media other than Task A and "every exchange open";
  - the motive of the 2016 closure, which remains inferred.
- **Independent read-only audit:** see [AUDIT.md](AUDIT.md).
  - The auditor re-derived claims 1–16 with its own scripts; all are supported, with no corrections.
  - It was cut off twice by the computer sleeping before writing AUDIT.md, so the main session compiled the
    file from its logs.
  - Not audited: wording, the subagent biology (BLAST, screen calls) and the Escher counts.

## 5. Files

- [TASK.md](TASK.md): decision, Stage 1 draft (history, inventory, trace, biology, options, criteria).
- `stage1/`: logs, `ko_um2.tsv` / `ko_um2_summary.json`, `evidence.md`, BLAST tables.
- `scripts/`: `stage1_um2.py`, `ko_compare_variants.py`, `inventory.py`, `trace.py`, `ext_models.py`,
  `yeastgem_loops.py`.
- `escher/`: the figure, viewer, scripts and logs.
- `$IYALI26_RESEARCH_ROOT/artifacts/tasks/um2_vacuolar_closures_20261028/`: the full Stage 1 JSON, the
  BLAST databases, the PubMed abstracts and iYli21 files (`stage1/`), and the Escher work folder (`escher/`).

**How the code and the method could be improved.**
- The knockout comparison applies option bounds in memory. A built candidate would need its own sweep, but
  none was built.
- "Freed" in this report means FVA capacity at 0 % of the optimum. A loop-aware FVA (loopless FVA or a
  thermodynamic constraint) would show directly that the C′ capacity is only futile. Test it by running
  cobra's `loopless` FVA on the C′ variant and checking that R795 keeps its range while every no-uptake
  range is 0.
