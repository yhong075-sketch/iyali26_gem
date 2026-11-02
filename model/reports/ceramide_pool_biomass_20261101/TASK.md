# Ceramide pool in `biomass_C` (2026-11-01)

**Status: approved 2026-11-01** (§0). The baseline state of the targets is in [before.json](before.json), dumped
read-only by [scripts/before_state.py](scripts/before_state.py).

## 0. Authorization

This round comes from the blocked-in-WT false-negative discussion (session `claude/blocked-wt-fn-discussion-e329eb`,
2026-11-01). The user's approvals, in order:

1. "Approve steps 1–2": a literature search and FVA on the bucket-3 reactions.
2. "Authorize an in-memory test of the sphingolipid ceramide-pool term", then "Approve as written" for the test plan.
3. "Close the evidence gaps".
4. "Plan the curation round".
5. "Approve as written" for the plan below, at the recommended coefficient of **5.4e-4 mmol/gDW**.

No commit or push without a separate instruction.

## 1. Starting state (verified this session)

| Item | Value |
| --- | --- |
| Worktree / branch / HEAD | `claude/blocked-wt-fn-discussion-e329eb` @ `b83e13a`, the same commit as `claude/combined-r24-trehalose-20261031` (PR #42; corrected 2026-11-02: PR #42 was already merged to main as `5a922d2` on 2026-11-01 06:54 UTC, so "not merged" was stale when written); clean before this folder |
| **Baseline** | `candidates/E5_…_aco_split_treh_biomass_r24_aco_20261031.xml`, SHA-256 `c1f70ac4198e443ce2bc089f534c021d4111d8c2e999888e5beae06cc7087540`; the newest candidate, **not a designated reference** |
| Baseline build argv | taken from its `.build.json`: 41 flags ending `--trehalose-biomass-removal --r24-homocitrate-aconitase` |
| Experimental list | `consensus_essential_genes.csv`, SHA-256 `1e887f5a…` |
| Medium / profile | Task A `6eca70fd…` + PO1f `35307853…` |
| Baseline knockout screen | [ko_screen_combined.tsv](../combined_r24_trehalose_20261031/ko_screen_combined.tsv), Gurobi 1e-9, certified-infeasible rule B. At 15 %, mtDNA + untargeted out of scope: 136 TP / 190 FN / 43 FP / 695 TN |

## 2. Targets and their current state ([before.json](before.json))

- **`biomass_C`** (the objective, maximised, bounds (0, 1000)): it has no sphingolipid term. The legacy lipid
  pseudoreaction R1371, which used to drain "complex sphingolipid", is closed.
- **10 ER ceramide species**, all present with formulas:
  - m386 and m391 (ceramide-1, C24/C26, dihydro, NO3);
  - m388 and m392 (ceramide-2, phyto, NO4);
  - m390 and m393 (ceramide-2′, α-OH dihydro, NO4);
  - m397 and m398 (ceramide-3, NO5);
  - m399 and m400 (ceramide-4, NO6).
- **Nothing drains these ceramides.** FVA on 2026-11-01 (97 LPs) showed that the whole sphingoid backbone is blocked
  in WT: R722 (`YALI1C04265g and YALI1F20637g`), R37 (`YALI1B22883g`) and the SUR2 and SCS7 reactions all have
  min = max = 0.
- The pool metabolite and the `xPOOL_CER_*` reactions do not exist.

## 3. Change (one change)

- **New pseudo-metabolite:** `sphingolipid_ceramide_pool[C_er]`, with no formula and charge 0.
- **10 new gene-less pseudo-reactions** `xPOOL_CER_<mid>`: ceramide `<mid>[C_er]` → pool, bounds (0, 1000).
  - One per species, so any mix of species can fill the pool. Dihydro, phyto and α-OH forms all count.
  - Pseudo-reactions with no formula, unbalanced by design.
- **`biomass_C` consumes the pool at 5.4e-4 mmol/gDW**, plus curation notes. `biomass_C` is not re-normalised;
  nothing else changes.

## 4. Evidence (tiers as checked 2026-11-01)

- **Coefficient: [Rupcic 1996](https://europepmc.org/article/MED/8997745)**, *Candida (Yarrowia) lipolytica* on
  methanol.
  - Total long-chain bases are 0.7 % of polar lipids; polar lipids are 52.3 % of lipid; lipid is 4.9 % of DCW
    (verified, abstract).
  - Arithmetic: 0.007 × 0.523 × 0.049 = 1.794e-4 g LCB/gDW; divided by 331.54 g/mol (C19-phytosphingosine),
    5.41e-4 mmol/gDW.
  - One LCB per sphingolipid molecule, so this is moles of sphingolipid (inferred).
- **Cross-check: [Ianutsevich 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6950712/)**, W29, glycerol, stationary
  phase.
  - "A minority of SL" (verified). Fig. 4A, which I viewed, reads at about 1.3–2 % of 11.2–12.98 mg membrane lipid
    per g. That gives about 0.15–0.25 mg/gDW, the same mass level.
  - As intact lipid at roughly 700–1000 g/mol, that is about 1.5–3.7e-4 mmol/gDW (inferred). The method likely
    undercounts IPC/MIPC.
- **Need:**
  - In *S. cerevisiae*, lcb1, lcb2, tsc10 and aur1 nulls are inviable (SGD, verified by a review agent).
  - Myriocin, which blocks the first step, inhibits *Y. lipolytica* growth at first (verified, abstract).
  - No *Y. lipolytica* deletion data were found for LCB1, LCB2 or TSC10.
- **Species choice:**
  - The *Y. lipolytica* sur2Δ strain is viable but grows slowly
    ([Shin 2025](https://europepmc.org/article/MED/39603335) abstract and the same group's patent US 20240271081,
    verified). Its sphingolipids shift to dihydro forms.
  - *S. cerevisiae* scs7 is viable (SGD).
  - So the pool must accept dihydro and non-hydroxylated ceramides. The Shin 2025 full text is closed access and was
    not read.
- **Provisional, like the r_4598 cofactor terms.** Neither yeast-GEM 9.1.1 nor iYali has a sphingolipid biomass
  term (verified by a review agent).

## 5. Expected effect

From the 2026-11-01 test at 5e-4: a scratch SBML, 1,169 LPs.
- **WT growth:** about 0.31205, down about 0.07 %.
- **Forced flux:** only R722 and R37, at growth × coefficient.
- **Knockout classes:** exactly 7 genes become essential at every cutoff.
  - 6 true positives: `YALI1C04265g`, `YALI1F20637g`, `YALI1B22883g`, `YALI1B26350g` (ELO-like),
    `YALI1F15781g` (PHS1-like) and `YALI1A05116g` (TSC13-like).
  - 1 new false positive: `YALI1A06453g`, the IFA38/KAR ortholog (reciprocal best BLAST).
- **What stays non-essential:**
  - SUR2 and SCS7.
  - Ceramide synthase `YALI1F09542g`, kept open by the cytosolic R2127 `YALI1F09542g or YALI1D05038g`.
- **Expected score**, 15 %, rule B, mtDNA + untargeted out of scope: 142/184/44/694.

The score is a development reference with user-defined negatives. It is not the reason for the change.

## 6. Implementation

- **Data:** `model/curation/biomass/ceramide_pool_biomass.json`.
- **Step:** `platform/scripts/gem_annotate/ceramide_pool_biomass.py`, with `apply_ceramide_pool_biomass()`. It fails
  closed unless:
  - the objective is `biomass_C` alone, maximised;
  - `biomass_C` (stoichiometry, bounds and note keys) equals the c1f70ac4 state;
  - the 10 ceramides have the pinned name, formula, charge and compartment;
  - the pool and the `xPOOL_CER_*` reactions are either all absent, or present exactly as curated (a no-op).
- **Flag:** `--ceramide-pool-biomass`, run after `--trehalose-biomass-removal`, the last step that edits `biomass_C`.
- **Test:** `platform/tests/test_ceramide_pool_biomass.py`, in memory. It covers apply, idempotence, the
  fail-closed cases and context rollback.

## 7. Acceptance criteria

1. Full pytest has the same failures as before the change, and the new tests pass.
2. The control build without the flag reproduces `c1f70ac4` byte for byte. Seeds 0 and 1 of the new build are
   byte-identical.
3. The whole-model diff against `c1f70ac4` shows only the 1 metabolite, the 10 reactions, the `biomass_C` pool
   coefficient and the `biomass_C` notes.
4. WT is optimal and passes the checks.
5. In the energy-cycle check (no uptake, `xMAINTENANCE` lower bound relaxed), max ATP is 0.
6. In the certified knockout screen, exactly the 7 genes of §5 change class, and nothing else.
7. An independent read-only audit supports the claims.

## 8. Compute budget and stop conditions

- At most **1,300 LPs**: WT about 2, the energy check about 4, the certified knockout screen about 1,125.
  Task A + PO1f, Gurobi 1e-9, with the real HOME.
- **Stop and report** if:
  - the control build differs from c1f70ac4;
  - the diff shows anything not in §3;
  - WT is not optimal;
  - the energy check gives ATP > 0;
  - any gene beyond the 7 changes class;
  - the cap is reached.

## 9. Out of scope (recorded, not changed)

- The SUR2/SCS7 consensus labels, against the *Y. lipolytica* sur2Δ viability.
- The R2127 ceramide-synthase OR rule.
- The KAR false positive: no backup reductase is identified in *Y. lipolytica*.
- The elongation-cycle chemistry: R1414, R1420, R1407 and R1427 are charge-imbalanced; R1421 and R1428 are
  elementally imbalanced because m1449 has no formula.
- The complex-sphingolipid copies (census UM1; handed to the disconnected-pathway session).
- The colleague's `lipid_unlump/missing_lipids.csv` ceramide formulas, which are one O lower than c1f70ac4. That is
  the colleague's line: review only.
