# UM7: nucleus-cytosol diffusion of ADP-ribose and diphosphate, candidate 5f54eb22 (2026-10-28)

## 1. Result

- **Candidate:**
  [`model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml`](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml),
  SHA-256 `5f54eb22086b06b2bd5e8715b89984060cd43c1a7fa21e33e5581d867ee54179` = `9d655c2b` (dolichol route, the
  newest candidate on `main` after PR #29) + `--um7-nuclear-diffusion`. Not a designated reference.
  Committed in `3720370` on branch `claude/hungry-mendel-f9f7df` (2026-10-28), not pushed.
- **Change:** two gene-less, reversible reactions in the "Transport, nucleus" group (SBO:0000185), the same kind
  as the 33 nuclear diffusions the model already has:
  - **R_ADPRtn** "ADP-ribose diffusion, nucleus": ADP-ribose `m1725[C_cy]` ⇌ `m818[C_nu]`;
  - **R_PPItn** "diphosphate diffusion, nucleus": diphosphate `m203[C_cy]` ⇌ `m210[C_nu]`.
- **Authorization** (TASK.md sections 0 and 11), from the user on 2026-10-28:
  - "Approve Option A", for the Stage 1 draft;
  - "Accept, restate gate", after the run stopped at the knockout gate ([STOP.md](STOP.md)).
- **Effect:**
  - The inherited nuclear NAD block becomes connected: nuclear NAD synthesis R569 / R579 → sirtuin step R567
    → ADP-ribose to the cytosolic pyrophosphatase R1912, and PPi to the cytosolic pyrophosphatase R419.
  - Blocked UM7 reactions (whole-model FVA at 0 %, Gurobi 1e-9): 15 → 6 with every exchange open, 15 → 8 on
    Task A. Outside UM7, R1912 (census UM171) is freed in both conditions and R576 (nicotinamidase) on Task A.
  - Nothing becomes blocked. WT growth is unchanged, with no knockout class change and no energy-generating
    cycle.
  - One knockout gains growth: `YALI1F37714g`, +1.2e-4 absolute (section 5).

## 2. Why (TASK.md section 3; tiers as given there)

- **Structural.** The census roots were nuclear ADP-ribose (RNC) and nuclear PPi (RNC), next to cytosolic
  ADP-ribose (RNP, UM171). This is the "same compound, dead end on each side of a membrane" signature of a
  missing transport. The reactions came unchanged from iYli21 (verified, TASK.md section 2). yeast-GEM has the
  same block and adds the gene-less ADP-ribose transport r_4603 (verified in the local yeast-GEM copy; model
  precedent, not evidence).
- **Mechanism, passive diffusion through the nuclear pore.** Timney et al. 2016 (in vivo, *S. cerevisiae*)
  define passive diffusion as the equilibration of "ions, small molecules, and small- to medium-sized
  macromolecules" between cytoplasm and nucleoplasm (verified). Applying this to ADP-ribose and PPi is
  **inferred**.
- **Nuclear NAD metabolism is real in yeast** (*S. cerevisiae*, verified):
  - Sir2 and Hst4 are nuclear (UniProt P06700, P53688, ECO:0000269). Hst1 is nuclear by a large-scale GFP
    screen only (P53685, which also notes cytoplasm). TASK.md section 3 cites P06700 alone for all three
    (corrected, TASK.md section 12).
  - Nma2 and Npt1 are concentrated in the nucleus; Qns1 is nuclear and cytosolic.
  - Ipp1 is nucleocytoplasmic.

  The *Y. lipolytica* orthologs on these reactions are reciprocal best BLAST hits (this session); their
  localization is **predicted only**. No *Y. lipolytica* localization data exist for any of these enzymes.
- **Why only these two links** (TASK.md sections 5–7):
  - Either link alone frees nothing.
  - An NAD link (Option A′) frees nothing more.
  - The homocitrate branch would need an invented mitochondrial carrier, and *Y. lipolytica* LYS1 is only
    predicted mitochondrial.
  - R92 stays closed: no reason for its closure is recorded, and opening it would only duplicate R90.

## 3. Implementation

- **Data:** [`model/curation/nuclear_transport/um7_nuclear_diffusion.json`](../../curation/nuclear_transport/um7_nuclear_diffusion.json)
  (`UM7-NUCLEAR-DIFFUSION-20261028`). It holds the two reactions with notes, evidence rows and tiers, and the
  reviewed identity of the four species.
- **Step:** [`platform/scripts/gem_annotate/um7_nuclear_diffusion.py`](../../../platform/scripts/gem_annotate/um7_nuclear_diffusion.py),
  `apply_um7_nuclear_diffusion`. It fails closed if:
  - the scope or status differs;
  - a species is missing or differs in name, formula, charge or compartment;
  - a dead-end end has gained a reaction (`reactions_before` on m818, m1725 and m210). This also rules out an
    existing link between a pair;
  - the group is missing or not unique;
  - the two ends of a link are not one compound in two compartments;
  - a reaction id exists with another definition (an id that already exists with *this* definition counts as
    applied, so a rerun is a no-op; TASK.md section 8 had said any existing id fails);
  - only one of the two reactions is present.

  A rerun, also on a saved and reloaded model, is a no-op.
- **Flag:** `--um7-nuclear-diffusion`, opt-in. It requires `--dolichol-route` and runs after it (`main.py`,
  `cli.py`).
- **Test:** [`platform/tests/test_um7_nuclear_diffusion.py`](../../../platform/tests/test_um7_nuclear_diffusion.py),
  6 tests and 24 subtests: spec vs task, 12 stale inputs, partial states, scope, the exact two edits with an
  SBML roundtrip, and the flag. It reads the baseline through a `Path`, because cobra 0.30 rejects str SBML
  paths of 260 characters or more on macOS (its Windows test, `"win" in platform`, matches `darwin`).
- **One implementation change against the draft:** section 8 planned a separate "a transport between the pair
  already exists" guard. It was dropped as unreachable: every pair has a nuclear end with a pinned reaction
  context, which catches any other link first. The test checks that case through the context guard.
- **Build command** (from `platform/`, miniconda Python 3.13.5, cobra 0.30.0):

  ```bash
  PYTHONHASHSEED=0 python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" --offline --no-solve --coq9-curation qcycle --energy-candidate E5 --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh --restore-selection-cleared-gprs --restore-selection-overwritten-gprs --restore-selection-cpa-ura2 --remove-selection-tkl-fragment-gpr --remove-r153-peptidase-gpr --r2041-forward --r2190-ribulose-product --r1964-gpr --gpi-chain-links --gpi-coproducts --tpp-requirement --thiamine-phosphate-chemistry --thiamine-proton-lock --heme-iron-biomass --r2193-iron-gpr --lipoylation-module --acyl-acp-chemistry --cofactor-biomass-terms --r2135-acyl-coa-pool --r1726-cpt-open --dolichol-route --um7-nuclear-diffusion --output-model ../model/candidates/E5_…_dolichol_route_v2_um7_nuclear_diffusion_20261028.xml
  ```

  The argv is the `9d655c2b` `.build.json` argv plus the new flag ([logs/build_argv.txt](logs/build_argv.txt)).

## 4. Verification (all this session)

| Check | Result |
| --- | --- |
| Builds | The control without the flag gives `9d655c2b…` byte for byte. `PYTHONHASHSEED=0` and `=1` both give `5f54eb22…` ([logs/sha256.txt](logs/sha256.txt)). 3 builds, the planned number. |
| Whole-model diff vs `9d655c2b` | 3 entries: R_ADPRtn added, R_PPItn added, group101 "Transport, nucleus" members ([model_diff_vs_9d655c2b.json](model_diff_vs_9d655c2b.json)). No species, gene, bound, GPR, formula, charge, annotation or note changed elsewhere. |
| Balance | Both reactions are mass and charge balanced (`check_mass_balance` = {}; ends are C15H21N5O14P2 −2 and H4O7P2 0). |
| WT (Task A + PO1f, Gurobi 1e-9, checked) | 0.3116608216 in both, +4.6e-15 relative. pFBA carries 0 through every blocked UM7 reaction and both links at the optimum; the already-active R568 / R580 carry 0.0010 in both models. |
| Energy-generating cycle | Every exchange (0, 1000), xMAINTENANCE 7.8625 relaxed: max R_NTP1 = 0 and max biomass = 0, both **optimal**, both models ([built_checks.json](built_checks.json)). |
| Whole-model FVA blocked sets (0 %, ≤ 1e-9) | Open condition 705 → 695; Task A 837 → 828. **Nothing newly blocked.** Freed (open): R567, R569, R579, R1017, R1082, R1212, R1214, R1267, R1270, R1912. Freed (Task A): R567, R569, R1017, R1082, R1212, R1214, R1267, R1912, R576. Exactly the Stage 1 prediction. |
| Knockout sweep (2 × 1,081 genes, as `PO1f_plasmid_LEU2` is excluded at runtime, + 2 WT = 2,164 LPs; checked, reset per knockout) | No class change at 1 / 5 / 10 / 15 %. At 15 %, TP / FN / FP / TN = 112 / 203 / 42 / 692 in both, with 16 unresolved (the same infeasible knockouts) and 326 in-model positives. No point failed the stale-point check; max knocked flux 0. One knockout changes by more than 1e-6 (section 5). Negatives are the user-defined class (model genes not on the list); the list was used in development, so this is a development reference. ([ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv)) |
| pytest | 314 → 320 passed; the same 19 failed and 14 errors before and after; the only difference is the 6 new tests ([logs/pytest_compare.txt](logs/pytest_compare.txt)). Caveat: the before run started before the `main.py` / `cli.py` edits but ran while they were made. |

Full scores at all cutoffs, the same in both models:

| Cutoff | TP | FN | FP | TN | Recall | Specificity | Precision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 % | 94 | 221 | 35 | 699 | 0.298 | 0.952 | 0.729 |
| 5 % | 100 | 215 | 37 | 697 | 0.318 | 0.950 | 0.730 |
| 10 % | 104 | 211 | 42 | 692 | 0.330 | 0.943 | 0.712 |
| **15 %** | **112** | **203** | **42** | **692** | **0.356** | **0.943** | **0.727** |

Out-of-scope positives (list genes with no model counterpart): 1,612 − 326 = 1,286, unchanged.

## 5. What it means

- **Mechanism vs biology.** The two links are what the structure needs: an ADP-ribose and a PPi outlet. The
  mechanism is an established general rule (NPC passive diffusion). The biology (nuclear sirtuins and NAD
  synthesis) is experimental in *S. cerevisiae* and inferred for *Y. lipolytica*.
- **The optimum does not change.** The links widen the flux space only, because no WT optimum needs nuclear NAD
  turnover.
- **The PNP knockout gain** ([STOP.md](STOP.md), verified):
  - In the knockout of `YALI1F37714g` (purine-nucleoside phosphorylase by KEGG EC 2.4.2.1; name not verified;
    13 reactions), growth goes from 0.3114518 to 0.3115725, KO/WT 0.99933 → 0.99972.
  - The route: spare tryptophan → kynurenine pathway → NaMN → deamido-NAD → nuclear NAD (R569) → sirtuin
    hydrolysis (R567) → ADP-ribose → R1912 → AMP + ribose-5-P, with nicotinamide recycled by R576. All of
    this runs at 0.0036.
  - Closing R_ADPRtn, R_PPItn, R567 or R1912 each restores the baseline knockout value exactly.
- **Why that happens (inferred).** R567 is a lumped "NAD + H2O → ADP-ribose + nicotinamide" step with no
  acetylated substrate. Once its products have outlets, it acts as an unconstrained NAD hydrolase. In yeast,
  sirtuin NAD use is coupled 1:1 to deacetylation, which makes O-acetyl-ADP-ribose (Tanner 2000). The user
  accepted the gain as a documented consequence ("Accept, restate gate"), and constraining R567 is a separate
  question.
- **Essentiality.** No predicted class changes. `YALI1A21156g` (QNS1 ortholog, on the list) stays essential,
  because its cytosolic R568 is required. The sirtuin, NMNAT and Ysa1 genes stay non-essential.

## 6. Limits

- **Evidence:** *S. cerevisiae* localization and a general nuclear-pore rule. No *Y. lipolytica* localization
  data; orthology by reciprocal BLAST only. Diffusion of ADP-ribose and PPi specifically is inferred.
- **R567:** its chemistry is a stand-in for sirtuin deacetylation and is not changed. The PNP knockout gain
  comes from that.
- **Still blocked:**
  - R396 and its transports (homocitrate route unresolved);
  - R92 (closed, reason unknown);
  - R579 / R1270 / R578 on Task A (no NMN source on that medium).
- **Pre-existing build message:** all three builds log a libSBML error while setting the notes on R304. The
  control build, which reproduces `9d655c2b` byte for byte, logs it too, so it predates this change and is
  not caused by it.
- **Cross-references:** the existing nuclear transports carry BiGG / MetaNetX ids; R_ADPRtn and R_PPItn carry
  only SBO:0000185. Adding ids was outside the approved scope; yeast-GEM's r_4603 is the obvious one for
  R_ADPRtn.
- **Not tested:** media other than Task A; the PNP knockout's real phenotype; a whole-model GLPK census rerun
  (the blocked sets here are Gurobi 1e-9 FVA, so the counts are not comparable to the GLPK census numbers in
  STATE.md, for example 708 vs 705).
- **Budget:** 3 builds against 3 planned. **About 60,000 LPs against about 20,000 planned (TASK.md section 8),
  an overrun of about 3×.** The whole-model FVA scans took most of them:
  - Stage 1 option scan: 4 variants × 2 conditions × about 4,640 ≈ 37,000;
  - built checks: 2 models × 2 conditions ≈ 18,500;
  - Stage 1 targeted LPs 2,142; sweep 2,164; figure 272; diagnosis about 20.

  These are cheap Gurobi LPs (minutes each scan). The planned figure had counted one whole-model scan per model,
  not two conditions per variant.

## 7. Next gates (each needs its own authorization)

- Constrain R567 (cap it, or rewrite it as deacetylation with an acetyl-protein source), if the user wants
  the sirtuin step to stop acting as a free NAD sink.
- The homocitrate branch (R396), if a *Y. lipolytica* localization or carrier becomes known.
- The HST4 ortholog `YALI1C31830g` (not in the model), a GPR question for R567.
- The next census target: UM2 vacuolar closures.

## 8. Audit

[AUDIT.md](AUDIT.md): an independent read-only reviewer, with its own SHA, diff and test checks and about 60 of
its own LPs on byte-identical copies.

- **Supported (9):** C1–C8, C11, C12. The reviewer reproduced the SHAs, the raw XML diff, WT, the EGC, 12 spot
  FVA checks, the `YALI1F37714g` knockout in both models and all four STOP.md controls.
- **Partially supported (2):**
  - C9: the user's wording cannot be seen by the reviewer.
  - C10: two citation weak spots, the Hst1 / Hst4 localization and the source of the Ysa1 IDA.
- **Unsupported:** none.

Its other findings were fixed in the text (TASK.md section 12, STOP.md step 4, this report, the Escher README)
and changed no model byte:
- stale paths and names in the Stage 1 draft;
- the cause of the long-path failure (cobra, not libSBML);
- "no pFBA flux through any UM7 reaction" (R568 / R580 carry 0.0010);
- the sweep is 2 × 1,081 genes, not 2 × 1,082;
- the guard accepts an identical existing reaction;
- the ATP import in the knockout runs through R1021, not R1017.

Two findings are recorded as limits (section 6): the pre-existing R304 notes error, and the missing
cross-references on the new reactions.

## 9. Escher figure

[escher/README.md](escher/README.md):
- `um7_nuclear_before_after.png` / `.svg`: capacity before | after, i.e. max |v| at FVA 0 % on Task A. pFBA
  is 0 on these reactions at the optimum in both models.
- `um7_nuclear_status.png`: the status view.
- `um7_nuclear_escher.html`: the self-contained viewer with five scenarios.

On the 23 drawn reactions, blocked on Task A goes 18 → 9 and flux-capable 3 → 14. The viewer was checked in
jsdom and in headless Chrome.
