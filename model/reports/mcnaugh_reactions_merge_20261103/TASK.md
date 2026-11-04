# TASK: add 17 reviewed mcnaugh_pass reactions and close the lumped steps they replace (2026-11-03)

## Authorization (user, this session, 2026-11-03)

- The user wants Andrew McNaughton's lipid unlumping merged: "yes, this is what important we wanna merge".
- Ownership: "Andrew owns it; go ahead".
- Scope: "All 56". Later narrowed to everything that is not a true duplicate or a reversal: "Approve, core screen".
- Biomass: on the two biomass reactions, "I think we can't skip", then "Both verbatim from Andrew".
- Chemistry:
  - "Copy his design, allow imbalance";
  - "Name it in the exception" for the C18:2 charge;
  - "Record it, change later" for the acyl-pool artefact.
- DGAT GPR: "Andrew's DGA1 or DGA2".
- Biology: "I think it's important to simulate the real biology".

## Baseline

- Model: `model/candidates/E5_..._ceramide_pool_bigg_annot_v2_20261103.xml`, SHA-256
  `9158e9f2a0d21f654f0a7919eed9f2bacba01aea8b315808d93ec14dd3af4129`, committed as `048e18ec`.
- Andrew's model: `iyali26_gem_lipid_unlump/mcnaugh_pass/iyali26.xml`, SHA-256 `e74a0f50a0a87700…`. It is read only.
- Python environment: `iyali26_gem/.venv`, Gurobi, checked solves at 1e-9.
- Medium and profile: Task A medium `po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) with PO1f `po1f_sd_leu.json`
  (`35307853…`).

## Review behind the scope

Andrew lists 56 reactions with no recorded counterpart in our model. Each was simulated on the medium above and
reviewed for biology by a reviewer, a skeptic and a tie-breaker.
- Duplicate review: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/mcnaugh_duplicate_biology_20261103/`.
- Review of the 29 to add: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/mcnaugh_reactions_review_20261103/`.

The outcome:

| Group | Reactions |
| --- | ---: |
| Duplicates of existing reactions, not added | 26 |
| Held (cytosolic COQ2, CoQ9 gate) | 1 |
| Not added: duplicates found by the second review | 10 |
| Not added: would reverse an authorized curation (ALCD2irm re-adds R132, G3PD1irm re-adds R349) | 2 |
| **Added** | **17** |

The 10 duplicates found by the second review are 4ABUTtm, ABUTt2r, 3HPCOAHYD, COMT2, COMT3, PAC2HYD and
DIPEP1-4_v.

## Targets

Curated data: `model/curation/mcnaugh_pass/mcnaugh_reactions_20261103.json`. The starting state is in `before.json`.

**Add, with Andrew's bounds and GPRs:**
- G3PAT160/161/180/181/182_lp, GPR YALI1C00230g (SCT1);
- AGPAT160-182_lp, GPR YALI1E22736g or YALI1F25951g;
- DGAT160-182_lp, GPR YALI1E38810g or YALI1D10264g;
- BIOMASS_C_core and BIOMASS_C_full, translated term for term, bounds 0..1000, objective 0.

All species and genes already exist in our model. The lipid-particle mapping is 1agp_lp = m577, pa_lp = m1705,
dag_lp = m1642 and tag_lp = m1640.

**Close** (bounds 0..0; the reactions stay):
- R353, R1846, R1771: the lumped Kennedy steps with generic acyl-CoA m570;
- R2165: cytosolic lumped TAG synthesis.

R1771 must close in the same change: left open next to the new DGATs, it loops and WT growth rises to 0.5127.

**Keep:** R351 (DHAP acyltransferase). Andrew's design has no counterpart for it, so this departs from his design and
goes on the list for him.

**Exceptions recorded in the reaction notes:**
- the R-group mass imbalance of the 15 lipid steps;
- a +4 charge on the three C18:2 steps, from the pre-existing m1735[C_lp] defect.

## Expected effect

- **Changes:** only the 17 new reactions and the bounds of the 4 closed ones. Species, genes, other reactions and the
  objective stay the same.
- **WT growth:** about 0.3124886, +0.139%. The rise is a recorded artefact of reversible acyl pools; it stays under
  the 0.5% gate.
- **Energy check:** ATP with no uptake = 0.
- **Knockout screen:** will change for genes on the new and closed reactions (SCT1, SLC1, YALI1F25951g, DGA1, DGA2,
  TGL3, TGL4).

## Acceptance criteria

1. The guarded `apply_mcnaugh_reactions`:
   - adds exactly the 17 and closes exactly the 4;
   - is idempotent;
   - fails closed without edits on a stale, partial, rewritten or out-of-scope state;
   - survives an SBML round trip.
2. A whole-model diff of v3 against v2 shows only those changes.
3. Balance: of the new reactions, only the 15 lipid steps (and the biomass pseudo-reactions) are unbalanced, as
   recorded.
4. FVA shows that every new reaction can carry flux.
5. WT growth and the energy check are as expected.
6. The certified knockout screen (biomass_C) is compared gene by gene with v2. The screen with BIOMASS_C_core as the
   objective is reported for information.
7. pytest: only the new tests differ from the before run.
8. An independent audit is run.

## Budget and stop conditions

**Budget:** one build, and a spare. Up to about 30 LPs for the diff, FVA, WT and energy checks, plus two certified
screens of about 1,100 LPs each.

**Stop if:**
- a guard trips on the baseline;
- the diff shows a change outside the targets;
- WT growth goes up by more than 0.5% or collapses;
- the energy check opens a cycle.

## Coordination

Three other sessions build from the same baseline in their own worktrees, and each touches `main.py`/`cli.py`:
- ribonucleotide reductase;
- folate/polyamine/glycerate;
- annotation corrections.

The compartment-name and gap-fill branches also exist. Combining any of them later needs one rebuild with all the
flags.
