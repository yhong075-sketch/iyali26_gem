# TASK: THF, FAD, riboflavin, NAD, NADP and CoA in biomass_C at the yeast-GEM values (2026-10-27)

## 1. Authorization

User, 2026-10-27, after reading the draft of this plan (the message as received ends mid-sentence, at the
300-character limit of the app's reply button; the TPP and lipoylation entries of STATE.md record the same
situation):

> Approved: implement the cofactor biomass-terms plan as written: --cofactor-biomass-terms on top of
> a1237821 (running after --acyl-acp-chemistry), THF m479, FAD m528, riboflavin m926, NAD m122, NADP m41
> and CoA m69 in biomass_C at the yeast-GEM 9.1.1 r_4598 values 6.34e-5, 1e-5, 9.9e-4, 2.65e-3, 5.7e

"The plan as written" is the draft of 2026-10-27 (shared notes folder, `cofactor_biomass_terms_TASK_draft_20261027.md`);
sections 2–10 below are that draft unchanged apart from the status line and the paths of this folder. The
draft's full sentence continued: "… 5.7e-4 and 1.9e-4 mmol/gDW; species identity check; phase-1 in-memory
supply checks with a stop on any failure; control and seed builds; whole-model diff; Gurobi 1e-9 two-model
knockout sweep with the knocked-flux check; blocked-set re-census; code-reviewer and independent audit;
REPORT.md and STATE.md entry. No commit or push."

**Why this change, in the paper's words.** The Ponce-de-León 2013 census of a1237821
([REPORT.md](../gap_metabolites_ponce2013_20261026/REPORT.md), section 5) found 72 blocked reactions that
are *coenzyme pseudo-gaps*: the cofactor is synthesised but never drained, so its biosynthesis cannot carry
flux although the cofactor itself is not a gap metabolite. The paper's remedy is a biomass term (or drain)
for the cofactor. `biomass_C` already does this for TDP (1.2e-6) and heme a (1e-6), both at yeast-GEM
values (STATE.md 2026-10-22). This task extends the same precedent to the six cofactors that yeast-GEM's
cofactor pseudoreaction prices and the census tested. The 16 false-negative genes on the affected reactions
are the essentiality question this change touches; they are not its justification (AGENTS.md section 2: no
change is made to improve a match rate).

## 2. Starting state (read 2026-10-27)

- **Branch / HEAD:** `claude/gap-metabolites-model-2c880a` at `aeba460` (main after PR #22).
  Worktree not clean at task start: `model/STATE.md` modified and
  `model/reports/gap_metabolites_ponce2013_20261026/` untracked (the census session's work; preserved).
- **Baseline model:**
  `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_20261023.xml`,
  SHA-256 `a12378210c7c66720d99fb7bf71d14cb7538ea6a3cf3e08919c85fb3c9c47d82`
  (re-hashed 2026-10-27). 2,316 reactions, 1,882 metabolites, 1,081 genes. Built 2026-10-23 at HEAD
  `4c4248d8` with the flags recorded in its `.build.json` (the last two are `--lipoylation-module
  --acyl-acp-chemistry`). It is the newest candidate on `main`, **not** a designated reference.
- **Reference model:** none designated (`model/candidates/README.md`).
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`
  (`IYALI26_RESEARCH_ROOT`; the build fails closed without it).
- **`biomass_C` of a1237821 (read with cobra):** 55 metabolites; the small terms are CoQ9 1e-4,
  lipoyl protein 1e-5, lipoyl-E2 1e-5, iron 3.04e-5, TDP 1.2e-6, heme a 1e-6. **None** of THF, FAD,
  FMN, riboflavin, NAD(P), CoA, PLP, dolichyl-P, siroheme, biotin or glutathione is present. The closed
  legacy biomass R1387 / R1710 carried riboflavin at 0.00099 (the yeast-GEM value).
- **before.json:** [before.json](before.json), written by [scripts/before_state.py](scripts/before_state.py)
  (biomass_C stoichiometry, notes, bounds; the six species with formula, charge, annotations and reaction
  lists; the route reactions of section 4).

## 3. Targets

| Term | Model species | Name / formula / charge in a1237821 | Cross-references in the model | yeast-GEM species | Coefficient (mmol/gDW) |
| --- | --- | --- | --- | --- | ---: |
| THF | `m479[C_cy]` | THF_C19H23N7O6, charge 0 | MNXM1128450, CHEBI:166569 | s_1487 THF [c] | 6.34e-5 |
| FAD | `m528[C_cy]` | FAD_C27H33N9O15P2, charge 0 | KEGG C00016, BiGG fad, MNXM1364149 | s_0687 FAD [c] | 1e-5 |
| riboflavin | `m926[C_cy]` | riboflavin_C17H20N4O6, charge 0 | KEGG C00255, BiGG ribflv, MNXM270 | s_1405 riboflavin [c] | 9.9e-4 |
| NAD | `m122[C_cy]` | NAD_C21H27N7O14P2, charge 0 | MNXM64096, CHEBI:13389 | s_1198 NAD [c] | 2.65e-3 |
| NADP | `m41[C_cy]` | NADP(+)_C21H28N7O17P3; formula field C21H29N7O17P3, charge 0 | KEGG C00006, BiGG nadp, MNXM5 | s_1207 NADP(+) [c] | 5.7e-4 |
| CoA | `m69[C_cy]` | coenzyme A_C21H36N7O16P3S, charge 0 | MNXM1094981, CHEBI:23354 | s_0529 coenzyme A [c] | 1.9e-4 |

- Exact yeast-GEM values: THF 6.34000025456771e-05, FAD 9.99999974737875e-06, riboflavin
  0.000989999971352518, NAD 0.00264999992214143, NADP 0.000569999974686652, CoA 0.000190000006114133.
  The rounded values above are used, as the TPP step used 1.2e-6.
- `biomass_C` gains six consumed terms and `cofactor_biomass_terms_*` notes. No reaction, metabolite,
  gene, GPR or bound changes. Cytosolic species only, as in yeast-GEM; mitochondrial and other copies
  are not drained (so the census's 51 is an upper bound, see section 7).
- Not included although in r_4598: NADH 1.5e-4, NADPH 2.7e-3, glutathione 1e-6 (the user named six
  cofactors; the glutathione drain freed no reaction in the census). FMN, PLP, dolichyl-P, siroheme,
  biotin and pantetheine-4'-P have no yeast-GEM value and stay out.

## 4. Facts established before any change (2026-10-27, reads only)

**Coefficients (verified).** `https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/main/model/yeast-GEM.yml`
fetched 2026-10-27; SHA-256 `2ddfaaa73ead41ab243127d1b34fa56db2c86052f6eabcd12bb4eea93947dc17`, identical
to the file read on 2026-10-22 for the TPP step; version 9.1.1; reaction `r_4598` "cofactor pseudoreaction".
They are *S. cerevisiae* consensus-model values, not *Y. lipolytica* measurements; provisional, as for TPP
and heme a. That cells contain THF, FAD, riboflavin, NAD, NADP and CoA is established/uncited.

**Species identity (partly verified at task start).** FAD, riboflavin and NADP carry KEGG and BiGG ids that
match the yeast-GEM species. THF, NAD and CoA carry only MetaNetX / ChEBI ids (CHEBI:166569, CHEBI:13389,
CHEBI:23354) that were not resolved at task start; the identity check in phase 1 opens them. The model
stores neutral formulas (charge 0) for all six; yeast-GEM stores charged forms. `biomass_C` is a
pseudo-reaction and is not balanced, so this does not block the task, but the NADP name/formula hydrogen
mismatch (C21H28 vs C21H29) is recorded as a finding.

**Medium supply on the benchmark condition (verified from the medium file and the model).** The Task A
medium (`model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd…`) opens pantothenate
R946, nicotinate R1268, riboflavin R1314 and folate R1127 at 1000, plus PABA R971. In the model,
extracellular folate reaches only `m1051[C_cy]`, which has no other reaction (R849 is its only one),
so medium folate is a dead end and THF must be made de novo; 5-formyl-THF uptake R991 exists but is
not in the medium. Riboflavin, pantothenate and nicotinate do enter the cell.

**Census evidence (LP, from the census raw outputs).** Draining the six cofactors (cytosol + mitochondrion
copies, every exchange open) unblocks 51 of the 72 pseudo-gap reactions: UM3 17, UM9 13 (the whole module),
UM14 5, UM33 4, UM7 3, UM40 2, and one each in UM10, UM98, UM108, UM113, UM115, UM131, UM180. Single
drains: THF 19, FAD 16, riboflavin 9, NADP 8, CoA 7, NAD 4.

**Route genes (model GPRs; `*` = on the consensus list; class = 15 % call in the 2026-10-23 sweep).**
Model reaction names; no gene name is verified.

| Route | Reactions (all blocked unless noted) | Genes |
| --- | --- | --- |
| Folate de novo | R377, R274, R14, R277, R278, R273, R749; R271 active | YALI1A10254g* non; YALI1A03069g* non; R273/R749: YALI1A09911g*, YALI1E29179g* (both non; Boolean rule recorded in before.json); R271: YALI1C17947g* ess |
| Riboflavin de novo | R378, R34, R12, R11, R699, R700 | YALI1E02885g*, YALI1B12042g*, YALI1D08339g*, YALI1A12284g* (all non); YALI1F16069g, YALI1C01520g (not listed) |
| Riboflavin → FAD | R697, R698, R306 | YALI1B02806g* non; YALI1D33889g* non |
| Pantothenate → CoA | R616/R1908/R1909, R647, R648, R615, R1092 (R265 mito) | YALI1D33780g* non; YALI1C09380g* non; R648 `YALI1E07601g and YALI1C08173g` (not listed); YALI1D24649g* non; R1092: YALI1F13202g (not listed) |
| NAD(P) | R565 `YALI1E28545g or YALI1E32908g` (not listed); R566/R572 add YALI1E21312g*; R578/R579/R580 NMNAT `YALI1E30400g or YALI1A20611g` (not listed); NAD synthase only as nuclear copies R568/R569 (YALI1A21156g*) | a cytosolic NAD synthase step was not found by name: phase 1 must show the cytosolic NAD supply route |

## 5. Implementation (curation-change skill, steps 2–3)

1. **Data:** `model/curation/biomass/cofactor_biomass_terms.json` (schema as
   `heme_iron_biomass.json`): `curation_id` COFACTOR-BIOMASS-TERMS-20261027, the authorization quote,
   baseline SHA, the six terms with name, formula, charge, compartment, reactions_before, coefficient,
   yeast-GEM species, source URL and file SHA, evidence tier, organism caveat, reviewer; one
   `accepted_before_state` = the exact `biomass_C` stoichiometry and note keys of a1237821.
2. **Step:** `platform/scripts/gem_annotate/cofactor_biomass_terms.py`, `apply_cofactor_biomass_terms()`.
   Fails closed if: the data file is stale (baseline SHA), any species is missing or has a different
   formula/charge/compartment, `biomass_C` composition differs from the accepted before-state, or any
   term is already present with a different coefficient. Returns `applied` / `already_correct`.
   Touches only `biomass_C` stoichiometry and notes.
3. **Flag:** `--cofactor-biomass-terms` in `cli.py` / `main.py`; requires `--lipoylation-module` (the last
   earlier step that edits `biomass_C`); runs after `--acyl-acp-chemistry`, as the final opt-in step, so
   the heme/iron guard's reviewed before-states are untouched.
4. **Test:** `platform/tests/test_cofactor_biomass_terms.py` on an in-memory model: applies on the
   a1237821 state; refuses a wrong before-state, a missing species, a changed formula, a stale SHA;
   second application is `already_correct` and adds nothing; nothing else in the model changes.
5. `python -m pytest` compared with the previous run test by test; then the `code-reviewer` subagent on
   the diff (platform/AGENTS.md).

## 6. Phases, budget and stop conditions

**Phase 1, in memory, before any build (Task A medium + PO1f, Gurobi 1e-9; about 40 LPs).**
Add the six terms to a copy of a1237821 and run:
1. WT with all six terms: optimal, growth > 0.
2. Each term alone (6 solves): optimal; record the growth cost per term.
3. pFBA at the optimum: the supply route of each cofactor (expected: THF de novo through R377 → R271;
   riboflavin through uptake R1314 or de novo; CoA through R946 → R616 → R647 → R648 → R615 → R1092;
   NAD through nicotinate R1268 → NMNAT → NAD synthase; NADP through R565). A route that runs through a
   reaction the census flagged as a compartment artefact (nuclear or envelope copies without transport)
   is a stop.
4. Vitamin uptakes R946, R1268, R1314 closed, all six terms: growth > 0 shows the de novo routes exist;
   if not, record which cofactor needs the vitamin (not a stop, but a limit for other media).
5. Energy-generating-cycle LP (max ATP with every exchange closed): 0, as in the census closure checks.
6. Identity check of THF, NAD and CoA cross-references (ChEBI / MetaNetX pages opened and recorded).
**Stop and report before building** on: a non-optimal or zero-growth WT in 1–2; an artefact route in 3;
an unresolved identity in 6. The user then re-scopes (for example drops a term); nothing is tuned.

**Phase 2, build (3 builds, no solves).**
```bash
python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" --offline --no-solve \
  <the a1237821 argv from its .build.json> --cofactor-biomass-terms \
  --output-model ../model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_20261027.xml
```
Control without the new flag must reproduce `a1237821…` byte for byte; `PYTHONHASHSEED=0` and `=1`
with the flag must give identical bytes. Record the command, output SHA-256 and `.build.json`.

**Phase 3, verification (skill step 4).**
1. Whole-model diff a1237821 vs new: only `biomass_C` (six new coefficients, the new notes) and the
   derived reaction lists of the six species may differ.
2. Solves on the Task A medium + PO1f, Gurobi 1e-9: WT; the six single-term controls repeated on the
   built file; the vitamin-closed control.
3. Knockout sweep: both models, all 1,081 genes, WT, reset per knockout, knocked-reaction flux check;
   classes at 1 / 5 / 10 / 15 %; the user-defined negative class stated; mtDNA genes out of scope.
   Every class change listed with its cause. A new false positive is reported, never fixed.
4. Blocked-set census on both models with the census scripts (open condition and Task A medium; about
   4 × 2 × 2,316 LPs): the reactions that become active, the reactions that become blocked (expected
   none), modules and gap counts, compared with 756 / 701 / 212.
5. Independent read-only audit of TASK.md, the diff, the logs and the deliverables.

**Budget:** about 40 + 10 LPs in phases 1 and 3.2, 2 × 1,082 knockout LPs, 4 FVA scans, 3 builds, pytest.
**Stop and report on:** a guard rejection, a diff outside `biomass_C`, seed-dependent bytes, a
non-optimal WT, a test regression, a phase-1 failure, or any class change whose cause cannot be traced
to one of the six routes.

## 7. Expected effect (inferred unless marked)

- **Blocked reactions (LP upper bound):** at most 51 fewer in the open condition; a cytosol-only term
  frees fewer where the census drain used the mitochondrial copy (NADP 8, NAD 4, CoA 7, THF 19 include
  mitochondrial reactions). Re-census gives the exact number.
- **WT growth** falls by the synthesis cost of the six pools (untested; NAD at 2.65e-3 and riboflavin at
  9.9e-4 are the large terms).
- **Knockouts on the Task A medium (inferred from the GPRs and the medium):**
  - THF: folate uptake is a dead end, so the de novo genes YALI1A10254g and YALI1A03069g should become
    essential (2 FN → TP); the R273 pair depends on its Boolean rule.
  - CoA: YALI1D33780g, YALI1C09380g, YALI1D24649g should become essential (3 FN → TP); the R648 AND pair
    YALI1E07601g / YALI1C08173g and the DPCK gene YALI1F13202g are not on the list and would become new
    false positives if no bypass exists. These are reported as such.
  - Riboflavin / FAD: riboflavin uptake is open, so the de novo genes likely stay non-essential (4 FN
    remain); the kinase YALI1B02806g and the FAD synthetase YALI1D33889g should become essential (2 FN →
    TP).
  - NAD / NADP: R565 has two non-listed genes in an OR, so no class change is expected there;
    YALI1E21312g (mitochondrial copies only) stays FN. The cytosolic NAD route is a phase-1 question.
- **Solver fringe:** FAD at 1e-5 puts its route at about 3e-6 flux; above the 1e-9 tolerance but within
  the range recorded as open problem 1 (STATE.md 2026-10-22). All sweeps on Gurobi 1e-9.

## 8. Acceptance criteria

1. Phase 1 passes with no stop.
2. pytest equals the previous run test by test plus the new test file passing; code-reviewer findings
   fixed or reported.
3. Control build reproduces `a1237821…`; seeds 0 and 1 identical.
4. Whole-model diff limited to `biomass_C` and the six derived reaction lists.
5. WT optimal on the Task A medium + PO1f; each cofactor supplied through a route with a mechanism.
6. Knockout class changes all explained by one of the six routes; false positives reported.
7. Blocked-set census: no reaction newly blocked; the freed set is a subset of the 51.
8. Independent audit verdicts recorded; disagreements reported, not argued away.

## 9. Environment

- Interpreter: the main checkout's `.venv/bin/python` (Python 3.13.5, cobra 0.30.0, optlang 1.8.3,
  pytest 9.0.2, gurobipy 13.0.1); Gurobi with the real `HOME` (licence memory).
- Condition: `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) and
  `state/strain_profiles/po1f_sd_leu.json` (`35307853…`).
- Census scripts: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/gap_metabolites_ponce2013_20261026/scripts/`.
- Large outputs go to `$IYALI26_RESEARCH_ROOT/artifacts/tasks/cofactor_biomass_terms_20261027/`.

## 10. Not in scope (findings only)

- The other 21 pseudo-gap reactions (FMN 11 overlapping FAD, PLP 8, dolichyl-P 6, siroheme 3, biotin 2,
  pantetheine-4'-P 4): no yeast-GEM value; a later decision.
- NADH, NADPH and glutathione terms of r_4598.
- Mitochondrial and other compartment copies of the six cofactors.
- Medium folate being a dead end in the model (R849 → m1051 only); the closed DHFR duplicates
  R1950 / R1951 (same gene as the active R271 / R272, 5 H+ instead of 3).
- The NADP name/formula hydrogen mismatch; the unresolved ChEBI ids of THF, NAD and CoA.
- Why the riboflavin de novo genes are consensus-essential in screens run on YNB, which contains
  riboflavin: an experimental-list question, not a model change.
- Other media, the lipid-unlump line, a benchmark score as a goal.
