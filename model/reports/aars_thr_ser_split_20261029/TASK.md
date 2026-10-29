# Cytosolic threonyl- and seryl-tRNA synthetase rules split from their mitochondrial paralogs (2026-10-29)

## 0. Authorization

**User, 2026-10-29**, answering "Which change should go through curation-change first?" after the OR-rule review
([../or_rule_review_r171_aars_20261029/REPORT.md](../or_rule_review_r171_aars_20261029/REPORT.md)):

> "Thr/Ser split first"

The option the user picked was described in that review (§8, candidate 1) as:
`R759` ← `YALI1F14882g`, `R760` ← `YALI1A09121g`, `R723` ← `YALI1F03908g`, `YALI1B07957g` kept in the model with no
reaction; `R673` (prolyl) unchanged until a presequence prediction is run.

**Gate (user's global instruction):** the plan in §5 is presented for approval before any data file, code, test or
build is written. Status: **approved** — user, 2026-10-29: "Approve as written" (reply button whose text was:
"Approved: implement the Thr/Ser split plan as written in model/reports/aars_thr_ser_split_20261029/TASK.md
section 5, on baseline 965ab947, flag --aars-thr-ser-split. No commit or push."). No commit or push in any case.

## 1. Starting state

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-prediction-e76be3` @ `ff5192d` (= `origin/main`); uncommitted: `model/STATE.md` (OR-review entry) and the two report folders of 2026-10-29 |
| Baseline model | `candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_r2171_1to1_20261028.xml`, SHA-256 `965ab947f2aadc2c69443e7a794358c502acfdf0892e78582278ff07a37429d7` (PR #31) |
| Why this baseline | The OR review ran on it. Its sibling `5f54eb22…` (UM7, PR #30) is also on main; both are `9d655c2b` + one flag and no build combines them. Building on `965ab947` keeps the diff to this change alone. The new flag does not depend on R2171 or UM7, so it can later be added to either line or to a combined build. |
| Baseline build argv | recorded in its `.build.json`: `--offline --no-solve --coq9-curation qcycle --energy-candidate E5 … --dolichol-route --r2171-pi-plc-1to1` (28 opt-in flags) |
| Reference model | none designated |
| Before-state of the targets | [before.json](before.json): rules, bounds, stoichiometry, species, annotations and notes of `R759`, `R760`, `R723` (and `R673`, untouched, as a guard) and the four genes |

## 2. Targets

| Reaction | Compartment | Before (= iYli21) | After |
| --- | --- | --- | --- |
| `R759` threonyl-tRNA synthetase, EC 6.1.1.3 | cytosol | `YALI1A09121g or YALI1F14882g` | `YALI1F14882g` (THS1-like; UniProt Q6C234 "Threonine–tRNA ligase, cytoplasmic") |
| `R760` threonyl-tRNA synthetase, EC 6.1.1.3 | mitochondrion | `YALI1A09121g or YALI1F14882g` | `YALI1A09121g` (MST1-like; UniProt Q6CHF9, location Mitochondrion matrix by automatic annotation) |
| `R723` seryl-tRNA synthetase, EC 6.1.1.11 | cytosol | `YALI1B07957g or YALI1F03908g` | `YALI1F03908g` (SES1-like; UniProt Q6C350 "Serine–tRNA ligase, cytoplasmic", UniRef50 cluster of SES1) |

`YALI1B07957g` (DIA4-like, the mitochondrial SerRS ortholog) loses its only reaction and stays in the model, unlinked,
as `YALI1F07747g` did on 2026-10-22. No mitochondrial seryl-tRNA reaction exists in the model, and adding one is out
of scope. Nothing else changes: no stoichiometry, bound, species, gene annotation or other rule.

## 3. Evidence (from the OR review, §4 and §7.2)

- **Verified:** UniProt names, lengths and locations; UniRef50 memberships; the screen calls (all four genes
  consensus-essential); the rules' iYli21 origin; reaction-level knockouts (`R759` or `R723` closed → growth 0).
- **Inferred:** that each pair is one cytosolic and one mitochondrial enzyme, from homology, length, N-terminal
  composition and the UniProt annotations (two of the locations are automatic ARBA annotations). No
  *Y. lipolytica* experimental localisation exists in the sources searched. AGENTS §5: an OR member needs evidence that
  it catalyses this step alone **in this compartment**; the mitochondrial-type paralogs have none for the cytosol.
- **Not evidence:** the recall gain; Yeast-GEM's rules (its `r_0995` has the same `SES1 or DIA4` defect).

## 4. Expected effect and acceptance criteria

Expected, from the in-memory scoring of the review (verified there on the same baseline):

- `YALI1F14882g` and `YALI1F03908g` FN → TP at 1, 5, 10 and 15 %; no other class change. 15 %: 112 / 203 / 42 / 692 →
  114 / 201 / 42 / 692; 16 unresolved unchanged.
- `YALI1A09121g` and `YALI1B07957g` stay FN (R760 carries no flux at the optimum; no mitochondrial Ser step).
- WT growth unchanged (0.3116608216); a GPR edit cannot change it.

Acceptance:

1. Control build without the new flag reproduces `965ab947…` byte for byte.
2. Whole-model diff of the new candidate against `965ab947…`: only the three rules, the added notes on those three
   reactions, and the reaction links of the four genes (`YALI1B07957g` loses `R723`; `YALI1A09121g` loses `R759`;
   `YALI1F14882g` loses `R760`). Anything else: stop and report.
3. The guarded step fails closed on every stale input in the test, is idempotent, and leaves nothing behind in a
   cobra context.
4. Full pytest: the same pass / fail / error set as before the change, plus the new tests passing.
5. Knockout sweep on both models: the class changes are exactly the two above, at every cutoff; no failing point.

## 5. Plan (for approval; nothing below is written yet)

1. **Data:** `model/curation/gpr/aars_thr_ser_split.json` (schema 1, assignment
   `AARS-THR-SER-SPLIT-20261029`): the authorization quote, and per reaction the before / after rule, name, bounds,
   stoichiometry, species, annotation, before-notes and the notes to add; the unlinked gene with its annotation and
   reactions-before; evidence rows with UniProt URLs, tier and reviewer.
2. **Step:** `platform/scripts/gem_annotate/aars_thr_ser_split.py`, `apply_aars_thr_ser_split(model)`, modelled on
   `r2193_iron_gpr.py`: every precondition of all three reactions and four genes is checked before any edit;
   idempotent with the notes as the record; notes restored on cobra context exit. No gene is added, so there is no
   hash-seed ordering issue.
3. **Wiring:** CLI flag `--aars-thr-ser-split`; called in `main.py` after the R2171 step. No dependency on another
   flag; the guards pin the before-state instead.
4. **Test:** `platform/tests/test_aars_thr_ser_split.py` on the baseline in memory: applied state, every stale-input
   case fails closed with no edits, idempotence, context rollback, SBML roundtrip, Boolean knockout logic per gene
   (no LP), flag parsing.
5. **pytest** before and after, compared test by test.
6. **Build** with the baseline argv plus `--aars-thr-ser-split` to
   `candidates/E5_…_dolichol_route_v2_r2171_1to1_aars_thr_ser_split_20261029.xml`, and a control without the flag to
   a scratch path.
7. **Verify:** whole-model diff; WT; knockout sweep of both models (Task A medium + PO1f, Gurobi 1e-9, reset per
   solve, knocked-flux check).
8. **Audit** by a fresh read-only subagent; REPORT.md; `model/STATE.md` entry; `candidates/README.md` row. Stop
   before any commit.

**Compute budget:** at most 2,300 LPs (2 WT + 2 × 1,081 knockouts + 30 target controls), 30 minutes. Stop
conditions: a failing control build, any unexpected diff, or any sweep point failing the checks; report, never
re-tune.

**Large outputs:** the two sweep TSVs go to `$IYALI26_RESEARCH_ROOT/artifacts/tasks/aars_thr_ser_split_20261029/`.
