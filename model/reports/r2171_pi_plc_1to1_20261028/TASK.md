# TASK: R2171 written 1:1 (Stage 2 of the hidden-source check) (2026-10-28)

**Authorization.** The user, 2026-10-28, after reading this draft (shared notes folder,
`r2171_stage2_TASK_draft_20261028.md`): "Approve Option A", meaning the sentence in section 9 (Option A, the
conventional output name). No commit or push. Sections 0–9 below are the draft unchanged apart from this paragraph
and the status line.

**Status: completed 2026-10-28 on branch `claude/lipid-source-leak-check` (base `origin/main` `0093c89`,
unchanged); candidate `965ab947`, [REPORT.md](REPORT.md), [AUDIT.md](AUDIT.md). Not committed.** Correction found by
the audit: section 5 item 5 says the dolichol `ko_compare.tsv` has a "`9d655c2b` side"; its sides are `ce9e1112` and
`a801aad0` (recorded there as the same optimization problem as `9d655c2b`). The baseline sweep was re-run instead
(REPORT section 6).
Baseline state: [before.json](before.json). Evidence below was read on 2026-10-28 from the model and from database
pages opened that day; no LP was run for the draft.

## 0. Why

Stage 1 ([REPORT.md](model/reports/lipid_source_leak_check_20261028/REPORT.md), branch
`claude/lipid-source-leak-check`, uncommitted) found that R2171 `PI + 100 H2O -> DAG + 100 Ins1P + 100 H+` is an
inositol source (observed in runs):
- INO1 knockout without inositol uptake grows at WT on four media (R2171 2.4e-7, IMPase 2.4e-5).
- Inositol demand 511–520 and Ins1P demand 1000 (the bound); closing R2171 gives 0.
- With all phosphate uptake closed, Ins1P still reaches 1000: R2171 also creates organic phosphate and carbon.
- Opening the DAG kinase R1720 makes the loop independent of R1726, so R2171 should be fixed first.

## 1. Starting state (read 2026-10-28)

- **Base:** `9d655c2b` (`9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c`), `origin/main`
  `0093c89`; built with the argv in its `.build.json` ending `--r1726-cpt-open --dolichol-route`. Not a reference.
- **Branch:** continue on `claude/lipid-source-leak-check` so Stage 1 and Stage 2 land together (Stage 1 files are
  uncommitted there). If the UM7 session or anything else moves `main` first, the base is re-pinned and stated.
- **R2171 now:** bounds (0, 1000), GPR `YALI1D13510g` (= YALI0D10813g), KEGG R03332, EC 4.6.1.13 / 3.1.4.11 /
  3.1.4.3, MNXR148100, no notes, group `group40`; identical in `model/source/iyali26.xml` and `iyli21.xml`.
  Charges: PI `m1393` 0, H2O 0, H+ +1, Ins1P `m813` −2, DAG `m1994` 0 → **charge −100**. Mass not checkable
  (PI `C11H17O13P*2`, DAG no formula).

## 2. Evidence (opened 2026-10-28 unless marked)

| # | Claim | Source | Tier |
| --- | --- | --- | --- |
| E1 | The model's own KEGG link R03332 is `PI + H2O <=> Ins1P + DAG`, 1:1:1:1. | [KEGG R03332](https://www.kegg.jp/entry/R03332) | verified |
| E2 | ModelSEED has the same 1:1 hydrolase (rxn06362, rxn39389, status OK); none of 8 "inositolphosphohydrolase" reactions has 100×. | [ModelSEED rxn06362](https://modelseed.org/solr/reactions/select?wt=json&q=id:rxn06362) | verified |
| E3 | EC 4.6.1.13 is the cyclic-phosphate-forming lyase (R03364: PI → Ins 1,2-cyclic P + DAG); KEGG's comment: "This enzyme is bacterial", animal activity is EC 3.1.4.11. | [KEGG EC 4.6.1.13](https://www.kegg.jp/entry/4.6.1.13) | verified |
| E4 | YALI0D10813g is the only *Y. lipolytica* member of KO K01771; the KO also lists *C. albicans* FGR22 / PLC2 and *Listeria* plcA. | [KEGG K01771](https://www.kegg.jp/entry/K01771), [yli:2910441](https://www.kegg.jp/entry/yli:2910441) | verified (orthology) |
| E5 | UniProt Q6C9J1 (YALI0D10813p): unreviewed, annotation score 1, 305 aa, one "phosphatidylinositol-specific phospholipase C X" domain (SMART SM00148), InterPro IPR000909 / IPR017946 / IPR051057; no function text. | [UniProt Q6C9J1](https://rest.uniprot.org/uniprotkb/Q6C9J1) | verified |
| E6 | The model has no inositol 1,2-cyclic phosphate species (KEGG C04299) (name and annotation scan). | model | verified (static) |
| E7 | The PIP2 phospholipase C R1740 (`YALI1C05570g`) is a separate reaction and is not touched. | model | verified (static) |

**Inferred:** the gene's real activity and product are unknown (X-domain only, unreviewed). Nothing supports the
100× coefficients; they look like R2145's lumped 0.01 scheme multiplied by 100.

## 3. Options

- **A (recommended). Write R2171 1:1:** `PI + H2O -> DAG + Ins1P + 2 H+`, as KEGG R03332 (the model's own link) and
  ModelSEED rxn06362. Two H+ balance charge under the existing charges (0 = 0 − 2 + 2); no charge edits (they
  ripple). GPR, bounds, name and annotations unchanged; notes record the old equation and the weak gene evidence.
  Removes the leak without new species and without removing a reaction.
- **B. Cyclic-phosphate lyase** (`PI -> DAG + Ins 1,2-cyclic P`, R03364), matching EC 4.6.1.13. Needs a new species
  and a route for it (a cyclic phosphodiesterase), none of which has *Y. lipolytica* or yeast support here. Not
  recommended without evidence.
- **C. Close R2171.** Removes the leak, but closing is a removal without evidence that the activity is absent, and
  the gene has a PI-PLC X domain. Not recommended.
- **Not in this stage:** AlphaFold of Q6C9J1 (pre-authorized for identity) could not show catalysis or product, so
  it would not decide A vs B.

## 4. Implementation (Option A; curation-change skill)

- **Data:** `model/curation/chemistry/r2171_pi_plc_1to1.json`: before and after equations; species guards (formula,
  charge, compartment, key annotations of `m1393`, `m32`, `m10`, `m813`, `m1994`); name, GPR, bounds, KEGG link;
  notes (`r2171_previous_equation`, evidence, limits, sources, authorization); evidence rows E1–E7.
- **Step:** `platform/scripts/gem_annotate/r2171_pi_plc.py`, `apply_r2171_pi_plc_1to1`, built like the R2135 step:
  scope and after-equation pinned in code; fails closed on a missing reaction or species, any species / name /
  GPR / bounds change, an unexpected equation, conflicting notes, the notes on the old equation, or the new equation
  without the notes; returns `applied` / `already_correct`; touches only R2171's stoichiometry and notes.
- **Flag:** `--r2171-pi-plc-1to1`, requires `--dolichol-route` and runs after it (`cli.py`, `main.py`).
- **Test:** `platform/tests/test_r2171_pi_plc.py` on the `9d655c2b` fixture, each guard matched by its message.
  The fixture is read through a `Path` object, because cobra refuses string paths of 260+ characters on macOS and
  this candidate's path is 277.
- **Output name:** by convention `..._dolichol_route_v2_r2171_1to1_<date>.xml` (about 289 characters as an absolute
  path). All loads in this task would then use byte-identical short-path copies with recorded SHA-256, as in Stage 1.
  A shorter name is possible but changes the convention; your call (section 8).

## 5. Verification and acceptance

1. **Builds:** the control without the flag reproduces `9d655c2b` byte for byte; `PYTHONHASHSEED` 0 and 1 with the
   flag give the same bytes, the seed written as line 1 of each log.
2. **Whole-model diff vs `9d655c2b`:** R2171 stoichiometry and notes only. R2171 charge balance 0.
3. **Built-file checks** (Task A + PO1f, Gurobi 1e-9): WT growth (expected unchanged; R2171 was 0 at the Stage 1 WT
   optimum); R2171 FVA; no-uptake max ATP and biomass 0 (xMAINTENANCE relaxed).
4. **Acceptance test = Stage 1's check on the candidate** (C1–C4; C5 dropped as invalid): inositol-removed growth 0
   and inositol / Ins1P demand 0 on all four; acyl and phosphate results still 0; the oleate positive control still
   grows; the in-memory R1720 variant on C1 also shows no inositol leak.
5. **Knockout sweep:** 1,082 knockouts on the candidate, Gurobi 1e-9, reset per knockout, knocked-flux check; the
   baseline is the `9d655c2b` side of `model/reports/dolichol_route_20261028/ko_compare.tsv` (same medium `6eca70fd`
   and profile `35307853`), reused only if the sweep code and settings match, otherwise re-run. Scores at 1, 5, 10
   and 15 % with the user-defined negative class, mtDNA out of scope. Expected: no class change on Task A, because
   the medium supplies inositol (inferred).
6. **Re-census** (GLPK, open and Task A) against the `9d655c2b` census of the dolichol task.
7. **Code review** (fresh agent, `.claude/agents/code-reviewer.md`), **independent read-only audit**, REPORT.md,
   `model/STATE.md` entry, and an **Escher before / after** of the Stage 1 maps F2 / F3 on the candidate.

## 6. Budget and stop conditions

- Builds: 3 (control, seed 0, seed 1), plus up to 3 more if review fixes touch the step.
- LPs: acceptance check ≤ 150; built-file checks ≤ 30; knockout sweep 1,082 (+ ≤ 50 fresh re-solves; + 1,082 only if
  the baseline cannot be reused); census 2 GLPK runs.
- **Stop** on: control not byte-identical or seeds differing; diff beyond R2171; the acceptance test still showing an
  inositol or Ins1P source; any new acyl or phosphate source; a WT change above 1e-6 or any knockout class change not
  traced to R2171; `main` moving the base.

## 7. Not in scope

The IPC / MIPC synthases (0.01 PI), R2165, R2187; opening R1720; a valid open condition; changing the context loader
to pass `Path` objects (the 260-character limit); the lipid-unlump species (colleague-owned); commit or push.

## 8. Decisions for you

1. Option A, B or C (A recommended).
2. Output file name: the convention (with short-path copies) or a shorter name.
3. Committing Stage 1 is separate; it can go into the same PR as Stage 2 when you say so.

## 9. Authorization sentence (Option A, convention name)

> Approved: Stage 2 Option A: write R2171 as PI + H2O -> DAG + Ins1P + 2 H+ via a guarded --r2171-pi-plc-1to1 step
> after --dolichol-route on 9d655c2b, conventional output name, with test, builds, diff, Stage 1 acceptance check,
> knockout sweep, re-census, review, audit, REPORT, STATE and Escher. No commit or push.
