# TASK (DRAFT, awaiting approval): delete the legacy biomass pseudo-reactions model-wide (2026-11-12)

Status: **decisions approved, waiting for the port** (2026-11-12). Nothing is built, changed or committed. The user
approved L1–L8 as recommended (section 0). Stage 1 computation, the code, the build and any commit each still need
their own approval. The investigation is static only (cobra reads, no LP, no FVA). Its scripts, outputs and
SHA256SUMS are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/legacy_biomass_cleanup_20261112/static_review/`
(below: `static_review/`).

## 0. Authorization so far

- User, 2026-11-12 (strict-sn port, D8): close newBiom, xBIOMASS, xLIPID, BIOMASS_C_core and BIOMASS_C_full in the
  port, and handle deletion of the legacy biomass family as a separate model-wide cleanup (relayed by the starting
  session; recorded in the port's draft TASK.md, D8).
- This round: investigate and plan only. No build, code edit, curated-data edit, LP or commit.
- User, 2026-11-12, on this draft: **"Approve as drafted"**, in answer to "Approve the L1–L8 recommendations in the
  cleanup TASK.md?". Read as every recommendation in section 11:
  - L1 (a): delete as in section 5;
  - L2 (a): R1371 goes too;
  - L3 (a): m1394 stays as a flagged orphan;
  - L4 (a): R1387/R1710 are deleted in the last step;
  - L5 (a): BIOMASS_C_core/full stay closed, and Andrew is asked about deleting them;
  - L6 (a): the R1387 note is kept and a cleanup note is added;
  - L7: Option A;
  - L8 (a): after the port, on its candidate.

  Not covered by this approval: the Stage 1 computation (section 8), writing the code, the build, sending the
  message to Andrew, and any commit or PR.
- **What blocks the next step:** L8 (a) puts Stage 1 on the port candidate. The port is not built, and its D0 is
  waiting on Andrew. If the port stays held, L8 (b) (run on 684f3752) needs the user's choice.

## 1. Starting state (checked 2026-11-12)

- Worktree `great-chebyshev-bc86db`, branch `claude/great-chebyshev-bc86db`, HEAD `efd7f780` (merge of PR #51), clean.
- Models read (static only, through short-path symlinks with the same SHA):
  - `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2`, the brief's model
    (`…_mcnaugh_rxns_v3_20261103.xml`, PR #47);
  - `684f3752764fb8f25219008966953863141724f34a451dffb77c235aa3cb2f62`, the newest candidate on `main`
    (`…_mcnaugh_rxns_v3_r2121_close_r1172_close_20261112.xml`, PR #51);
  - `model/source/iyali26.xml` (`5c8c199e…`), for the source bounds.
- **The family is identical in e584151f and 684f3752**: bounds, stoichiometry size, SBO, groups and notes keys
  (`static_review/outputs/family_static.txt`). Both have 2,346 reactions, 1,883 species, 1,081 genes.
- Reference model: none designated (`model/candidates/README.md`).
- Pending port: `model/reports/strict_sn_lipid_port_20261110/TASK.md`, uncommitted, in worktree
  `confusion-matrices-iyali21-7e63ca` (branch `claude/strict-sn-lipid-port-2342e4`). D0 (Andrew) is still open.

## 2. Question

Can the legacy biomass pseudo-reactions be removed from the built model, which ones, and what has to change in the
builder, tools, tests and guard data so that the build, the earlier guards and the exporters keep working?

They are modelling constructs, not biochemistry, so the AGENTS.md deletion rule (evidence that the reaction does not
exist in the organism) does not apply to them as it does to an enzyme step. The reason to remove them is model
clarity, not a score: a memote gain is not a reason and is not claimed.

## 3. Static facts (verified this session unless marked)

### 3.1 The family

| ID | Name | Source bounds | Candidate bounds | Products | Can it carry flux? (static) |
| --- | --- | --- | --- | --- | --- |
| newBiom | newBiom | 0..1000 | 0..1000 | biomass_ m1401 | every substrate has an open producer |
| xBIOMASS | xBIOMASS (iYali objective) | 0..1000 | 0..1000 | m1401; eats protein_ m1726 + lipids_ m1727 | yes, with xAMINOACID and xLIPID |
| xAMINOACID | Amino acid composition | 0..1000 | 0..1000 | protein_ m1726 (only consumer: xBIOMASS) | yes |
| xLIPID | Lipids pool | 0..1000 | 0..1000 | lipids_ m1727 (only consumer: xBIOMASS) | yes |
| R1372 | yeast 5 biomass pseudoreaction | **0..0** | 0..1000 | m1401; eats lipid_ m1400 | **no: m1400's only producer is R1371, closed** |
| R1371 | lipid pseudoreaction | 0..0 | 0..0 | lipid_ m1400 (consumers R1372, R1387, R1710 only) | no (closed; also eats ergosterol ester m1394, which nothing makes) |
| R1387 | yeast 6 biomass pseudoreaction | 0..0 | 0..0 | m1401 + uncharged tRNAs; eats m1400 | no |
| R1710 | yeast 8 biomass pseudoreaction | 0..0 | 0..0 | same as R1387, plus heme a | no |
| BIOMASS_C_core | Andrew's core biomass | absent | 0..1000, objective 0 | ADP, Pi, H+, tRNAs; **no m1401** | every substrate has an open producer |
| BIOMASS_C_full | Andrew's full biomass | absent | 0..1000, objective 0 | same; **no m1401** | every substrate has an open producer |

Not in scope and kept: `biomass_C` (sole objective, coefficient 1), `R1373` "growth" (the only consumer of m1401),
`xMAINTENANCE` (NGAM, 7.8625..1000), the `xPOOL_*` pools, `TRNA_BIOMASS_*` and `R_DM_GPI`.

Sources: `family_static.txt`, `substrates.txt`, `pseudo_684f.txt`.

- **None of the ten has a GPR**, so no gene loses a reaction and no knockout call can change through the GPR.
- **R1372 is open only because our builder opens it.** `fix_biomass_reaction` (`platform/scripts/gem_annotate/biomass.py`,
  run as "Priority 3" from `main.py:326`) sets it to 0..1000. Its lipid input m1400 has no open producer, so it cannot
  carry flux. The brief's "R1372 (open)" is nominal.
- **R1371 is a member of the family the brief did not list.** It only feeds m1400, which only R1372, R1387 and R1710
  consume.
- **BIOMASS_C_core and BIOMASS_C_full are pure sinks.** They do not make m1401, so their flux leaves the model
  instead of reaching R1373.
- **Species that are left without any reaction** (`dead_ends_684f.txt`, `pseudo_684f.txt`):
  - deleting the nine of the brief: protein_ m1726 and lipids_ m1727; lipid_ m1400 keeps only the closed R1371;
  - adding R1371: lipid_ m1400 and ergosterol ester m1394 too. m1394 is used by R1371 only.
- **No other species loses its last open producer or consumer** when the five port drains or all nine go
  (first-order incidence; flux effects through longer chains are not checked, see section 8).
- **Objective:** `biomass_C` alone in both candidates. None of the family is in the objective.

### 3.2 Andrew's own design

His strict-sn generator `platform/tools/lipid/lp_sn12_candidate.py:421` **closes** newBiom, xBIOMASS and xLIPID by ID;
it does not delete them. BIOMASS_C_core/full come from his mcnaugh_pass handoff and were merged "Both verbatim from
Andrew" (user, 2026-11-03; `model/STATE.md`, PR #47 entry).

## 4. Where the IDs are used (static read of the code, 2026-11-12)

The key point: **almost every dependency is upstream in the build.** Each step's guard runs while the reactions still
exist. A deletion step that runs **last** leaves them untouched. Deleting early (in the source or before those steps)
would break each of them.

| Place | What it does with the IDs | Effect of a last-step deletion | Effect of an early deletion |
| --- | --- | --- | --- |
| `patches.py` `couple_trna_charging_to_biomass` (template `R1387`, excludes `R1387`/`R1710` from the charging-reaction search) | builds the 20 tRNA pairs. Called from `split_trna_charging_from_biomass` (`main.py:696`, canonical or `--trna-biomass-mode split`) and from the `couple_trna_biomass` patch operation (not used: the patch table holds only EG-GPR-001) | none: runs before | fails (`KeyError` → `ValueError`) unless the template moves into curated data |
| `metadata_reaction_selection.json` | carries R1372's stoichiometry and the split-tRNA `biomass_C`; applied at `main.py:721` | none | R1372 row becomes "reaction missing" |
| `biomass.py` `fix_biomass_reaction` | opens and diagnoses R1372 (Priority 3) | none; it already tolerates a missing R1372 (logs a warning) | none |
| `main.py:483–518`, `add_sbo_terms.py:29` | regexes that assign SBO terms | none | none |
| Guard data: `trehalose_biomass.json`, `cofactor_biomass_terms.json`, `heme_iron_biomass.json`, `tpp_requirement.json` (`reactions_before/after`, `notes_before` with `experimental_trna_biomass_template: R1387`) | pin the state at their own step | none: their steps run before | each guard fails |
| `biomass_C` note `experimental_trna_biomass_template = "R1387"` | provenance text | the note points at a removed ID (decision L6) | same |
| **`platform/tools/export_bigg_ids.py:82–85`** | **raises "Crosswalk entities not in the model"** if any row of Andrew's SHA-pinned crosswalk is missing | **breaks**: the crosswalk lists R1371, R1372, R1387, R1710, newBiom, xBIOMASS, xLIPID, xAMINOACID (→ `PROT_BIOMASS`, "identical") and m1394, m1400, m1401, m1726, m1727 | breaks |
| `tools/validate_energy_candidates.py:270` | records `flux.get('R1372')` | `None` instead of 0 (R1372 always carried 0, being blocked) | same |
| `tests/test_closed_energy.py:17` | asserts `newBiom` closed by `close_model` | none: it loads its own pinned baseline (`reports/dipeptide_energy_audit_20260924/config.json`) | none |
| `tests/test_mcnaugh_reactions.py:71` | checks the curated JSON lists BIOMASS_C_core/full | none: reads the JSON, and the mcnaugh step still adds them | none |
| `tests/test_cofactor_biomass_terms.py`, `test_heme_iron_biomass.py`, `test_trehalose_biomass.py`, `test_ceramide_pool_biomass.py`, `test_dolichol_route.py` (objective `R1710` / `xBIOMASS` in negative tests) | load their step's pinned baseline file | none | none on the baselines |
| `tests/test_reaction_selection.py:157` | checks R1372's selected stoichiometry on its own fixture | none | none |
| `tools/lipid/lp_sn12_candidate.py:421` (Andrew's) | closes newBiom/xBIOMASS/xLIPID by ID | fails only if run on a cleaned candidate; it runs on his base `bc2aac8f` | same |
| `bigg_crosswalk_annotation_decisions_20261103.csv` | no row for the family (the grep hits are MetaNetX IDs containing "R137…") | none | none |
| `bigg_export_policy_20261103.json` | m1726 "PROTEIN" verdict row | the row's entity is gone; harmless unless the exporter requires it | same |

## 5. Per-reaction recommendation

| ID | Recommendation | Reason |
| --- | --- | --- |
| newBiom | **delete** | open alternative composition into R1373; no notes, no origin record; superseded by `biomass_C` |
| xBIOMASS | **delete** | iYali objective, superseded; open into R1373. The iYali numbers stay in the immutable source (the CL:PC ratio decision of 2026-11-12 cites iYali, not our copy) |
| xAMINOACID, xLIPID | **delete** with xBIOMASS | their only outlets are m1726/m1727 → xBIOMASS |
| m1726 protein_, m1727 lipids_ | **delete** | orphaned |
| R1372 | **delete** | structurally blocked; open only through our Priority-3 step |
| R1371 (added to scope) | **delete** (L2) | feeds only the three yeast 5/6/8 reactions |
| m1400 lipid_ | **delete** with R1371 | orphaned |
| m1394 ergosterol ester | **keep as an orphan species and flag it** (L3) | a real compound class; its only reaction is the pseudo-reaction. The missing steryl-ester route is a gap for another line, not a reason to delete the species |
| R1387, R1710 | **delete, last step only** (L4) | closed since the source; only use is the tRNA template during the build, which runs before. The source file keeps both |
| BIOMASS_C_core, BIOMASS_C_full | **keep closed (as the port does); delete only if Andrew agrees** (L5) | Andrew's content, merged verbatim by the user's choice; useful as an in-memory comparison objective (screen 97/229/29/709 at 15 % on 2026-11-04). Their equations stay in `mcnaugh_reactions_20261103.json`, so deletion can be undone by the data |
| m1401 biomass_, R1373, xMAINTENANCE, biomass_C | **keep** | after the cleanup `biomass_C` is m1401's only producer, so R1373 flux equals growth |

"Keep closed" and "delete" give the same LP for every closed member. The difference is the file: IDs, notes and
group membership.

## 6. Mechanism

### Option A (recommended): one opt-in step that runs last

- Flag `--legacy-biomass-cleanup`, module `platform/scripts/gem_annotate/legacy_biomass_cleanup.py`, called after
  `--r1172-hmgcoa-transport-close` (the current last step, `main.py:967–973`) and after the port when that lands.
- Data: `model/curation/biomass/legacy_biomass_cleanup_20261112.json` with:
  - the authorization text and the curation ID (`LEGACY-BIOMASS-CLEANUP-20261112`);
  - per target, the pinned before-state: stoichiometry, bounds, GPR (empty), SBO, groups and notes keys. Two
    allowed bound states for the port's five (open with the port off, 0..0 with it on), and nothing else;
  - the species to remove (m1726, m1727, m1400) and the species to keep as orphans (m1394, if L3 = keep);
  - the expected after-state: objective `biomass_C` only; m1401 made only by `biomass_C` and used only by R1373.
- Guard, fail closed before any change:
  - every target exists and matches its pinned state; objective is `biomass_C` alone;
  - no target has a gene; no target is in `retained_reactions.json`;
  - the species set that would be orphaned equals the listed one exactly.
- Apply: `model.remove_reactions(targets, remove_orphans=False)` (cobra 0.30.0 also removes them from groups), then
  remove the listed species. Idempotent: an already-clean model returns `already_applied`.
- Record a cleanup note (L6) and the removed equations in the build audit, so the file names what went.

### Option B (not recommended): remove early or in the source chain

It would need a curated tRNA template in place of R1387, and rewritten guards and `notes_before` in four biomass
data files. It would also break the reproducibility of every earlier candidate's argv. No gain over A.

## 7. Changes each choice needs (Option A)

| Item | Delete (A) | Keep closed only |
| --- | --- | --- |
| New step module, data file, flag (`cli.py`), call in `main.py` | yes | the closures already come from the port |
| New test `platform/tests/test_legacy_biomass_cleanup.py`: guard refuses a changed target, an extra objective, a gene on a target, an unexpected orphan; idempotence; SBML round trip; both port states | yes | no |
| `tools/export_bigg_ids.py`: accept crosswalk rows whose entity is listed as removed in the cleanup data, report them "removed by LEGACY-BIOMASS-CLEANUP-20261112", keep refusing any other missing entity; plus a test. Andrew's crosswalk is not changed (SHA-pinned) | yes | no |
| `bigg_export_policy_20261103.json` m1726 row | leave; the exporter change reports it | no |
| `tools/validate_energy_candidates.py` R1372 column | leave (`None` = absent; it was always 0) | no |
| `patches.py` tRNA functions, `biomass.py`, `main.py` regexes, `add_sbo_terms.py` | no change | no change |
| Four guard JSONs (trehalose, cofactor, heme/iron, TPP) | no change: they describe their own step | no change |
| Existing tests | none expected to change; confirmed by a before/after suite run (about 25 min; baseline on a short-path archive copy) | none |
| `model/candidates/README.md`, REPORT.md, `model/STATE.md` | yes | yes (port's) |

## 8. What has to be checked (proposed Stage 1, needs approval as computation)

Inputs: 684f3752 (or the port candidate once built), Task A medium
`model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` + PO1f profile, Gurobi 1e-9, checked solve status.
In memory; no file written to `model/`.

1. WT growth with the family removed vs present. Expected identical: removing a sink cannot raise a maximum, and
   none of the family makes a `biomass_C` precursor (inferred). Any drop is reported, not tuned.
2. FVA on every reaction, before vs after: list reactions that become blocked. The drains could be masking blocked
   reactions through longer chains, for example lipid species whose only flux outlet is xLIPID → xBIOMASS → R1373.
   By AGENTS.md §5, any newly blocked reaction is a gap to report, especially on essential genes.
3. ATP with no carbon uptake stays 0 (EGC check, with the xMAINTENANCE lower bound relaxed and status checked).
4. Certified knockout screen at 1, 5, 10 and 15 %, rules A and B, only if step 2 shows a change on a gene-carrying
   reaction. Otherwise a growth equality is enough, since no GPR changes.

Budget: about 2 LPs plus one FVA per state; 1 screen at most. Stop if WT changes or a reaction with a gene becomes
blocked, and report.

Then the build (one build plus a spare), whole-model diff (only the listed reactions, species and group memberships
removed, plus the cleanup note), independent read-only audit, REPORT.md and a `model/STATE.md` entry, as the
curation-change skill describes.

## 9. Interaction with the strict-sn port

- The port closes newBiom, xBIOMASS, xLIPID, BIOMASS_C_core and BIOMASS_C_full (D8 b). It runs right after
  `--mcnaugh-reactions` and pins the live bounds of its targets at that point. The cleanup runs last, after the port,
  so the port's guard sees the reactions as they are today. The two lines do not block each other.
- The cleanup guard accepts both states of the five: open (port off) and 0..0 (port on).
- Under the port, xLIPID's generic lipid terms have no source (port draft §4.3, inferred). Deleting it there is
  LP-neutral by construction. Without the port, Stage 1 step 2 is what shows whether xLIPID masks anything.
- Order: recommended **after the port**. Then the cleanup's Stage 1 runs on the port candidate, the one that will be
  the chain's base. If the port stays held on D0, the cleanup can run on 684f3752. Only the `cli.py`/`main.py` flag
  lines would need a merge.
- The CL-in-`biomass_C` follow-up (after the port, with Andrew) does not use any family member.

## 10. Does Andrew need to agree?

- **newBiom, xBIOMASS, xLIPID, xAMINOACID, R1371, R1372, R1387, R1710:** no. They come from our source model, are not
  in his ownership paths (`platform/tools/lipid/`, `model/curation/lipid_unlump/`, `model/candidates/lipid_unlump/`),
  and his generator touches only our model's closures. Recommended courtesy: tell him, since his tool closes three of
  them by ID and his crosswalk lists them.
- **BIOMASS_C_core/full: yes, recommended.** AGENTS.md does not list `model/curation/mcnaugh_pass/` under his
  ownership (verified by reading AGENTS.md), so the formal rule does not force it. But they are his reactions; the
  user chose "Both verbatim from Andrew"; and the port's D0 already needs his agreement. Ask in the same message.
  Until he answers, keep them closed (L5).

## 11. Decisions for the user

| # | Decision | Options | Recommendation |
| --- | --- | --- | --- |
| L1 | Go ahead with a cleanup at all | (a) delete as in section 5; (b) keep everything closed and only document it; (c) stop | (a). It is LP-neutral (to be shown in Stage 1); the gain is one biomass definition, no accidental alternative objective, and four fewer orphan pseudo-species. Plain (b) is acceptable: it gives the same LP |
| L2 | R1371 (not in the brief) | (a) delete with R1372; (b) keep closed | (a): it feeds only the deleted reactions |
| L3 | m1394 ergosterol ester | (a) keep as an orphan species, flag it as a gap; (b) delete | (a): a real compound; deleting it needs the biology argument the rule asks for |
| L4 | R1387, R1710 | (a) delete in the last step; (b) keep closed as in-file templates | (a): the template is used only earlier in the build, and the source keeps both |
| L5 | BIOMASS_C_core/full | (a) keep closed (the port), ask Andrew about deletion; (b) delete now; (c) leave open if the port does not land | (a) |
| L6 | `biomass_C` note `experimental_trna_biomass_template = R1387` | (a) keep and add a `legacy_biomass_cleanup` note naming the removed IDs; (b) rewrite it | (a): history stays readable; no earlier guard is affected |
| L7 | Mechanism | Option A (last step) or Option B (early) | A |
| L8 | Order with the port | (a) after the port, on its candidate; (b) now, on 684f3752 | (a); (b) if D0 stays open |

Also needs authorization when the time comes: Stage 1 computation (section 8), the build, and any commit or PR.

## 12. Limits of this draft

- Static only. "LP-neutral" and "no newly blocked reaction" are inferred until Stage 1 runs.
- The dependency map comes from a grep of `platform/` and `model/curation/` for the IDs plus reading each hit. A
  script that builds an ID by string concatenation would be missed.
- The research workspace and notebooks outside `platform/` were not searched for code that reads these IDs.
- Whether `split_trna_charging_from_biomass` runs in the E5 argv (canonical build) or the split comes from the
  metadata selection was not resolved. Both happen before `main.py:721`, so the conclusion holds either way.
