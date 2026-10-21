# R2041 (D-arabinokinase): direction review and restore plan

Read-only. No model, data or code was changed, nothing was built, and no LP was solved. Scope, inputs
and SHAs are in [TASK.md](TASK.md).

Evidence tiers:
- **verified:** checked this session (2026-10-19) and the output is shown here;
- **recorded:** taken from a project record and not re-run;
- **provable:** follows mathematically from verified or recorded facts;
- **inferred:** reasoned, with no single source.

## 1. Summary

- **The loop the patch describes does not exist in a56f3649 (verified + provable).**
  - D-arabinose 5-phosphate (`m1884`) has exactly two reactions: R2041 and the isomerase R2042. There
    is no A5P phosphatase in a56f3649 or in the source model.
  - The only closed cycle through reverse R2041 also runs the xylulokinase R803. Its net
    stoichiometry is exactly zero, so it is an ATP-neutral internal loop, not an energy-generating
    cycle.
  - The closed-input ATP dissipation test already run on a56f3649's baseline `0ac2addf` (R2041
    reversible there too) returned 0 in C_cy and in C_mi (recorded). a56f3649 is `0ac2addf` minus
    R132 (verified), and removing a reaction cannot raise an LP maximum. So the same test on
    a56f3649 is ≤ 0, that is 0 (provable).
- **The case for forward-only rests on thermodynamics, not on a loop.**
  - eQuilibrator gives ΔrG'° = −24.7 ± 4.7 kJ/mol, K'eq = 2.1 × 10⁴ (pH 7.5, pMg 3, I = 0.25 M,
    298.15 K; verified). Because the reaction has 2 substrates and 2 products, ΔrG'm (1 mM) is the
    same.
  - Reverse flux needs [A5P]/[D-arabinose] > K'eq × [ATP]/[ADP]. Even at the 2σ edge
    (−15.3 kJ/mol, K'eq ≈ 480) and ATP/ADP = 1, that is a 480-fold excess of the phosphate ester
    (derived).
  - Enzymology is weak but consistent. ENZYME and Rhea write the kinase direction. Rhea links no
    UniProt sequence. No *Y. lipolytica* enzyme is known, and R2041 has no GPR. The single primary
    reference (Volk 1962) could not be opened.
- **Recommended change (needs approval): direction only.**
  - Set R2041 to [0, 1000] and keep the selection's 2-H⁺ equation. That equation is formula- and
    charge-balanced with the stored species.
  - Do not bring back the patch's 1-H⁺ equation. With the current neutral ATP/ADP it leaves an
    H −1 / charge −1 residual, against `platform/AGENTS.md` (chemistry).
  - Use a post-selection direction step modeled on the existing R1931 and R1159 steps
    (`model/curation/direction/*.json` plus a guarded `apply_*`), not a selection amendment (§6).
- **Expected effect: none on growth (inferred).** All 109 recorded R2041 flux values in 51 report
  files are 0. Those are single optima on several models; FVA has not been run. The LP plan in §6.3
  (10 core, 2 optional) would confirm this.

## 2. R2041 in each state

R2041 is `ATP + D-arabinose ⇌ ADP + D-arabinose 5-phosphate + n H⁺`, in C_cy, EC 2.7.1.54. It has
no GPR in every state.

| State | Bounds | n H⁺ | Residual with stored species | Where |
| --- | --- | --- | --- | --- |
| Source `iyali26.xml` | [−1000, 1000] | 0 | (not checked) | `model/source` |
| Builder patch `fix_d_arabinokinase_direction_and_proton` (2026-08-14, `badc826`) | [0, 1000] | 1 | H −1, charge −1 (stated in the patch comment) | selection "before" |
| Selection 2026-09-10 (`metadata_reaction_fields_preserve_trna_20260910`) | [−1000, 1000] | 2 | none (verified) | a56f3649, c98961e9, 51c6fabd (verified) |

Stored species (verified, a56f3649):
- ATP `m141` is C10H16N5O13P3, charge 0;
- ADP `m143` is C10H15N5O10P2, charge 0;
- D-arabinose `m454` is C5H10O5, charge 0;
- A5P `m1884` is C5H9O8P, charge −2;
- H⁺ is +1.

H balance with 2 H⁺: left 16 + 10 = 26; right 15 + 9 + 2 = 26. Charge: 0 = 0 + (−2) + 2.

The patch's 1 H⁺ is the right count only when ATP and ADP are the charged microspecies (ATP⁴⁻,
ADP³⁻). The patch says so itself: `remaining_gate` is the "ATP/ADP connected-component microspecies
migration". That migration has not happened. The 2-H⁺ equation is therefore the balanced one for
this model, and the proton count must be re-reviewed if the migration ever happens.

In a56f3649 the patch's note survives as `metadata_previous_curated_reaction_correction`. The
reaction also carries the selection's own note: "User-selected model version; not new validation of
chemistry or native gene function" (verified). The selection row is therefore not evidence for
reversibility.

## 3. Static trace in a56f3649

Every reaction that uses D-arabinose (`m454` C_cy, `m1193` C_ex), A5P (`m1884`), D-arabinitol
(`m1823`), D-xylulose (`m1015`), D-ribulose 5-P (`m96`) and D-xylulose 5-P (`m938`) was listed. Only
the ones that matter for R2041 are shown.

| Reaction | Equation (stored IDs) | Bounds | GPR | EC |
| --- | --- | --- | --- | --- |
| R2041 | ATP + D-arabinose ⇌ ADP + A5P + 2 H⁺ | [−1000, 1000] | — | 2.7.1.54 |
| R2042 | Ru5P ⇌ A5P + 2 H⁺ | [−1000, 1000] | — | 5.3.1.13 |
| R1063 / R1062 | D-arabinose transport / exchange (secretion only) | [−1000, 1000] / [0, 1000] | — | — |
| R2190 | NAD + D-arabinitol ⇌ NADH + D-arabinose | [−1000, 1000] | YALI1F03452g | 1.1.1.250 (§7) |
| R254, R255 | D-arabinose → D-arabinono-1,4-lactone (NAD / NADP) | [0, 1000] | YALI1F38466g / 9 genes | 1.1.1.116 / — |
| R253 | lactone + O₂ ⇌ dehydro-lactone `m453` + H₂O₂ | [−1000, 1000] | YALI1A22228g | 1.1.3.37 |
| R1964 | D-xylulose + NADH ⇌ NAD + D-arabinitol | [−1000, 1000] | — | — |
| R1980 / R1981 | D-arabitol transport / exchange (secretion only) | [−1000, 1000] / [0, 1000] | — | — |
| R803 | D-xylulose + ATP → ADP + Xu5P | [0, 1000] | YALI1F14583g | 2.7.1.17 |
| R714 | Ru5P ⇌ Xu5P | [−1000, 1000] | YALI1C16621g | 5.1.3.1 |

Findings (verified unless marked):

- **No phosphatase.** `m1884` takes part in exactly R2041 and R2042, in a56f3649 and in the source
  model. Neither D-arabinose, Ru5P nor Xu5P has a phosphatase among its reactions.
- **Lactone branch is dead.** `m453` has one reaction (R253), so R253 cannot carry steady-state flux.
  R254 and R255 (lactone only via R253) are therefore blocked as well. This is inferred from
  connectivity; FVA was not run.
- **The one closed cycle through reverse R2041** runs R2041 reverse, R2190 reverse, R1964 reverse,
  R803 forward, R714 reverse and R2042 forward. Each direction is allowed by the stored bounds:

  ADP + A5P + 2 H⁺ → ATP + Ara; Ara + NADH → arabinitol + NAD; arabinitol + NAD → xylulose + NADH;
  xylulose + ATP → ADP + Xu5P; Xu5P → Ru5P; Ru5P → A5P + 2 H⁺.

  The summed stoichiometry is exactly zero. The ATP that reverse R2041 makes is spent by R803, so the
  cycle makes no ATP. It is an internal (type III) loop, which matters for loopless or flux-sampling
  analyses, but it is not an energy-generating cycle.
- **With R2041 forward-only** that cycle cannot run in either direction. Running it forward would need
  Xu5P back to D-xylulose, but R803 is forward-only and no Xu5P phosphatase exists.
- **What reverse R2041 does allow in the open network (inferred):** Ru5P → A5P → ATP + D-arabinose,
  with D-arabinose (or D-arabinitol through R2190) secreted. That is "ATP from a sugar phosphate
  ester", which thermodynamics forbids (§5). From glucose through the oxidative PPP it is ATP-neutral
  per glucose, because hexokinase pays 1 ATP. It wastes carbon. An LP would use it only when that
  helps the objective, for example in a knockout that leaves PPP carbon with no other exit. That has
  not been tested.

The patch's stated rationale, "the D-arabinose-phosphate isomerase / phosphatase loop", is therefore
**not reproduced** on this network. It may have been general reasoning or about another model state.
The commit (`badc826`, 2026-08-14) records no LP witness. Its focused test
`tests/test_d_arabinokinase_correction.py` existed on the old line (`a2a998a`). It is absent from
`285eaee`, dropped at the merge `65b81b8` (2026-10-07) (verified with `git cat-file` / `git log -m`).

## 4. Existing LP evidence (no new solves)

**Closed-input ATP dissipation on `0ac2addf` (recorded).**
- Source: `coq9_respiratory_package_20261008/controls.json`, runner `run_controls.py`, protocol
  `tools/diagnose_closed_energy.close_model`.
- Protocol:
  - every single-sided reaction (exchange, sink, demand) and every biomass reaction is set to
    [0, 0];
  - the forced maintenance flux is released;
  - the solver maximizes xMAINTENANCE (C_cy) or a balanced mitochondrial ATP drain (C_mi);
  - context: PO1f overlay, Gurobi with 1 thread, Presolve 0, tolerances 1e-7.
- Result: `ATP_closed_C_cy` and `ATP_closed_C_mi` both optimal with objective 0.0.
- R2041 in `0ac2addf` is [−1000, 1000] with 2 H⁺ (verified).

**Carry-over to a56f3649 (verified + provable).** A static whole-reaction comparison found:
- only R132 is missing in a56f3649;
- no other reaction differs in bounds, stoichiometry or GPR;
- no metabolite differs in formula or charge.

This agrees with the recorded `r132_mito_adh_20261012/whole_model_diff.json`. Deleting a reaction
fixes one flux variable to 0, so the feasible set can only shrink. With the same closure and context,
a56f3649's closed-input ATP maximum is therefore ≤ 0, hence 0.

The overlay used in that test does two things (verified, profile `po1f_sd_leu_accrispr_v1`):
- it sets R612 (URA3 OMP decarboxylase) to [0, 0];
- it swaps R45's GPR to the plasmid LEU2.

So the recorded zero is for the overlaid model. A cycle that needs R612 in the plain model is not
excluded by that record. It is unlikely (inferred), because R612 is an irreversible decarboxylation
whose CO₂ would have to be re-fixed with no exchange. Optional LPs P1 and P2 (§6.3) close this gap.

**Recorded R2041 fluxes (verified scan).**
- Search: every `"R2041": <number>` in `model/reports/**/*.json|csv|tsv`.
- Result: 109 values in 51 files, all 0.
- These are single optima on several earlier models (E5 validation, CoQ candidates, screens,
  diagnostics), not FVA, and not a56f3649 WT. They show the reverse direction was not used in those
  optima. They do not show it is unusable.

## 5. Direction evidence (AGENTS.md §5)

| Line of evidence | Finding | Tier |
| --- | --- | --- |
| Enzyme definition | ENZYME [EC 2.7.1.54](https://enzyme.expasy.org/EC/2.7.1.54): "D-arabinose + ATP = D-arabinose 5-phosphate + ADP + H(+)" | verified |
| Rhea | [RHEA:24588](https://www.rhea-db.org/rhea/24588): the same equation, with 0 UniProt enzymes linked. KEGG [R01573](https://www.kegg.jp/entry/R01573) links the bidirectional RHEA:24591 | verified |
| Characterized enzyme | KEGG [EC 2.7.1.54](https://www.kegg.jp/entry/2.7.1.54) cites only Volk 1962, J Biol Chem 237:19–23 (PMID 13926593 as given by KEGG; the PubMed page was not opened; *Propionibacterium pentosaceum*). The full text (doi:10.1016/S0021-9258(18)81352-7) returned HTTP 403, then an Elsevier redirect stub | **unread**: no claim about its reversibility data |
| Organism | R2041 has no GPR. No *Y. lipolytica* D-arabinokinase is known to this project | verified (model); "not known" is not proof of absence |
| Thermodynamics | [eQuilibrator](https://equilibrator.weizmann.ac.il/search?query=C00002+%2B+C00216+%3D+C00008+%2B+C01112) (KEGG C00002 + C00216 = C00008 + C01112): ΔrG'° = ΔrG'm = −24.7 ± 4.7 kJ/mol, K'eq = 2.1 × 10⁴ at pH 7.5, pMg 3.0, I = 0.25 M, 298.15 K. This is a component-contribution estimate, not a measurement | verified (page read 2026-10-19 23:46 UTC) |
| Conditions (derived) | Reverse requires Q = [ADP][A5P]/([ATP][Ara]) > K'eq. K'eq at −24.7 / −20.0 / −15.3 kJ/mol is 21,245 / 3,190 / 479 (RT = 2.479 kJ/mol). For any ATP/ADP ≥ 1, reverse needs A5P/Ara > 479 even at the least favourable 2σ value | derived; *Y. lipolytica* ATP/ADP, A5P and D-arabinose levels were not looked up |
| Mechanism vs biology | The mechanism (ATP loop) is absent (§3–4). Biology: forward-only is supported by thermodynamics and by the kinase's definition. Physiological use of the reverse direction is unsupported | inferred from the rows above |

Not used as evidence: the reaction name, other models' bounds, the selection's choice, and any
match-rate effect.

## 6. Plan for approval

### 6.1 Options

| Option | What changes | For | Against |
| --- | --- | --- | --- |
| **A (recommended). Post-selection direction step, R1931/R1159 pattern** | New `model/curation/direction/r2041_direction.json` plus a guarded `apply_r2041_direction`, run after the selection. Bounds only | Same mechanism as the two existing reviewed direction fixes the selection would otherwise undo. Leaves the audited selection file (`abfc4a94…`) unchanged. Keeps the balanced 2-H⁺ equation | One more step and file |
| B. Selection amendment, R_NTP1 pattern (2026-09-14) | In `metadata_reaction_selection.json`: R2041 `fields` → `["stoichiometry"]`, an `excluded_metadata_fields` row for its bounds, an `amendments` record with the previous SHA `abfc4a94…`, `field_counts.bounds` 6 → 5 | Same final model state | Changes the selection SHA that the 2026-10-19 reviews pin. Always on. Splits the provenance: the bounds come from a patch note the selection marks "previous" |
| C. Extend `selection_gpr_restore.py` (2026-10-19 pattern) | New change set in its registry | Reuses the registry | That module never writes bounds (it only compares them), so it needs new code. It is uncommitted in two worktrees, so merges could conflict. It is named and noted as a GPR restore |
| D. Keep the selection state | Nothing | No work | A thermodynamically disfavoured reverse direction and an internal loop stay available |
| Not offered: re-apply the patch's 1-H⁺ equation | — | — | Re-creates the H −1 / charge −1 residual |

### 6.2 Option A in detail

These are the steps of the `curation-change` skill. Nothing below has been done.

1. **Data:** `model/curation/direction/r2041_direction.json`, laid out like `r1931_direction.json`.
   - `curation_id`: `R2041-FORWARD-<date>`.
   - `status`: `user_authorized`, with the user's words quoted.
   - `before_bounds` [−1000, 1000], `after_bounds` [0, 1000].
   - The exact 2-H⁺ `stoichiometry`, the five `species` (formula, charge, compartment), `gpr` "",
     `ec_code` 2.7.1.54.
   - `evidence`: the ENZYME, Rhea and eQuilibrator URLs; `review`: this report.
   - `notes`:
     - `r2041_direction_curation`;
     - `r2041_previous_bounds`;
     - `r2041_direction_evidence`: "EC 2.7.1.54 kinase direction; eQuilibrator ΔrG'° −24.7 ± 4.7
       kJ/mol, K'eq 2.1e4";
     - `r2041_direction_limits`: orphan reaction, no *Y. lipolytica* enzyme, Volk 1962 unread, the
       patch's loop rationale not reproduced; the 1-H⁺ equation in the previous note applies only
       after ATP/ADP microspecies migration, and the proton count must then be re-reviewed;
     - `r2041_direction_authorization`.
2. **Code:** a guarded `apply_r2041_direction(model, spec=None)`, as in `r1159_direction.py`.
   - Fail closed unless all of these match: schema, ID, status, reaction ID, name, EC, empty GPR,
     bounds in {before, after}, exact stoichiometry, exact species, and no conflicting notes.
   - It sets bounds and notes only, and a rerun is a no-op.
   - A species guard makes it fail closed if ATP/ADP are migrated, so the proton question cannot be
     skipped silently.
3. **Wiring:** after `apply_r1159_direction` in `main.py`. Two choices for the user:
   - (i) **opt-in** flag (for example `--r2041-forward`), recommended for now. Builds without the
     flag keep reproducing a56f3649, c98961e9 and 51c6fabd byte for byte.
   - (ii) **always on**, like R1931 and R1159. Every later build changes.
4. **Focused test** on an in-memory model:
   - applies on the selection state;
   - is idempotent;
   - fails closed on stale bounds, stoichiometry, species or note;
   - flag off leaves the model unchanged.
5. **Build** a new file, for example `model/candidates/E5_coq9_resp_pkg_r132_r2041_forward_<date>.xml`.
   - Command: a56f3649's recorded command plus the flag.
   - Control: the same build without the flag must reproduce `a56f3649…`.
6. **Whole-model diff** against a56f3649: only R2041's lower bound, its reversibility marker and its
   new notes.
7. **Then:** `pytest` compared test by test, the authorized LPs (§6.3), an independent read-only
   audit, REPORT.md and a STATE.md entry. No commit or push without an explicit instruction.

**Base (the user's choice).**
- Building in this worktree on `285eaee` gives a one-change diff against a56f3649. This is the
  recommended base.
- Stacking on the uncommitted GPR fixes (51c6fabd, `clever-kare-7d2e1a`) waits for that work to be
  committed or merged.

### 6.3 LP checks (need authorization; none run)

Context:
- The context of `coq9_respiratory_package_20261008/run_controls.py`: the
  `coq_candidate_checks_20261005/manifest.json` config, PO1f SD-Leu, Gurobi with 1 thread, Presolve 0,
  tolerances 1e-7, 60 s per solve.
- M0 = a56f3649 (R2041 reversible). M1 = the new candidate (forward-only).

| # | LP | Model | Expected | Why |
| --- | --- | --- | --- | --- |
| 1–2 | Closed-exchange max ATP dissipation, C_cy (xMAINTENANCE) | M0, M1 | 0, 0 | The asked-for test with (M0) and without (M1) R2041 reverse flux. M0 is provably ≤ 0 already. This reproduces it directly |
| 3–4 | The same, C_mi drain | M0, M1 | 0, 0 | Same |
| 5–6 | WT growth | M0, M1 | equal (inferred) | The change must not alter WT |
| 7–8 | FVA of R2041 (min, max) at 100 % WT | M0 | min < 0 would show reverse flux is in the optimal space | Whether the reverse direction can be used at optimum |
| 9–10 | FVA of R2041 (min, max) at 100 % WT | M1 | expect 0, 0 (blocked, inferred) | Documents that the reaction is inactive after the change |
| P1–P2 (optional) | Closed C_cy / C_mi without the PO1f overlay | M0 | 0, 0 | Closes the R612 gap in §4 |

A "closed ATP test with forced reverse R2041" is **not** proposed. Fixing R2041 < 0 only shrinks the
feasible set, so its objective is ≤ 0 by construction. Its feasibility is already shown statically
by the zero-sum cycle in §3.

- **Budget:** 10 LPs, or 12 with P1–P2.
- **Stop conditions:**
  - any closed objective > 1e-7: stop and report the witness flux;
  - WT M0 ≠ M1 beyond 1e-9 relative: stop and report R2041's role;
  - a non-optimal status: record the raw value. It is not biological evidence.
- **Optional, separate authorization:** a full single-KO screen of M1 (1,074 LPs), compared with the
  recorded a56f3649 screen at 1/5/10/15 %.
  - Purpose: to report call changes, not to justify the change.
  - Negatives are the user-defined class (model genes not on the experimental list). The screen is
    a development reference, not independent validation.

**Acceptance criteria:**
- the closed tests are 0 on M1;
- WT is unchanged;
- the whole-model diff is limited to R2041;
- the focused test and the unchanged suite pass;
- the audit finds no undeclared change.

## 7. Limits and items outside this review

- **R2190 EC versus product (affects another line; candidate).** R2190 carries EC 1.1.1.250, and its
  GPR is YALI1F03452g. ENZYME [1.1.1.250](https://enzyme.expasy.org/EC/1.1.1.250) is
  "D-arabinitol + NAD(+) = D-ribulose + NADH + H(+)" (verified). The model writes D-arabinose as the
  product.
  - Either the product or the EC annotation is wrong. Which one has not been checked: the local entry
    may point to a different enzyme.
  - R2190 is the step that links D-arabinose to the arabinitol node, and so to the internal loop
    in §3.
  - Its source audit is in the colleague-owned `model/curation/lipid_unlump/` file; review and
    recommendation only.
- **R2041 and R2042 are both orphans (affects the orphan line).** ENZYME
  [5.3.1.13](https://enzyme.expasy.org/EC/5.3.1.13) places the A5P isomerase in the Kdo pathway of
  bacterial lipopolysaccharides and plant cell walls (verified). Whether either reaction exists in
  *Y. lipolytica* was not examined. Removing them would need evidence that they do not exist, which
  is not part of this proposal.
- **Not covered:**
  - knockout-specific use of reverse R2041 (only the optional screen would show it);
  - builds without the E5 chain;
  - measured *Y. lipolytica* metabolite ratios.
- **Interpreter:** the main checkout's `.venv` (Python 3.13.5, cobra 0.30.0). Gurobi licence
  messages appeared at load; no solve was run.

## 8. Commands and scripts

The scripts are in the session scratchpad, not in the repository. They can be added under
`scripts/` here on approval.

| Script | What it does |
| --- | --- |
| `trace_r2041.py` | R2041 fields, annotation, notes and species |
| `trace_net.py` | Arabinose-related metabolites |
| `trace_rxns.py` | Every reaction of the listed metabolites, with bounds, GPR, EC and element/charge residual |
| `compare_models.py` | Whole-reaction comparison of `0ac2addf` and a56f3649 |
| `cycle_sum.py` | The summed stoichiometry of the cycle, and R2041 in c98961e9 and 51c6fabd |

External pages were saved with `curl` to the scratchpad:
- `enzyme_2.7.1.54.txt`, `enzyme_1.1.1.250.txt`;
- `rhea_24588.tsv`, `kegg_ec_2.7.1.54.txt`;
- `eq_r01573_kegg.html`.

The recorded-flux scan was an inline Python walk over `model/reports` for the pattern
`"R2041": <number>`.

## 9. Proposed STATE.md entry (not written)

Superseded the same day: the user approved Option A, and the combined review and implementation entry is
in this worktree's STATE.md (see [r2041_forward_20261019](../r2041_forward_20261019/REPORT.md)). The text
below is kept as written.

Not written. Several worktrees already hold diverging uncommitted STATE.md heads (see the
`selection_gpr_review_20261019` hand-over). Text to add when merged:

> ## 2026-10-19: R2041 (D-arabinokinase) direction review, investigation only
>
> Delegated read-only review; no model, data or code change, no LP. The 2026-09-10 metadata selection
> reverts the builder patch `fix_d_arabinokinase_direction_and_proton` on R2041 to [−1000, 1000] with
> 2 H⁺ (balanced with stored species; the patch's 1 H⁺ leaves H −1/charge −1 until ATP/ADP microspecies
> migration). In a56f3649 the patch's "isomerase/phosphatase ATP loop" does not exist: A5P has only
> R2041/R2042, and the one closed cycle through reverse R2041 (R2190, R1964, R803, R714, R2042) sums to
> zero (verified). The recorded closed-input ATP test on 0ac2addf (= a56f3649 + R132) is 0 in C_cy and
> C_mi, so a56f3649 is 0 (provable). Forward-only is supported by thermodynamics (eQuilibrator ΔrG'°
> −24.7 ± 4.7 kJ/mol, K'eq 2.1e4; verified). Plan awaiting approval: bounds-only post-selection step
> `R2041-FORWARD-<date>` (R1931/R1159 pattern), 10 LPs. Side finding: R2190 EC 1.1.1.250 product is
> D-ribulose, model writes D-arabinose (candidate). Report:
> [r2041_direction_review_20261019](reports/r2041_direction_review_20261019/REPORT.md).
