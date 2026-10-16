# Task D: complex I transhydrogenation module (holistic combinatorial test)

Written 2026-10-13, before any scored run. Investigation and diagnostics only; no model change.

## Authorization (user, 2026-10-13)

- "now do the task D": Task D of `claude_code_prompt_w29_tasks_20261012.md` (header §0–2 plus Task D),
  with the Task C handoff (`model/STATE.md` 2026-10-13; Task C REPORT §12.4).
- Plan approved in the conversation, with these decisions:
  - **Grouping:** hierarchical. Tier 1 is 2^k over at most 12 mechanism groups. Tier 2 is a 2^m test of the
    members of each group that sits in a minimal cut set, with the other groups of that set held closed.
  - **Checkpoint:** after discovery, freeze the part table here and continue into the 2^k runs. Stop and ask
    only if a group is a judgment call rather than one clean mechanism.
  - **Stop rule:** if every minimal cut set costs more than 2 % of WT growth, stop and report. The evidence
    for the parts and the maps are finished; step 5 (13C ratios) waits for the user.
  - **Step 4:** any capacity or regulatory assumption is proposed in the report; the ASSUMPTION condition
    file is written only after the user chooses.
- **Allowed (prompt §0):** local cobrapy solves, new condition files (only as above), reports under
  `model/reports/`, a `model/STATE.md` entry.
- **Not allowed:** building or editing any model XML; editing existing media, evidence or benchmark files;
  deleting files; git commit, push, merge or branch switch; cluster jobs; any change made to improve a
  match rate.

## Starting state (verified 2026-10-13)

| Item | Value |
|---|---|
| Worktree / HEAD | branch `claude/task-d-78f16d` @ `4d97b828ac2b0dc95d592060605af75069705c09` (main with PRs #4 and #5), clean |
| Model | `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` (Task B candidate, R132 removed), SHA-256 `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`, tracked at HEAD |
| Task A medium | `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv`, `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| 13C caps (ASSUMPTION, diagnostic) | `model/conditions/assumptions/w29_13c_upper_caps_20261012.json`, `df1e33a8aa72c2d1b83129cb6a0e0d5e6eb283e3aec8feeecd844bcc2dc4257b`; "module" glycine-cleavage form only (Task C Amendment 1) |
| R825 reversible (ASSUMPTION) | `model/conditions/assumptions/r825_reversible_20261013.json`, `e45cab424c30ee9aec00823cb67392742ac6c36ba7091664249765a2533efcaa`; step 5 only |
| Overlay (diagnostic) | `holistic_overlay.py`, copied byte-identical from Task C (`bf6c7ce54d8338c365d98c82ca345ff2bd8c6055b38d402e8144a16b4abf30c4`) |
| PO1f profile / experimental list / mtDNA scope | `35307853…` / `1e887f5a…` / `16b3e9b4…` (as in Task C) |
| W29 MFA | `$IYALI26_RESEARCH_ROOT/reference/w29_worland2024/glucose_central_fluxes.csv` `844615d0…`; repository ↔ paper link **inferred, pending user confirmation** |
| Research root | `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`, passed on every command (the session hook reports the variable unset) |
| Python | main checkout `.venv/bin/python` (cobra 0.30.0, Gurobi 13.0.1, Python 3.13.5) for LPs; conda Python with matplotlib for maps |

Every LP script checks these SHAs before it runs.

**Pre-plan exploration (disclosed).** Before the plan, about 1,050 LPs ran in the session scratchpad on Task C's
helpers. They reproduced the escape qualitatively and showed multi-step routes beyond the pair list. None of
those numbers is used as a result; the scripts below re-run everything.

## Configurations

| Name | Medium | Overlay | Caps | Use |
|---|---|---|---|---|
| **A** (primary) | Task A medium + PO1f runtime profile, glucose uptake ≤ 2.414 | holistic (diagnostic) | 13C module form | steps 1–4 |
| **V2** (sensitivity) | W29-type defined V2, built in memory as in Tasks B/C, glucose fixed at 2.414 | holistic | 13C module form | steps 1–4 |

- Objective `biomass_C`. Complex I KO: R2062 bounds (0, 0).
- The caps are WT-state 13C bounds. They are applied to WT and KO alike, as in Task C §5.2; this is a
  diagnostic, not a screen condition.
- "R132 off" is the base model itself (Task B removed R132).

## Step 1: escape (fixed before running)

`escape.py`, configurations A and V2, WT and KO:
- growth, pFBA vector and FVA at the optimum (fraction 1.0);
- node balances on cytosolic NADH (m123), cytosolic NADPH (m40) and matrix NADH (m30), plus matrix NADPH
  (m176): every reaction's signed contribution (pFBA) and its FVA range, per 100 glucose;
- a map of the escape with the prompt's nodes.

## Part list (rules fixed before running)

`parts.py` builds the round-0 pool, in configuration A:
1. the pair scan: same compartment, identical stoichiometry apart from NAD(H) vs NADP(H), H+ ignored, same
   or opposite orientation; run on the base model and again after the overlay;
2. every reaction carrying NAD(H) or NADP(H) in two compartments (transport);
3. reactions that hold both an NAD and an NADP species (kinases), reported with the pool;
4. the earlier parts R361, R1931, R1936, R328, R691 and R1310;
5. the handoff parts: {R412, R1172} as one part, and R1708 as its own part.

**Pruning.** A reaction is dropped if FVA at any growth ≥ 0 (fraction 0) gives [0, 0] in both WT and KO, or
if the overlay closes it. A pair is dropped if either member is dropped: neither a cycle nor a cofactor
choice is then possible.

**Closure = escape direction only (diagnostic).** For each reaction, compare its FVA range at the KO optimum
with its range at the WT optimum (same configuration, same earlier closures):
- the escape direction is a direction the KO optimum can use (|bound| > 1e-6 per mmol gDW⁻¹ h⁻¹) and the WT
  optimum does not require;
- the closure sets that direction to 0. Example: R550 reverse, R1385 reverse; R549 and R329 whole because
  their only direction is the escape one;
- a pair is closed when every member that has an escape direction is closed in that direction;
- a reaction whose KO use is in a direction the WT also requires has no escape direction. It is a
  "capacity" step: it is not closed in the main test; its whole-reaction closure cost is reported;
- **exceptions, fixed by the user:** {R412, R1172} and R1708 are closed whole, as in Task C's module test;
- no enzyme is removed from any XML.

**Discovery rounds** (`discovery.py`; at most 6):
- Close every part found so far (WT and KO). Solve the KO and write the four node balances, plus the list of
  reactions that carry ≥ 5 per 100 glucose in the KO pFBA vector in a direction the WT optimum does not
  require.
- Name the route that took over (by stoichiometry and compartment). Add it as a new part with the closure of
  one or more of its WT-free steps, choosing the step with the smallest single-closure WT cost.
- Earlier parts and pool reactions with no escape direction in any round are dropped as "not relevant",
  with their FVA shown.
- Discovery stops when: the KO/WT with every part closed is < 0.15; or the KO escape has no WT-free step
  carrying ≥ 5; or 6 rounds have run.

**Grouping.** Parts with one mechanism form a tier-1 group (for example, "cytosolic NAD/NADP identical-
stoichiometry pairs"). k ≤ 12 groups. If a group would not be one clean mechanism, stop and ask.

The frozen part table is appended below before the 2^k runs.

## Step 2: combinatorial test (fixed before running)

`combinatorial.py`:
- **Tier 1:** every subset S of the k groups, configurations A and V2. Two LPs per subset (WT and KO, each
  with S closed).
  - KO/WT is the KO growth over the WT growth under the same closures.
  - WT cost = 1 − WT(S) / WT(none closed), same configuration.
  - A **hit** is: WT feasible, and KO/WT < 0.15, or KO with no solution. A KO with no solution is reported
    separately as **unresolved, not lethal** (benchmark rule).
  - A **minimal cut set** is a hit with no proper subset that is a hit. KO/WT is not monotone in S, because
    the WT also falls, so every minimal hit is listed.
  - **Accepted** = minimal cut set with WT cost ≤ 2 %.
- **Leave-one-out:** from the full set (all groups closed), reopen each group singly, then each part singly.
- **Tier 2:** for each minimal cut set (at most the 5 cheapest by WT cost, configuration A), and each group in
  it with m ≥ 2 parts: all 2^m subsets of its parts, with the other groups of the set closed. Report the
  minimal part sets that keep the hit, and their WT cost. V2 re-checks the resulting part sets only.

## Step 3: evidence

For every part in a minimal cut set (and for the enzymes the prompt names: MTD1, ADE3, GDH1/GDH3/GDH2,
ALD6/ALD4/ALD5, IDP1, IDH1/IDH2), record:
- the *Yarrowia* gene(s) from the GPR, and the name/function per AGENTS.md §4;
- the *S. cerevisiae* ortholog by local Smith-Waterman (Task B's implementation, SHA-pinned), with
  characterisation papers opened before they are cited;
- compartment, physiological direction and thermodynamics, kept as separate questions (prompt rule 3).

Also check *S. cerevisiae* Erg13 localization against a primary source (Task C handoff; a comparison, not
*Yarrowia* evidence). Literature search may use a read-only subagent; every DOI is opened before it is cited.

## Step 4

If complex I essentiality needs a real enzyme closed or a direction reversed without evidence, the report
states the capacity or regulatory assumption it needs and proposes the ASSUMPTION condition file. The file is
written only after the user chooses.

## Step 5 (only if steps 2–4 end without a stop)

13C ratio floors × glucose uptake: PDH R693 ≥ 0.866629, CS R236 ≥ 0.282, SDH R740 ≥ 0.151, G6PDH R325 ≥ 0.503
(the exact MFA lower bounds are read from the pinned MFA file). All 16 subsets, R825 as supplied and with the
ASSUMPTION file (every such result labelled), WT and KO, configurations A and V2. A ratio set is accepted only
if WT loses ≤ 2 %.

## Maps

PNG, one panel per condition, a value on every edge (per 100 glucose, pFBA vector), removed steps marked:
- the step-1 escape node map;
- WT, complex I KO, and the KO under each accepted set (or, if none is accepted, the cheapest hit, labelled
  "not accepted").

Nodes: G6P, oxPPP, NADPH (cyt), MTHFD pair, NADH (cyt), NDH2, pyruvate, PDH, TCA, NADH (mito), complex I, CoQ,
plus any route that enters a minimal cut set.

## Budget and stop conditions

- **Budget:** about 25,000 LPs; hard stop at 40,000 (each script records its LP count).
- **Stop if:**
  - any input SHA differs from the table above;
  - WT is infeasible in a configuration;
  - k cannot be kept ≤ 12 with clean mechanism groups (ask);
  - every minimal cut set costs more than 2 % of WT growth (stop and report; step 5 waits);
  - the LP budget is reached.

## Deliverables

TASK.md, REPORT.md (claims tagged verified / inferred / recalled), CSV tables, maps (PNG), scripts, an
independent read-only audit, a `model/STATE.md` entry. Large outputs go to
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/complexI_transhydrogenation_20261012/`. No commit or push.

## Frozen part table (written 2026-10-14, after discovery, before any 2^k run)

### Discovery outcome (configuration A; verified, `discovery_round1..4.json`)

| Round | Parts closed (cumulative) | WT | KO | KO/WT |
|---|---|---:|---:|---:|
| 0 (no closure; `escape.json`) | 0 | 0.282216 | 0.233223 | 0.8264 |
| 1 | 9 parts, 12 reactions | 0.280795 | 0.151502 | 0.5395 |
| 2 | 15 parts, 18 reactions | 0.280795 | 0.122409 | 0.4359 |
| 3 | 17 parts, 20 reactions | 0.280795 | 0.114031 | 0.4061 |
| 4 | 22 parts, 25 reactions | 0.280795 | 0.021640 | 0.0771 |

Discovery stopped at round 4 by its first rule (KO/WT < 0.15). The 0.50 % WT cost is the HMG-CoA part,
closed whole by the user's instruction.

**Rule refinements used during discovery (disclosed).** Every candidate step tied at WT cost 0
(`step_costs_round1..3.json`), so "smallest WT cost" never decided. The tie-breaks I applied, in order:
1. an earlier or pool part that is a WT-free step of the route (R328, R1931);
2. the single step that feeds every branch of the route (R2065);
3. the step that makes or consumes the cofactor (R1968, R2230, R2119, R533 reverse);
4. the first committed step (R244, R2281);
5. an exchange over its transport (R1148).

**Dropped as not relevant** (pool reactions with no escape direction in rounds 0–4): R134, R330, R691, R361,
R135, R2070, R487, R738, R1264, R1265, R1593, R1594. Dropped pairs (a member blocked or closed by the
overlay): R130/R1377 and R_ALCDH_nadp_hi, R2069/R131, R254/R255, R534/R_MDHy, R737/R738. These are capacity
steps or idle; they are not closed in the test.

**Added after the last discovery run:** `pyruvate_secretion` (R1310, earlier part, escape direction first
seen in round 4) and `NADPH_transport_em` (R1590/R1591, FVA-only escape direction in round 4).

### Amendment 1 (user, 2026-10-14): grouping

Group 10 is a family of seven routes, not one mechanism. The user chose to keep 12 tier-1 groups, with
group 10 labelled a family and resolved by tier 2 (2^7) inside any cut set it enters.

### Tier-1 groups (k = 12) and their parts

Closures are in `parts_definition.json` (SHA-256 recorded in the combinatorial output).

| # | Group | Parts (closure) | Mechanism |
|---|---|---|---|
| G1 | cyt_NADP_NAD_pairs | MTHFD (R549 fwd, R550 rev); GDH (R329 fwd) | identical-stoichiometry NAD/NADP pairs forming cytosolic NADPH → NADH cycles |
| G2 | ALD_cofactor_choice | ALD_cyt (R1378 fwd); ALD_mito (R136 fwd) | acetaldehyde dehydrogenase cofactor choice (cytosol: NADH instead of NADPH; matrix: NADPH instead of NADH) |
| G3 | matrix_citrate_import | citrate_import_mito (R817 rev); citrate_import_mito_831 (R831 fwd) | carrier-mediated citrate import into the matrix, feeding matrix NADP-IDH |
| G4 | matrix_IDH_reductive | IDH_mito_reductive (R1385 rev) | matrix NADP-IDH reductive carboxylation |
| G5 | matrix_NADH_malate_export | matrix_MDH_reductive (R533 rev) | matrix NADH leaves as malate (matrix MDH reverse) |
| G6 | R1708 | R1708 (whole) | succinyl-CoA:acetate CoA-transferase (user: own part) |
| G7 | HMG-CoA_route | R412 + R1172 (whole) | matrix acetyl-unit export as HMG-CoA (user: one part) |
| G8 | serine_disposal | serine_cycle (R511 fwd); serine_transsulfuration (R244 fwd) | serine → pyruvate disposal closing the serine-pathway NADPH → NADH cycle |
| G9 | intermediate_secretion | glyoxylate_secretion (R1148 fwd); pyruvate_secretion (R1310 fwd) | secretion of a partly oxidised intermediate through a no-GPR exchange |
| G10 | aa_amine_catabolism (**family**) | GABA_shunt (R328); GSA_P5C (R1931); leucine_degradation (R1968); isoleucine_degradation (R2281); polyamine_catabolism (R2065); malonate_semialdehyde_cycle (R2119); R1936 | amino-acid/amine synthesis–degradation routes that make cytosolic or matrix NAD(P)H |
| G11 | mito_FAS_cycle | mito_FAS_cycle (R2230 fwd) | matrix fatty-acid synthesis + peroxisomal β-oxidation |
| G12 | cofactor_transport_membranes | NADP_transport_lp (R1628 fwd, R1629 rev); NADPH_transport_em (R1590 fwd, R1591 rev) | no-GPR NADPH/NADP+ exchange between cytosol and a membrane compartment |

Tier 1 runs on A and V2 with this table. Tier 2 resolves the parts of each multi-part group inside the five
cheapest minimal cut sets (configuration A), and V2 re-checks the resulting part sets.

## Post-audit runs (2026-10-14)

The independent audit found two gaps, and both runs below were added to close them:
- **Capacity-step closure costs:** this file required them, but they had not been run → `capacity_costs.py`.
- **Shuttle-cap thresholds:** they had been read from a coarse grid → `shuttle_cap.py fine`.

Both are diagnostics. No pre-registered rule or result changed.

## Amendment 2 (user, 2026-10-14): step 4 decision

The user chose **option (d)** of REPORT §7. No ASSUMPTION condition file is written. The record stands that the
model needs a matrix-to-cytosol NADH-export assumption to reproduce complex I lethality. Not authorized here:
- options (a)–(c);
- any model change;
- the follow-up evidence review of R831's citrate/malate mode.
