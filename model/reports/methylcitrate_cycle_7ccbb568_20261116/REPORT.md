# Methylcitrate cycle on 7ccbb568: R_OBDHm added, R2075 closed (2026-11-16)

Scope, authorization, starting state and acceptance criteria are in [TASK.md](TASK.md). The audit is in
[AUDIT.md](AUDIT.md). Scripts are in [scripts/](scripts/), logs in [logs/](logs/) and the figures in
[figure/](figure/). Large outputs are in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/methylcitrate_cycle_7ccbb568_20261116/`:
- the control and seed-1 builds;
- the sweep TSV;
- the pytest junit files and logs.

Nothing is committed or pushed.

## Answer in short

- **Built.** `model/candidates/E5_7ccbb568_methylcitrate_20261116.xml`, SHA-256
  `ad901ddb4a2cfe5b09eb29f9c4b6630c0ee77cf93cef5be5aff3c7dfc954c279`. It is built from baseline `7ccbb568` with one new
  opt-in flag, `--methylcitrate-cycle`, and changes only two things:
  - **R_OBDHm (new):** matrix 2-oxobutanoate + CoA + NAD → propionyl-CoA + CO2 + NADH + 4 H+, bounds (0, 1000). Its
    gene rule is copied from R693 (pyruvate dehydrogenase) and marked INFERRED.
  - **R2075:** the gene-less methylmalonyl-CoA mutase, closed to (0, 0) and kept in the model.
- **Verified this round (reproduced from the built files):**
  - **Builds.** The control build is byte-identical to 7ccbb568. The two candidate seeds are byte-identical.
  - **Diff.** The whole-model diff has 9 entries, all on the targets. R_OBDHm is mass- and charge-balanced.
  - **WT growth.** 0.31110053592594944 (7ccbb568: 0.311100535925948).
  - **Blocked set.** R25, R95 and R552 now reach 7.5987 (Task A, no growth floor). R2074, R2075, R2076, R2077 and
    R2112 stay at 0/0. R834 reaches 3.3108. The route observed for it is side route B, which is left as is; that it
    is the only route is inferred.
  - **Small capacity rises elsewhere.** The WT maximum maintenance flux rises from 63.333 to 63.372. The maxima of
    R693, R412, R741, R83 and R410 rise by less than 0.1 %.
  - **Energy check.** 0 for both `R_NTP1` and `biomass_C`.
  - **Full certified knockout sweep (1,081 genes).** No gene changes status or class at 1, 5, 10 or 15 %, under
    rule A or B. The scores are identical to the recorded baseline sweep. The largest change in knockout growth is
    1.4e-14.
- **Not settled by this round:**
  - the rule of R_OBDHm (inferred);
  - the absence of a methylmalonyl-CoA mutase (inferred from annotation and homology);
  - side route B.

## 1. What changed and how it was built

| Item | 7ccbb568 | ad901ddb |
| --- | --- | --- |
| R_OBDHm [C_mi] | absent | m42 + m73 + m27 → m78 + m44 + m30 + 4 m28, (0, 1000); rule `YALI1E31947g and YALI1F27556g and YALI1D26360g and YALI1B13179g and YALI1D30618g` (= R693); subsystem and group "Propanoate metabolism" (`group60`); SBO:0000176, Rhea 33223; four notes, one saying the rule is INFERRED |
| R2075 [C_cy] | (0, 1000), no notes, no gene | (0, 0), three notes: why, previous bounds, closure evidence (INFERRED) |
| The five PDH genes | as in [before.json](before.json) | each also linked to R_OBDHm |

**Source of the step.** The step and its data are byte-identical copies from the other line's methylcitrate-cycle round
([TASK.md §1](TASK.md)):
- `platform/scripts/gem_annotate/methylcitrate_cycle.py` (`79dabdd0…`);
- `model/curation/methylcitrate_cycle/methylcitrate_cycle_20261116.json` (`17367e6f…`).

Three fields therefore name the other line: `baseline_model_sha256` (cbf2ef95), `report` and the module docstring.

**Code added on this line:**
- the `--methylcitrate-cycle` flag in `cli.py`;
- the parameter, call and pass-through in `main.py`; the call runs after `--reaction-xref-cleanup`, last of all;
- `platform/tests/test_methylcitrate_cycle_7ccbb568.py`: the other line's eight checks, run on the tracked
  7ccbb568 file.

**What the step pins before any edit:**
- the two targets;
- R693;
- the nine species;
- four species' reaction lists;
- the five PDH genes' reaction links;
- the group;
- the absence of any reaction with the same chemistry.

If any of these differs, it raises before making any edit. It also raises when only one of the two edits is
present. A rerun on the built file returns "already_correct".

**Build command** (from `platform/`; the full argv is in [scripts/args_7ccbb568.txt](scripts/args_7ccbb568.txt), driven
by [scripts/run_builds.zsh](scripts/run_builds.zsh)):

```bash
PYTHONHASHSEED=0 /opt/miniconda3/bin/python -m scripts.gem_annotate <the 47 arguments of 7ccbb568's .build.json> \
  --methylcitrate-cycle --output-model "$WT/model/candidates/E5_7ccbb568_methylcitrate_20261116.xml"   # WT = worktree root
```

| Build | Seed | Output SHA-256 | Check |
| --- | --- | --- | --- |
| Control (no new flag) | 0 | `7ccbb56834a37ad4c75e434868e9e1b9e1f6a657a11ca17ee0a93dbe56232146` | `cmp` identical to the tracked 7ccbb568 |
| Candidate | 0 | `ad901ddb4a2cfe5b09eb29f9c4b6630c0ee77cf93cef5be5aff3c7dfc954c279` | formal file in `model/candidates/` |
| Candidate | 1 | `ad901ddb…` | `cmp` identical to seed 0 |

The build record `E5_7ccbb568_methylcitrate_20261116.build.json` shows:
- HEAD `af7c5432`;
- `input_unchanged` and `source_unchanged_during_build` both true;
- 2,348 reactions, 1,883 species, 1,081 genes;
- selection status "applied".

The dirty status lists only:
- this round's files;
- the r1172 investigation folder;
- `model/STATE.md`.

All three builds log one libSBML error: "Setting notes on sbase" for R304. That error comes from the 7ccbb568
argv, not from this step: the control build, which reproduces 7ccbb568 byte for byte, logs it too.

## 2. Whole-model verification

[scripts/model_diff.py](scripts/model_diff.py) compares every reaction, metabolite, gene, group, compartment and the
model record. It gives 9 entries ([model_diff.json](model_diff.json)):
- R2075: bounds (0, 1000) → (0, 0); notes {} → the three curation notes.
- R_OBDHm: present only in the candidate.
- `YALI1E31947g`, `YALI1F27556g`, `YALI1D26360g`, `YALI1B13179g`, `YALI1D30618g`: reaction list + R_OBDHm, nothing
  removed.
- `group60` ("Propanoate metabolism"): members + R_OBDHm, nothing removed.

No metabolite (formula, charge, compartment, annotation, notes), compartment, objective coefficient or other reaction
changes.

The text diff of the two XML files has 3 hunks, 51 lines
([logs/xml_line_diff_vs_7ccbb568.txt](logs/xml_line_diff_vs_7ccbb568.txt)):
- the R2075 line (upper bound `cobra_0_bound`) and its notes;
- the new R_OBDHm block;
- one group member.

**Balance.** `check_mass_balance()` on the built R_OBDHm returns {} (mass and charge). The 4 H+ follow the model's
storage of propionyl-CoA m78 as the tetra-anion (C24H36N7O17P3S, −4) next to neutral CoA m73 (C21H36N7O16P3S, 0). This
is a protonation convention of the model, not chemistry; the reaction's notes say so. No formula or charge changed, so
no other reaction's balance changed.

## 3. Mechanism and evidence

The mechanism is the investigation's ([../r1172_side_effect_20261116/REPORT.md](../r1172_side_effect_20261116/REPORT.md)
§3). This round re-solves the LP results on the built files; the static scan was not re-run in this round's scripts.

**Why the nine were blocked:**
- **Proven (static identity, investigation).** No matrix reaction makes or destroys the CoA moiety. 34 matrix
  reactions in 7ccbb568 conserve it and 9 transports move it. The audit's re-scan finds 35 in ad901ddb; the new one,
  R_OBDHm, also conserves it. Net import of CoA moiety into the matrix is therefore 0 at every steady state.
- **No exit after the closure** (investigation; observed in the LP with two solvers). After R1172 was closed, no
  reaction could take CoA moiety out of the matrix:
  - R2077, R2083 and R1291 can only import (statically verified);
  - R2084 and R2118 are blocked (LP);
  - R2025 only passes on propionyl-CoA that came in on R2113 (LP).
- **Each of the nine needs CoA moiety to come in:**
  - R25 frees imported propionyl-CoA's CoA;
  - the methylmalonyl route ends in a succinyl-CoA import;
  - R834 imports CoA.

**Why R_OBDHm frees R25, R95 and R552.**
- Observed in the LP this round; the reasoning is a consequence of the identity.
- R_OBDHm makes propionyl-CoA from matrix CoA, and R25 returns that CoA. The loop closes inside the matrix, with no
  transport and no CoA exit.
- In the forced route (R25 = 1, minimum total flux, [figure_flux.json](figure_flux.json)) R_OBDHm = 1.00. R834,
  R2025, R2113 and R2077 are all 0.

**Why R2075 is closed.**
- With R_OBDHm alone, side route A opens (in-memory LP in the investigation, §6): R2113 out → R2112 backwards →
  R2074 → R2076 → R2075 → R2077 in.
- It runs through four gene-less steps. R2112 runs as a reductive carboxylation, against the enzyme's described
  direction.
- Closing R2075 stops it. The five methylmalonyl reactions are 0/0 in the built candidate (LP, this round).

**Side route B is left as is** (no user decision; the curation data records it). It is the route observed for R834.
The route was seen in the investigation's minimum-flux solve, and R2025's export capacity (3.3108) equals R834's
maximum on ad901ddb. That no other route exists is inferred: no R834 route solve was run on ad901ddb.
- R834 imports CoA (3.31 at most);
- R_OBDHm turns it into propionyl-CoA;
- the gene-less R2025 exports the propionyl-CoA;
- R2024 runs backwards in the cytosol (AMP + PPi → ATP), and propanoate leaves.

The energy check is 0, so this does not create free ATP. But it is a CoA exit through a transport for which no
evidence exists that acyl-CoAs cross the inner membrane this way (inferred).

**External evidence**, with tiers. No new fetch was made in this round; the tiers say where each was read.

| Claim | Source | Tier |
| --- | --- | --- |
| No *Y. lipolytica* gene for EC 5.4.99.2 (human control `hsa:4594`) | [KEGG yli link](https://rest.kegg.jp/link/yli/ec:5.4.99.2), [hsa link](https://rest.kegg.jp/link/hsa/ec:5.4.99.2) | verified in this session's investigation (fetched 2026-11-16); not re-fetched |
| *S. cerevisiae* ICL2, a 2-methylisocitrate lyase, is in the mitochondrial matrix | [PMID 11092862](https://europepmc.org/article/MED/11092862) | verified (abstract) in the investigation; re-read by the audit |
| Rhea 33223 = 2-oxobutanoate + NAD+ + CoA = propanoyl-CoA + CO2 + NADH, no EC | [Rhea query](https://www.rhea-db.org/rhea?query=2-oxobutanoate%20propanoyl-CoA&columns=rhea-id,equation,ec&format=tsv&limit=20) | verified by the other line and by this round's audit (TSV row read); not opened by me |
| 2-oxobutanoate goes to propionyl-CoA via PDH or pyruvate oxidase (authors' statement, in a *Y. lipolytica* paper) | [PMC6987463](https://europepmc.org/article/PMC/PMC6987463) | read by the audit (full text, introduction); a general pathway statement, not a *Y. lipolytica* experiment |
| Rat PDH and BCKDH both decarboxylate 2-oxobutyrate | [PMID 3718468](https://pubmed.ncbi.nlm.nih.gov/3718468/) | read by the audit (abstract); mammalian |
| *Y. lipolytica* phd1 deletion cannot grow on propionate | [PMC5991449](https://europepmc.org/article/PMC/PMC5991449) | read by the audit (full text) |
| R_OBDHm's rule = R693's PDH rule | user decision ("just infer") | inferred, no *Y. lipolytica* evidence |
| R2075 has no enzyme in *Y. lipolytica* | KEGG above; orphan gene search 2026-10-19 (no W29 hit to MMUT or MutB at E ≤ 1e-5) | inferred: absence by annotation and homology, no experiment |

## 4. Targeted solves (Task A + PO1f, Gurobi 1e-9, reset per solve, targeted check)

Script: [scripts/targeted_solves.py](scripts/targeted_solves.py). Outputs:
- [targeted_control.json](targeted_control.json): the control build, 7ccbb568, 26 LPs;
- [targeted_candidate.json](targeted_candidate.json): 28 LPs.

All solves were optimal, and WT passed the targeted check.

| Quantity | 7ccbb568 | ad901ddb | Expected (TASK §2) |
| --- | --- | --- | --- |
| WT growth | 0.311100535925948 | 0.31110053592594944 | unchanged within 1e-9: yes (Δ 1.4e-15) |
| R25, R95, R552 max (no floor) | 0 | 7.598721858878638 each | 7.5987: yes |
| R_OBDHm max | absent | 7.5987218588786725 | — |
| R2074, R2075, R2076, R2077, R2112 min/max | 0 / 0 | 0 / 0 | 0/0: yes |
| R834 max | 0 | 3.310840167653775 | 3.3108: yes |
| Energy check: max `R_NTP1`, max `biomass_C` | 0, 0 | 0, 0 | 0: yes (only `xMAINTENANCE` had a positive lower bound) |
| KO/WT of `YALI1E31947g`, `YALI1F27556g`, `YALI1B13179g`, `YALI1D30618g` | 0.98106 | 0.98106 | same class |
| KO/WT of `YALI1D26360g` | 0 | 0 | same class |

**Context from the figure solves** ([figure_capacity.json](figure_capacity.json), 48 LPs;
[figure_flux.json](figure_flux.json), 5 LPs). The largest flux rises as follows; R1172 stays at 0/0:

| Reaction | 7ccbb568 min / max | ad901ddb min / max |
| --- | --- | --- |
| R963 (2-oxobutanoate transport, gene-less) | −2.94 / 2.67 | −2.94 / 7.60 |
| R515 (threonine deaminase) | 0 / 2.94 | 0 / 7.60 |
| R2025 (propionyl-CoA transport; negative = out of the matrix) | −2.67 / 0 | −3.31 / 0 |
| R2024 (propanoate:CoA ligase; negative = run backwards) | −2.67 / 0 | −3.31 / 0 |
| R693, R412, R741, R83, R410 (largest absolute flux) | 15.802, 69.339, 65.638, 12.386, 1.7980 | 15.812, 69.387, 65.684, 12.394, 1.7993 |

The WT maximum maintenance flux (`wt_maintenance_ceiling`, sweep summaries) rises from 63.33341 (353cefdc) to 63.37226
(ad901ddb). The certificate ceilings stay at 4.828. The energy check is 0, so this is not free ATP. A likely
reason, inferred and not tested: R_OBDHm adds a matrix NADH source from 2-oxobutanoate.

- At the growth optimum (pFBA), total flux is the same in both models (284.50119). The drawn reactions carry the same
  flux, and the cycle is 0, so growth does not need it. Eight undrawn reactions (R1379, R1381, R2081, R2109, R334,
  R624, R632, R773) differ: an alternative optimum of equal total flux, not a change of route.
- In the forced route, 2-oxobutanoate comes in on R963, which already carried flux in 7ccbb568; it is not newly
  active.

## 5. Full certified single-gene sweep

[scripts/run_sweep.zsh](scripts/run_sweep.zsh) runs `ko_screen_certified.py` (`313d5d01…`, unchanged) with
`checked_solve.py` (`69c42690…`) on the candidate. Inputs and their SHAs are in
[logs/sweep_inputs_sha256.txt](logs/sweep_inputs_sha256.txt); the summary is in
[sweep_candidate_summary.json](sweep_candidate_summary.json).

**The run.** 1,081 genes, 1,123 LPs (1,082 sweep, 1 WT ceiling, 40 certificate).
- All passed the targeted check.
- 20 knockouts are infeasible, and all 20 are certified maintenance shortfalls.

**The comparison** ([scripts/compare_sweeps.py](scripts/compare_sweeps.py),
[sweep_compare_vs_7ccbb568.json](sweep_compare_vs_7ccbb568.json)) is against the recorded sweep of 353cefdc
(`sweep684_candidate.tsv`, `49b4f27b…`). 353cefdc is LP-identical to 7ccbb568 (static check in the investigation).
- Same gene set.
- No status or certificate change.
- No class change in any of the 8 class columns (rules A and B, each at 1/5/10/15 %).
- The largest change in knockout growth is 1.39e-14 (`YALI1F29751g`).

**Scores**, identical to the baseline at every cutoff and scope. Development reference: model genes not on the
consensus list count as non-essential (user-defined negatives); mtDNA genes are out of scope. Rule B, mtDNA scope:

| Cutoff | TP | FN | FP | TN | Recall | Precision | MCC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 % | 117 | 209 | 37 | 702 | 0.359 | 0.760 | 0.405 |
| 5 % | 124 | 202 | 39 | 700 | 0.380 | 0.761 | 0.419 |
| 10 % | 128 | 198 | 44 | 695 | 0.393 | 0.744 | 0.417 |
| **15 % (primary)** | **145** | **181** | **44** | **695** | **0.445** | **0.767** | **0.465** |

- **Positive coverage.** 326 experimental positives are model genes in scope (TP + FN). The
  [summary](sweep_candidate_summary.json) carries both scopes, with and without the untargeted genes.
- **Meaning.** The build changes no call. R95's and R552's listed-essential genes keep their calls. A gene is
  essential in the model when its knockout stops growth. The cycle was blocked before and is not used at the optimum
  now, so it cannot affect these calls.

## 6. Tests

Root `.venv`, Python 3.13, cobra 0.30.0, pytest 9.0.2, `--continue-on-collection-errors`. Both runs use short-path
copies:
- **Before:** `git archive af7c5432` in `/private/tmp/claude-501/mc7_before`.
- **After:** the same archive plus every modified or untracked file
  ([logs/pytest_after_extra_files.txt](logs/pytest_after_extra_files.txt)), in `/private/tmp/claude-501/mc7_after`.

The new test alone, in the worktree, passed: 8 tests, 43 subtests, 48.7 s.

| Run | Passed | Failed | Errors | Skipped | Subtests passed |
| --- | --- | --- | --- | --- | --- |
| Before (`af7c5432`) | 493 | 3 | 18 | 24 | 1,348 |
| After | 501 | 3 | 18 | 24 | 1,391 |

[logs/pytest_compare.txt](logs/pytest_compare.txt) compares the junit files test by test. The only differences are the
8 new `test_methylcitrate_cycle_7ccbb568` tests, all passed (+43 subtests). The standing failures and errors are the
same before and after:
- 3 failures: `test_dfba_new_fn_essentiality`, `test_lipid_moiety_ledger`, `test_vacuole_supply`;
- 18 errors: 14 in `test_lp_sn12_candidate` and 4 `gem_review_lit` collection errors.

The environment, the totals and the junit SHAs are in [logs/pytest_env.txt](logs/pytest_env.txt).

The after-copy was made at 19:22, from the files listed then. That list includes its own temporary `.z` name, which
is harmless. Every code file in the copy is identical to the worktree; the audit checked this.

## 7. Figures

Hand-laid schematic, same layout as the investigation's. Script: [scripts/draw_schematic.py](scripts/draw_schematic.py).
It reads only this round's solve outputs from the two built files.

| File | What it shows |
| --- | --- |
| [figure/methylcitrate_status.png](figure/methylcitrate_status.png) | Largest flux each drawn reaction can carry, before / after. Green: freed inside the matrix (R_OBDHm, R25, R95, R552). Red: freed only through side route B (R834). Dashed: closed (R1172, R2075). |
| [figure/methylcitrate_after_route.png](figure/methylcitrate_after_route.png) | After, with the cycle forced (R25 = 1): R963 → R_OBDHm → R25 → R95 → R552; no transport of CoA or acyl-CoA |
| [figure/methylcitrate_before_wt.png](figure/methylcitrate_before_wt.png), [figure/methylcitrate_after_wt.png](figure/methylcitrate_after_wt.png) | pFBA at the growth optimum. The drawn fluxes are the same before and after; the cycle is 0 in both |

## 8. Genes

Names are not verified in this round. The roles are the model's rules: a model or GPR assignment is not experimental
validation.

| Gene | Model role | Evidence |
| --- | --- | --- |
| `YALI1E31947g`, `YALI1F27556g`, `YALI1B13179g`, `YALI1D30618g` | subunits in R693's 5-gene AND rule (pyruvate dehydrogenase); now also R_OBDHm | GPR only; no name verified here |
| `YALI1D26360g` | in R693 and 12 other reactions: branched-chain 2-oxo acid dehydrogenase (R1968, R1969, R1970, R2281), 2-oxoglutarate dehydrogenase (R613, R614), R1976 (named glutaryl-CoA dehydrogenase, with the 2-oxoglutarate dehydrogenase rule) and glycine cleavage (R357, R360, R361, R364, R365); the subunit these rules share. Now also R_OBDHm | GPR only; no name verified here |
| `YALI1F03803g` (R95), `YALI1F39620g` (R552), `YALI1F20788g` (R834) | genes of the freed reactions, on the consensus list | experimental list; all three are FN at 15 % before and after (KO/WT 1.0) |
| `YALI1E01019g`, `YALI1E03300g` (R25), `YALI1C24124g` (R552, OR partner) | genes of the freed reactions, not on the consensus list | TN before and after (KO/WT 1.0) |

## 9. What this means, method notes, and proof levels

**What it means.** The methylcitrate cycle is no longer a dead branch in this line's model.
- R25, R95 and R552 carry flux inside the matrix. R95 and R552 have listed-essential genes.
- R834, also listed, carries flux only through side route B.
- No call, score or WT growth changes. Some capacities rise slightly (§4).
- The cycle is reachable, not required. On Task A medium the optimum does not use it, so the change cannot be
  checked against these knockouts.
- A propionate-medium test would be the real check, but the line has no such condition.

**Could the method be better?**
- **Per-compartment conserved-moiety scan.** The blockage came from a conserved-moiety constraint that no single
  reaction shows. A scan for moieties confined to one compartment, built from the left null space restricted to that
  compartment's species, would find such traps before a closure is built. Test it by checking that it flags the
  matrix CoA moiety on 7ccbb568 and not on the model before R1172 was closed.
- **Step tests.** The other line's step test depends on that line's build chain. A test that reads a tracked
  candidate, as here, would travel between lines.

**Proof levels:**
- **Provable:** the CoA-moiety identity (static); the diff is limited to the targets (whole-model and text diff); the
  balance of R_OBDHm.
- **Observed in runs:** every LP value above, and the sweep. All were optimal or certified.
- **Tested:** the step's guards, by the focused test (38 rejection subtests and 5 PDH-knockout subtests, plus
  preservation, roundtrip and context checks). The audit added four guard probes, all of which failed closed.
- **Untested:**
  - the rule of R_OBDHm;
  - the absence of a mutase;
  - growth on propionate;
  - whether acyl-CoAs cross the inner membrane on R2025, R2113 or R2077;
  - the biology of side route B.

## 10. Limits and next gates

- **Inferred rule.** R_OBDHm's rule is INFERRED (user decision). No *Y. lipolytica* experiment places this step, its
  catalyst or its compartment. The branched-chain 2-oxo acid dehydrogenase complex is the other candidate and is not
  used.
- **R2075's closure** rests on absence by annotation and homology. No experiment shows absence.
- **Side route B** stays open. Whether R834 should be active only through it is a separate decision. The gene-less
  acyl-CoA transports were reviewed in a separate session, "Review gene-less acyl-CoA transports across the inner
  membrane", on branch `claude/lucid-visvesvaraya-a40bbb` (uncommitted).
- **Copied data.** The data and module are byte-identical copies from an uncommitted round on another line. If that
  round changes them before it lands, this copy must be re-synced. When the lines merge, `--methylcitrate-cycle` will
  be defined on both.
- **Baseline sweep.** The comparison sweep ran on 353cefdc, which is LP-identical to 7ccbb568 by a static check, not
  by a re-run on 7ccbb568.
- **Not done:**
  - a propionate condition;
  - GLPK on the built file (the investigation's GLPK check was in memory);
  - FVA at 100 % of WT growth.
- **Next gates:**
  - the user's review of this candidate;
  - a decision on side route B;
  - a commit or push, which needs a separate instruction.

## 11. Compute and commands

LPs: 54 targeted + 53 figure + 1,123 sweep = **1,230** (budget 1,300). Three builds, no LP. One full pytest suite
before and one after; the LPs that the existing tests run inside those suites are not counted here.

```bash
# from model/reports/methylcitrate_cycle_7ccbb568_20261116/
/opt/miniconda3/bin/python3 scripts/extract_args.py ../../candidates/E5_353cefdc_xref_20261115.build.json scripts/args_7ccbb568.txt
/opt/miniconda3/bin/python3 scripts/before_state.py ../../candidates/E5_353cefdc_xref_20261115.xml before.json
zsh scripts/run_builds.zsh
/opt/miniconda3/bin/python3 scripts/model_diff.py ../../candidates/E5_353cefdc_xref_20261115.xml ../../candidates/E5_7ccbb568_methylcitrate_20261116.xml model_diff.json
# from platform/, PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts
/opt/miniconda3/bin/python $TASK/scripts/targeted_solves.py MODEL.xml MEDIUM.csv PROFILE.json OUT.json   # control, candidate
/opt/miniconda3/bin/python $TASK/scripts/figure_flux.py CONTROL.xml CANDIDATE.xml MEDIUM.csv PROFILE.json $TASK/figure_flux.json
/opt/miniconda3/bin/python $TASK/scripts/figure_capacity.py CONTROL.xml CANDIDATE.xml MEDIUM.csv PROFILE.json $TASK/figure_capacity.json
zsh $TASK/scripts/run_sweep.zsh ad901ddb4a2cfe5b09eb29f9c4b6630c0ee77cf93cef5be5aff3c7dfc954c279
/opt/miniconda3/bin/python3 $TASK/scripts/compare_sweeps.py $RR/artifacts/tasks/folinic_acid_module_20261112/sweep684_candidate.tsv $A/sweep_candidate.tsv $TASK/sweep_compare_vs_7ccbb568.json
/opt/miniconda3/bin/python3 $TASK/scripts/draw_schematic.py $TASK/figure
```
