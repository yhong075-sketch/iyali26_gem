# R2135 acyl-CoA pool remap and R2145 1:1 rescale, candidate e9d35857 (2026-10-27)

## 1. Result

- **Candidate:**
  [`model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_20261027.xml`](../../candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_20261027.xml),
  SHA-256 `e9d358574eeb4369ecd39f63e5a2ad6c113c4681232c6504cb1ef42c4fe89b1b` = `00e8e397` +
  `--r2135-acyl-coa-pool`. Not a designated reference. Not committed.
- **Authorization:** "Approve Option A" (user, 2026-10-27), meaning Option A of the Stage 1 draft, whose
  section 12 gives the full sentence (TASK.md section 0).
- **Two edits, both inherited unchanged from iYli21:**
  - **R2135** "AC-CoA transport": acetyl-CoA `m68[C_cy]` ⇌ acyl-CoA `m1984[C_er]` becomes ER-membrane
    pooled acyl-CoA `m567[C_em]` ⇌ `m1984[C_er]`.
  - **R2145** (PA phosphatase, GPR `YALI1D35593g`): `0.01 m1993 + H2O → 0.01 m1994 + Pi` becomes
    `m1993 + H2O → m1994 + Pi`.
- **Effect:**
  - **Leak removed.** The model no longer makes long-chain acyl-CoA from acetyl-CoA, which had let the
    opened R1726 loop replace fatty-acid synthase (`../kennedy_closures_20261027/STOP.md`).
  - **Phosphate source removed.** R2145 no longer releases phosphate from water.
  - **Growth:** WT on Task A + PO1f goes 0.3120049 → 0.3116824 (−0.10 %).
  - **Knockouts:** no gene changes class at any cutoff, with or without R1726 opened on top.

## 2. Why (evidence; tiers as in TASK.md section 3)

- **The species are mismatched.** `m1984` is MetaNetX MNXM1136895 "acyl-CoA", the same id and Reactome id
  as the pooled acyl-CoAs `m567` / `m570` / `m790`. Its only consumers, R2134 and R2168, are KEGG R01496 /
  R02251, which use acyl-CoA C00040 (C22H35N7O17P3SR). `m68` is acetyl-CoA C00024 (verified).
- **Conclusion:** linking acetyl-CoA to `m1984` is a species-identity error (inferred: the local entry points
  to a different entity). The name "AC-CoA" most likely meant acyl-CoA (inferred).
- **R2145's 0.01 coefficients** match ModelSEED rxn09551's lumped "yeast-specific" PA / DAG. The record is
  verified; that iYli21 copied R2145 from it is inferred from R2145's `seed.reaction` annotation. ModelSEED
  lists rxn09551 as an alias of BiGG DAGPYP_SC (iMM904 / iND750) and flags it mass-imbalanced itself (found
  by the audit).
  In this model the same species are linked 1:1 to per-molecule PA and DAG (R2144, R2169), so every
  dephosphorylated unit released 100× its phosphate. The new form is the per-molecule EC 3.1.3.4 reaction
  catalysed by the PAH1 ortholog `YALI1D35593g` (KEGG K15728; orthology, verified).
- **How the baseline used the leak** (verified, Gurobi 1e-9, built files, [built_checks.json](built_checks.json)):
  - TAG for biomass was made through PA → R2144 → R2145 → DAG → R2169 → R2168, with acetyl-CoA from R2135
    as the third chain.
  - R2135, R2168 and R2145 had FVA minima above 0 at 99.9 %.
  - R2145 at 0.0593 supplied about 0.059 of the cell's phosphate. Phosphate uptake was 0.0259; it is 0.0846
    on the candidate.

## 3. Implementation

- **Data:** [`model/curation/chemistry/r2135_acyl_coa_pool.json`](../../curation/chemistry/r2135_acyl_coa_pool.json).
  It holds both reactions' before and after equations, the species guards (formula, charge, compartment,
  MetaNetX ids, `m1984`'s reaction set), the notes and the evidence rows.
- **Step:** `platform/scripts/gem_annotate/r2135_acyl_coa_pool.py`, `apply_r2135_acyl_coa_pool`.
  - **Scope pinned in code:** the after-equations, plus a check that the species list covers every id.
  - **Fails closed** on:
    - missing reactions or species;
    - any species, name, GPR or bounds change;
    - an unexpected equation;
    - conflicting notes;
    - this curation's notes on the old equation;
    - the new equation without the notes;
    - a half-applied pair.
  - **Returns** `applied` / `already_correct`, and touches only the two reactions' stoichiometry and notes.
- **Flag:** `--r2135-acyl-coa-pool` (`cli.py`, `main.py`). It requires `--cofactor-biomass-terms`, runs after
  it, and runs before `--r1726-cpt-open`.
- **Test:** `platform/tests/test_r2135_acyl_coa_pool.py`, 4 tests and 29 subtests, each guard matched by its
  message.
- **Code review:** a fresh general-purpose reviewer following `.claude/agents/code-reviewer.md` found no
  bugs, and one medium-low fail-closed gap. A data file with an unknown after-state id raised `KeyError` after
  R2135 had already been rewritten; a swapped after-state was accepted. It also raised three low items: notes
  on the old equation were accepted, guard tests were missing, and `m1984`'s reaction context was not
  pinned. All were fixed (after-equations pinned in code, species coverage, notes guard, context guard, 11
  more subtests).
- **Test bug found while fixing the review:** `Model.copy()` shares each metabolite's annotation dict with
  the original. The identity subtest mutated the fixture and went undetected only because the subtests after
  it fail earlier. The test now replaces the dict.
- **Not changed:** the reviewer's out-of-scope note that `--r1726-cpt-open` should require this flag. That is
  a change to the paused R1726 work and needs the user's approval.
- **pytest:** see section 4.

## 4. Verification

| Check | Result |
| --- | --- |
| pytest | 274 → 278 passed. The difference is exactly the 4 new tests; the same 11 failures and 18 errors (missing `httpx2`, lipid / vacuole tests) as before ([pytest_compare.txt](logs/pytest_compare.txt)). "Before" is the Kennedy task's last full run on the same code state. |
| Builds | Control without the flag = `00e8e397…` byte for byte. `PYTHONHASHSEED=0` and `=1` (set on the command line; the build record does not store the variable, so only the file names carry it) with the flag = `e9d35857…` (round 1, before the review fixes). The final candidate = `e9d35857…` again, so the fixes changed no model byte ([sha256.txt](logs/sha256.txt), [sha256_round1.txt](logs/sha256_round1.txt)). `.build.json`: HEAD `31e7aa4`, dirty (this task's and the Kennedy task's uncommitted files), input and source unchanged, 2,316 / 1,882 / 1,081. |
| Whole-model diff | 4 entries: R2135 stoichiometry, R2135 notes (+6 keys), R2145 stoichiometry, R2145 notes (+6 keys); no existing note changed ([model_diff_vs_00e8e397.json](model_diff_vs_00e8e397.json)) |
| WT | 0.3120049202 → 0.3116823733, optimal, checked (residual 1.8e-14) |
| No-uptake LPs | Max R_NTP1 and max biomass_C are 0 and optimal on both models, with `xMAINTENANCE` 7.8625 relaxed |
| Route at the optimum | R2135 / R2168 / R2145 carry 0 (FVA at 99.9 % [0, 6.0e-3]: still usable, not chosen); R2134 still [0, 0]; Pi uptake 0.0259 → 0.0846 |
| Mass balance | Not checkable: `m567`, `m1984`, `m1993` and `m1994` have no formulas. cobra reports R2135 `{}` and R2145 `{H: 1, O: 3, P: 1}`, the same as before, from the formula-less lipids. Charge is 0 everywhere. A limit, not a pass. |
| Knockout sweep | ([ko_compare_summary.json](ko_compare_summary.json), [ko_compare.tsv](ko_compare.tsv))<br>- 2 × 1,082 LPs, Gurobi 1e-9, reset per knockout, knocked-flux check.<br>- No failing point; largest residual 1.5e-11, largest bound violation 1.8e-10.<br>- **No class change at 1 / 5 / 10 / 15 %**: TP / FN / FP / TN at 15 % stay 111 / 204 / 42 / 692 (MCC 0.383).<br>- The same 20 infeasible knockouts (16 after the mtDNA exclusion) remain unresolved; the 16 are excluded from the confusion counts.<br>- Positive coverage: 326 consensus genes are in the model; 111 + 204 + 11 unresolved positives.<br>- Development-reference scores; user-defined negatives (model genes not on the list); mtDNA out of scope. |
| Stacked R1726 check | ([ko_sweep_r1726_stacked_summary.json](ko_sweep_r1726_stacked_summary.json))<br>- R1726 opened in memory on the candidate: WT 0.3116824, unchanged, so R1726 adds no growth.<br>- 1,082 knockouts, no failing point, **no class change versus the candidate** at any cutoff.<br>- FAS `YALI1B19844g` / `YALI1B25427g`, desaturase `YALI1C07638g` and ligase `YALI1D22124g` stay lethal. |
| Re-census (GLPK) | ([census_compare_open.json](census_compare_open.json), [census_compare_taskA_po1f.json](census_compare_taskA_po1f.json))<br>- Open condition: 722 → 720 blocked. R1343 and R1344 are freed: the GLPK fringe pair the cofactor task found active at Gurobi 1e-9 on both models, not an effect of this change.<br>- Task A: 852 → 852.<br>- No newly blocked reaction; gap metabolites unchanged.<br>- The audit's Gurobi 1e-9 FVA on both models gives R1343, R1344 and R970 a maximum of 1.47e-4; R970 stays GLPK-blocked on both. |
| Budget | TASK §10 allowed 3 builds and about 20 built-file LPs. **Used: 6 builds** (control and both seeds, repeated after the review fixes) **and about 70 built-file LPs** (both models with and without R1726). Both sweeps ran within budget. |

## 5. What it means

- **The baseline's growth was 0.10 % too high.** Part of its TAG was built with acetyl-CoA. On the candidate,
  the pFBA point makes TAG through the lipid-particle DGAT R1771 with real pooled acyl-CoA. A plain FBA optimum
  can use R2135 → R2168 instead (7.3e-5, found by the audit): the optima are not unique, and both routes now
  use real acyl-CoA.
- **Phosphate:** the baseline drew about 70 % of its phosphate from R2145 (pFBA point). Phosphate is not
  limiting on Task A, so this did not change growth, but it would matter on a phosphate-limited medium
  (inferred, untested).
- **Kennedy R1726:** with this fix, opening R1726 no longer changes any class and adds no growth on Task A.
  The Kennedy task can be resumed on top of this candidate (§7).

## 6. Limits

- **The remap target is a judgement.** `m567[C_em]` is the ER-membrane pooled acyl-CoA. A pool-to-pool link
  between the ER membrane and the ER is a species-identity correction, not a measured transport (inferred).
  The ER-membrane pool rather than the lipid-particle pool `m570` was chosen because `m1984` is in `C_er`.
- **Mass balance cannot be checked** for any of the four lipid species.
- **The remaining 0.01-scaled reactions are unchanged:** the IPC and MIPC synthases, R2171 and R2187
  (TASK.md section 11). The IPC synthases carry no flux on Task A; the MIPC synthases and R2187 were not
  checked.
- **Coverage:** only the Task A medium + PO1f. Gene identities are orthology (KEGG), not experiment. The
  knockout classes are development-reference scores.
- **Not done:** no Escher maps; no other media; no memote.

## 7. Next gates

1. **Kennedy R1726 (paused).** It can be re-run on `e9d35857` with `--r2135-acyl-coa-pool
   --r1726-cpt-open`. Expected, per the stacked check: no class change and no growth gain; UM12 unblocked
   (to be re-censused). Making `--r1726-cpt-open` require `--r2135-acyl-coa-pool` is part of that change.
   It needs the user's approval.
2. The other 0.01-scaled reactions of R2143–R2187, and the two TAG and two DAG representations of the
   block.
3. The DAG kinase closure R1720.

## 8. Audit

See [AUDIT.md](AUDIT.md).

## 9. Large outputs

`$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2135_acyl_coa_pool_20261027/` holds:
- builds: control, seed 1, round 1;
- census runs;
- the stacked sweep table `ko_sweep_r1726_stacked.tsv`;
- the Stage 1 variant log.
