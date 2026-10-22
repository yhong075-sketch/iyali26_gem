# Independent audit: heme a + iron(2+) in biomass_C (2026-10-22)

A fresh read-only reviewer, with no history of this work, was given TASK.md, the diff, the logs, the
deliverables and claims C1–C11. It made no repository changes and used 18 LPs (limit 30) and 5
read-only web requests. Its scratch scripts are in this session's scratchpad (not kept). Its report
is reproduced below, lightly reformatted. The main session's responses are in
[REPORT.md §8](REPORT.md#8-independent-audit).

## Verdicts

| Claim | Verdict | Evidence checked |
|---|---|---|
| **C1** Baseline SHA | **Supported** | `shasum` gives cf0825b8…686b7c. It matches the baseline `.build.json` output.sha256 and `model/candidates/README.md`. The file is committed at d43e8b9 and unmodified. |
| **C2** Control build reproduces the baseline | **Supported** | `.scratch/.../control_no_flag.xml` has SHA cf0825b8, and `cmp` against the baseline shows no difference. Its argv has no `--heme-iron-biomass`. |
| **C3** New SHA, seed identity, argv | **Supported** | New file is 7f836c3a…16539d; `heme_iron_seed1.xml` is byte-identical (`cmp`). The argv equals the baseline argv plus `--heme-iron-biomass` and the output path; argv[0] differs only by worktree path. Every source and data hash recorded in the new `.build.json` matches the files on disk. **Limit:** the seed values appear only in commands.txt, not in `.build.json`. |
| **C4** Diff is only the two terms and seven notes | **Supported** | A raw text `diff` of the two XMLs shows exactly 9 added lines and nothing removed or changed: 7 `heme_iron_biomass_*` note lines, plus `m1635[C_cy]` 1e-06 and `m948[C_cy]` 3.04e-05 in `R_biomass_C`'s list of reactants. The reaction, species and gene counts are equal. |
| **C5** Guard behaviour | **Partially supported** | Only `biomass_C` is edited; every raise comes before any edit; the step is idempotent; the flag requires `--coq9-biomass-dilution` (main.py:126), which itself requires E5 (main.py:114). **Gaps, each shown by an in-memory probe on cf0825b8** — the step still returns "applied" when: R1189 is deleted; R863 or R302 is closed; R383 is deleted; the objective is moved off `biomass_C`; an existing `biomass_C` note value is changed; R1710 is opened. The data file's `baseline_model_sha256` is recorded but never enforced, and coefficient sign and magnitude are not validated. |
| **C6** pytest unchanged plus 4 new passes | **Supported** | Outcome lines extracted from both logs: identical except 4 new PASSED in `test_heme_iron_biomass`. FAILED/ERROR messages identical (11 failed, 17 errors, 31 skipped in both). Subtests 536 → 553 (= 12 + 5). Warnings 7 → 8, the extra one from the new test. |
| **C7** Coefficients and sources | **Supported** | Archive `SHA256SUMS` checks OK. Git blob SHA-1s of the archived TSV and YAML equal the blobs GitHub lists at commit 2d594ae1; the YAML says version 9.1.1. `r_4598`: s_3714 heme a [c] at −1e-6; `r_4599`: s_0924 iron(2+) [c] at −3.04e-5 (stored at float32 precision). Ferroheme b (s_0712) is in none of the 10 pseudoreactions. The TSV cites "Forster2003" and DOI 10.22161/ijeab/2.2.2. In `model/source/iyali26.xml`, R1710 has bounds (0,0) and `m1635[C_cy]` at −1e-6. The primary sources were not opened by anyone. |
| **C8** Growth and closures | **Supported** | Fresh GLPK contexts: WT new 0.3124262956, baseline 0.3124267446. pFBA: R1189 −9.8102e-6, 60.9 % of the cap. Closing R1189, R2193, R302, R383 or R1709 each gives optimal with \|growth\| ≤ 7.8e-18. |
| **C9** Class changes and scores | **Supported** | Recounted from `ko_compare.tsv` with both corrections (reproduced): 13 changes, identical at 1/5/10/15 %, all non-essential → essential, 9 positives and 4 negatives, each gene's only reaction on the heme route or R2193. No other gene changes. At 15 %, corrected new 92/219/37/694/16 and baseline 83/228/33/698/16. |
| **C10** HEM4 artifact, HEM3 time limit | **Supported for HEM4/HEM3; "other genes" only partly settled** | Reproduced: WT solve, then the YALI1F13039g knockout in the same context gives "optimal" at WT growth; R792's forward variable has ub = 0 but primal 3.124e-7, the only bound violation above 1e-12. A fresh knockout gives 0.0; a fresh YALI1F34051g knockout also gives 0.0. That the violation was accepted because of GLPK's scaled tolerance is inferred. **Other genes:** reactions with 1e-9 < \|v\| < 1e-5 in WT/pFBA optima of both models: only the 11 heme/iron reactions and YALI1D22124g's acyl-CoA ligases (R1445/R1451/R1452), a gene essential in both sweeps. This is a screen, not proof, because sweep solves warm-start from the previous knockout. **Check that would rule it out:** verify every variable respects its bounds after every knockout solve, or re-solve every non-essential call fresh or with Gurobi at FeasibilityTol 1e-9, in both models. |
| **C11** Recorded findings | **Supported** | R2193's only gene YALI1F07747g: ncbigene 2908631 / yli:2908631; KEGG NAME "Phosphorus acquisition-controlling protein", Pfam HLH / bHLH_IRO3. Iron species have formula None and charge 0. R302 {charge +2, H −2, Fe +1}; R727, R384, R383 as stated. **New finding:** the gene's UniProt cross-reference A0A1D8NM31 comes back as "deleted" from UniProt REST, and KEGG's entry for 2908631 links to a different accession, Q6C2U3. |

## Issues

1. **Medium — numerical fragility caused by the change.** Heme-route fluxes of 1.6–3.1e-7 sit at
   GLPK's primal tolerance. In the change's own evaluation this already produced 2 of the 13 calls
   wrong or unresolved. Any later GLPK warm-start sweep, or any FVA or blocked-reaction tool that
   treats |v| < 1e-6 as zero (inferred), will misread this route. The acceptance criteria do not
   cover it.
2. **Medium — the budget statement contradicts itself.** "All within budget" while the first
   sweep's LP count is "not recorded"; possibly about 2,100 more LPs (inferred). The script was
   changed and the sweep rerun with no recorded check-in (AGENTS.md §7).
3. **Medium-low — guard gaps** (C5).
4. **Low — evidence not kept:** the stack sample, a script and output for TASK.md's reproduction
   table, raw KEGG responses for `gene_identities.tsv`, and the command that made the `.outcomes`
   files (which also contain `ERROR cobra.io.sbml` log lines).
5. **Low — unrecorded concurrency:** the after-change pytest run overlapped the three builds, the
   diff and the solves. Outcomes were unaffected.
6. **Low — `hem4_tolerance_check.py`:** the docstring gives Gurobi's default FeasibilityTol as 1e-6,
   but the recorded value is 1e-7. The residual metric measures |S v| only, not the bound violation
   that is the actual artifact.
7. **Low — compartment convention not stated:** cytosolic heme a makes the gene-less R1709
   mandatory. Biologically heme a stays in mitochondrial cytochrome oxidase (established/uncited).
   No gene result depends on this.
8. **Low — the score gain is not validation:** YALI1B19490g (adrenodoxin-NADP+ reductase) is hit
   only through heme a. In *S. cerevisiae*, Arh1/Yah1 essentiality is usually attributed to Fe–S
   cluster biogenesis (recalled, uncited), which the model does not contain.
9. **Low — TASK.md wording:**
   - "a fall of about 1e-6" — the actual fall is 4.5e-7;
   - "Heme b → heme o → heme a (R302, R384, R383)" — R302 makes heme b;
   - "the newest candidate".
10. **Info — REPORT §7 overstated:**
    - "A stale value can only show up as WT growth" — it can be an intermediate value;
    - "only the heme route gained new near-tolerance fluxes" was not checked in the report;
    - the section was labelled "(verified)" but contains inferred statements.
11. **Info — open items:**
    - the authorization quote cannot be checked against files;
    - the STATE.md entry was not yet written;
    - REPORT §8 was pending.
