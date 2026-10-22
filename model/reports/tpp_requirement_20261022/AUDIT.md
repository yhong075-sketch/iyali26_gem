# AUDIT: TPP requirement (TPP-REQUIREMENT-20261022)

**Date:** 2026-10-22. **Auditor:** an independent, read-only subagent with no conversation history.
It edited no repository file, ran no builder and made no commit. This is its report as returned,
lightly reformatted. The responses are in [REPORT.md §8](REPORT.md).

The audit covered the first candidate, `12efa141…`. After the audit, the guard and the note texts
were fixed and the candidate was rebuilt (REPORT §8).

**What it read:**
- `TASK.md`, `AGENTS.md`
- `model/curation/thiamine/tpp_requirement.json`
- `platform/scripts/gem_annotate/tpp_requirement.py`, plus the `cli.py` and `main.py` diffs. The stored
  `logs/code_diff.patch` equals the working-tree diff once index lines are ignored.
- `platform/tests/test_tpp_requirement.py`
- Both models and every `.build.json`, including the three under `.scratch/tpp_requirement_20261022/`
- All report outputs, logs and scripts
- The W29 gene map, `consensus_essential_genes.csv` and `orphan_gene_search_20261019/db_provenance.json`

**What it ran:**
- A raw XML diff of the two models.
- SHA-256 of all models.
- Regenerated `before.json`, `after.json` and `tpp_requirement.json` into scratch. All three are
  byte-identical to the stored files.
- The focused pytest: 5 passed, 18 subtests.
- About 20 diagnostic LPs (GLPK, plus Gurobi for cross-checks).
- A full re-run of 2 × 1,074 single-gene knockouts.
- Guard probes on in-memory copies.
- KEGG REST GETs.
- A streamed read of the yeast-GEM YAML. It was hashed and parsed but not saved.

`REPORT.md` appeared during the audit. It was not audited.

## Claims

| # | Verdict | Evidence |
|---|---|---|
| C1 | **supported** | `diff` of the two XML files gives exactly 3 hunks: (1) R1342 `fbc:upperFluxBound` goes from `cobra_default_ub` to `cobra_0_bound`, plus 7 `tpp_requirement_*` notes; (2) 7 `tpp_requirement_*` notes on biomass_C; (3) `<speciesReference species="M_m975__91__C_cy__93__" stoichiometry="1.2e-06">` inside biomass_C's listOfReactants. Nothing else differs: no other reaction, species, gene product, group, parameter, objective or model-level note. Counts are 2311/1876/1074 in both. The objective is identical. m975 reactions gain only biomass_C. |
| C2 | **partially supported** | Recomputed hashes: `control_no_tpp_flag.xml` and `control_no_tpp_flag_after_fix.xml` = cf0825b8…; the candidate and `tpp_seed1.xml` = 12efa141…. The candidate, seed1 and after-fix control `.build.json` files record identical `source_sha256` and `data_sha256` (tpp_requirement.py 6847ef77, JSON 7b1fa25c, HEAD abe7e40). **Unverified:** the seed values. PYTHONHASHSEED=0/1 appears only in `commands.txt`; no build log or build.json records the environment. The evidence also depends on `.scratch/` (git-ignored), which disappears with the worktree. |
| C3 | **partially supported** | Holds as claimed: all validation runs before any edit (lines 33–81); the focused test's 13 fail-closed cases and 5 scope cases pass; the new model file gives `already_correct`; the baseline file gives `applied`; the in-build unescaped notes are accepted. **Probes it wrongly accepts:** (i) R933 made reversible plus any mitochondrial TPP source is `applied`; (ii) R753 made reversible and gene-less plus a TTP source is `applied`; (iii) a spec with R1342 `bounds_after` [0,500] is `applied` and leaves R1342 at (0,500); (iv) a spec whose biomass `after_stoichiometry` also changes m1000 is `applied` with only TPP added, and re-applying then fails "precondition differs"; (v) a non-biomass objective is `applied`. **It can leave partial edits:** a spec with R1342 `bounds_after` [10,0] raises ValueError *after* biomass_C has already gained −1.2e-6 TPP (with no notes). |
| C4 | **supported** | `cli.py` uses `store_true`, so the default is False (test `parse_args([])`). `main.py:126` raises "defined on the GPI co-products candidate" when `--gpi-coproducts` is missing. `main.py:760-764` is the block immediately after the `gpi_coproducts` block. |
| C5 | **supported** | Re-queried KEGG REST: `link/yli/K00946` and `link/sce/K00946` return HTTP 200 with a 1-byte empty body, as do `link/yli/ec:2.7.4.16` and `link/yli/R00617`. `link/yli/K00949` returns `yli:2911892`. `get/yli:2911892` gives SYMBOL YALI2_E01340g and NCBI-ProteinID XP_504226. Gene map row: `XP_504226.1,YALI1_E25443g,YALI2_E01340g,…,YALI1E25443g,True,R752;R753`; the assembly is GCF_001761485.1. yeast-GEM streamed from raw.githubusercontent: 3,301,690 bytes, SHA-256 2ddfaaa7…dc17 (matches), version "9.1.1". `r_4598` "cofactor pseudoreaction" has `s_1475: -1.20000004244503e-06`. s_1475 is TDP, C12H16N4O7P2S, charge −2, kegg C00068, bigg thmpp. `r_4041` consumes `s_4205` (cofactor) at 1, so the value is per gDW. KEGG `get/cpd:C00068` = Thiamine diphosphate. |
| C6 | **supported** (stored outputs; BLAST not re-run) | `sprot_titles.tsv` has 27 titles starting "Thiamine-monophosphate kinase". Every `GN=thiL` entry is in that set. Log: 27 queries, 152 forward hits, 0 W29 proteins at E≤1e-5. CSV: min E = 0.006 (O27447 → YALI1F27904g, qcov 0.269). Control: P35202 → XP_504226.1 (YALI1E25443g) at E 1.78e-24. Selection is by title, because the title table carries no EC numbers. |
| C7 | **supported** | GLPK re-run: WT 0.31242674464567 (base) and 0.31242674050093 (new). New WT pFBA: R934 = R752 = 3.749e-7. THI80 KO −1.7e-14 (stored 1.7e-30); Gurobi KO 0.0. R1340 closed gives 0.3124265729 via R970 → R413/R645 + R1344 → R754 → R1337 → R752. KO with R1342 reopened gives 0.3124267372 (KO/WT 0.99999999). |
| C8 | **supported** | Full independent re-run of both models × 1,074 genes. Max deviation from the stored TSV is 0 for both models, with no status mismatch. The only class change is YALI1E25443g (non-essential → essential) at 1/5/10/15%. Next-largest growth changes: YALI1D25037g (R934 transporter) 1.7e-7, then 5.0e-9. The 20 infeasible genes are the same set on both models. Caveat: these 20 are unresolved on both, so the 1e-6 statement covers 1,054 genes. |
| C9 | **supported** | Direct LP max of R752: GLPK and Gurobi both 3.749120886e-7 (= 1.2e-6 × growth). Stored FVA max was 2.907e-7. The auditor's FVA re-runs gave 3.745e-7 (GLPK) and 3.700e-7 (Gurobi): all below the true maximum and not reproducible between runs. |
| C10 | **supported** | `check_mass_balance` for R1342 gives {C −6, H −5, N −1, O −4, P −1, Cl −3} on both models, with all charges 0. biomass_C raises "No elements found in metabolite m401[C_cy]" (chitin `H2O(C8H13NO5)n`). The new term m975 is itself formula-less, so biomass stays non-computable even if chitin were fixed. |
| C11 | **supported** | No `<>&'"` in any of the 4 note texts × 2 targets or in the authorization. Notes read back from the new XML equal the spec verbatim for both reactions. |
| C12 | **supported** (with a provenance caveat) | `pytest_after_results.txt`: "11 failed, 226 passed, 31 skipped, 17 errors, 554 subtests". Before: "11 failed, 221 passed, …, 536 subtests". The per-test diff shows only the 5 new tests as PASSED (+18 subtests). Caveat: the "before" run lasted about 10:34:00–10:40:00, but `cli.py` and `main.py` have mtime 10:37:42. The edits landed during the run (module imports at collection likely predate them). |

## Findings

### (a) Scientific overstatements

1. **`tpp_requirement.json:163`, also written into biomass_C's model notes.** The rationale says "No reaction consumed cytosolic TPP (m975)", and TASK.md says the same twice. Both are literally false. R753 (ATP + TPP → ADP + TTP) and R933 (TPP cy → mi) consume m975[C_cy]. They carry no flux only because m976 and m1116 have no other reactions. Correct wording: "no reaction drains cytosolic TPP to a sink."
2. **TASK.md, "No TMP kinase in either yeast (verified 2026-10-22)".** What was verified is narrower: there is no KEGG K00946 assignment in yli (the CLIB122 genome) or sce, and no ThiL BLAST hit in W29 at E≤1e-5. That is not the absence of TMP kinase activity. The data file's tier text states this limit correctly; the TASK.md heading drops it.
3. **json evidence tier, "absence of a ThiL orthologue verified by … BLAST".** The positive control, THI80 P35202, is a TPK and belongs to a different family from ThiL (PurM-like). It shows that the database and pipeline work. It does not show that a divergent eukaryotic ThiL would be detected. Low severity.
4. **json targets, "TPP (m975)".** m975 has no formula or cross-reference, and its name "TDP_" is ambiguous. Its identity as thiamine diphosphate is inferred from context: R752 = KEGG R00619, R753 = R00616, R933 is named "thiamine diphosphate transport". The data file does not label this identification as inferred.

No overstatement found on: the yeast-GEM coefficient, the KEGG queries, the gene ID mapping, the consensus-list statement, or the gene naming.

### (b) Outside the authorized scope or process

1. **Stop condition not observed.** TASK.md says "Stop and report on: a guard rejection". The seed-0 and seed-1 builds were rejected ("biomass_C notes outside this curation differ"). Within about 1m45s the guard was changed (html.unescape tolerance added), the curation JSON was regenerated, and the rebuild started. No file records that the user was told before the rebuild. The fix appears technically correct, but it loosens a gate (AGENTS.md §7).
2. **Compute budget overrun.** TASK.md sets "about 30 diagnostic solves". After authorization, `solves.json` records `lp_count` 83 (60 of them FVA) and `r752_max_check` adds 6 LPs (3 with Gurobi, a solver not named in the budget). `reproduce_figure` (53 LPs) ran before the change. The knockout count (2,150) is within budget.
3. **Large outputs in the worktree.** Three 14 MB control and seed builds are in git-ignored `.scratch/`, not under `$IYALI26_RESEARCH_ROOT/artifacts/tasks/tpp_requirement_20261022/` as AGENTS.md §3 requires. `logs/sha256.txt` points at `../.scratch/…`, so the C2 evidence is lost when the worktree is cleaned.
4. **Within scope.** Only the two authorized targets changed. R1342 is closed, not deleted. HEAD is still abe7e40. No `model/STATE.md` entry exists yet.

### (c) Code defects (`platform/scripts/gem_annotate/tpp_requirement.py`)

1. **Partial edit on a malformed spec.** biomass_C is edited before R1342's bounds are set. With R1342 `bounds_after` [10, 0], setting the bounds raises after biomass_C already carries TPP at −1.2e-6 with no notes. In the builder the exception aborts before any write, but a caller that catches the error gets a half-edited model.
2. **Spec values are not validated.** The scope check validates keys and status only. `bounds_after` [0, 500] is "applied" and R1342 stays open. Biomass `bounds_after` and non-TPP `after_stoichiometry` entries are silently ignored, and the next re-apply fails "precondition differs". Fix: assert R1342 `bounds_after == [0,0]`, biomass `after == before + {TPP: coeff}`, and `bounds_after == bounds_before`.
3. **The TPP-neighbour check covers only reaction IDs.** The bounds, GPR and stoichiometry of R753 and R933 are not checked, nor the m1116/m976 neighbourhoods. If an upstream step makes R933 reversible and adds a mitochondrial TPP source, the guard returns `applied`, and the THI80 knockout would be bypassed through mitochondria. Fix: pin R753 and R933 like R752, and pin the reaction lists of m976[C_cy] and m1116[C_mi].
4. **Objective not checked.** If biomass_C is not the objective, the change is applied and has no effect. Minor.
5. **Escape handling is asymmetric.** `_unescaped` uses `html.unescape` (the HTML5 entity set: it accepts `-&#62;` and legacy entities without a semicolon; harmless). The conflict check compares raw values. If a future curation text contains `'` or `"`, a file-read model holds `&apos;`/`&quot;` and the conflict check raises on a correct model. Today only `make_curation_json.py`'s assert prevents this.
6. **Test gaps.** No case covers defects 1–4. Step ordering is verified only by reading the code.

### Limit (inferred, not a defect)

The TPP coefficient (1.2e-6) is only about 12× the GLPK and Gurobi feasibility tolerance (1e-7). C9 shows that tolerances already distort FVA at this flux scale. A solver exploiting the TPP-row tolerance could, in principle, return a THI80-KO growth up to about 1e-7/1.2e-6 ≈ 0.083 (KO/WT ≈ 0.27), which would flip the 15% call. Neither GLPK nor Gurobi did here (KO = −1.7e-14 and 0.0). Suggested check: re-solve the KO with feasibility tolerance 1e-9, or report the m975 row residual in the KO solution.
