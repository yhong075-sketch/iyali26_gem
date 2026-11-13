# Folinic-acid (5-formyl-THF) module: R748 to SHMT, R990 closed, cytosolic R_THFAT added, 2026-11-12

**Current candidate: `353cefdc` on main's `684f3752` (R2121 + R1172 closures); see §8.** The sections below document the first round on `e584151f`.

Scope, authorization and acceptance criteria: [TASK.md](TASK.md). Independent audit and its resolution:
[AUDIT.md](AUDIT.md). Not committed.

## 1. What changed and why

**Candidate.** `model/candidates/E5_coq9_resp_pkg_r132_..._mcnaugh_rxns_v3_folinic_acid_v2_20261112.xml`,
SHA-256 `85c623c86eeb4b818348a81f923f0ee4c0f9fe4b8bcee43b50cbd5e9554af834`. It is the e584151f build plus the new
opt-in flag `--folinic-acid-module`.

**Baseline.** `e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2`. Code at HEAD `154eb64c` plus this
round's uncommitted files.

**Superseded first build.** The first build (`..._folinic_acid_20261112.xml`, `d4832b29…`) differs from v2 only in
the text of two notes, which the audit asked to reword. It is the same optimisation problem; the byte-level diff is
four lines ([logs/xml_text_diff_v2_vs_v1.txt](logs/xml_text_diff_v2_vs_v1.txt)). All LPs below were run on
`d4832b29`. It is kept until the user decides.

**Build command.** Run from `platform/`. [scripts/run_build.zsh](scripts/run_build.zsh) sets `PYTHONHASHSEED`; seeds 0
and 1 give the same bytes, and the seed is in the log file names, not in the `.build.json`.

```bash
python -m scripts.gem_annotate $(cat ../model/reports/folinic_acid_module_20261112/scripts/args_e584151f.txt) \
  --folinic-acid-module --output-model ../model/candidates/<name>.xml
```

**Why.** The review (session notes, linked in TASK.md) found three problems:
- Cytosolic 5-formyl-THF `m159` had no producer. The exchange R991 is export-only, so the FAU1-type cyclo-ligase R64
  was blocked. That holds while R991's lower bound is 0, as in every medium under `model/conditions/`; it is proven
  from the bounds and the control FVA confirms it (§3).
- The mitochondrial 5-formyl-THF source R748 carried the glycine-cleavage T protein `YALI1F04082g`.
- R990, a gene-less "5-formyl-THF isomerase", has exactly the stoichiometry of R64 + R543.

The user called the module a gap to fix and authorized F1 + F2 + F4.

| Entity | Before (e584151f) | After |
| --- | --- | --- |
| R748 [C_mi] methenyl-THF + H2O → 5-fTHF + H+ | `YALI1F04082g`; EC 2.1.2.10; "tetrahydrofolate aminomethyltransferase" | `YALI1D28569g or YALI1E19557g` (as R359); EC 2.1.2.1; "5,10-methenyl-THF hydrolysis (SHMT side activity)"; `PROTEIN_CLASS` 2.1.2.1 plus 3 evidence notes |
| R990 [C_cy] ATP + 5-fTHF + H2O → ADP + Pi + 10-fTHF, gene-less | (0, 1000) | (0, 0) plus 3 notes; kept |
| R_THFAT [C_cy] (new) | — | m160 + H2O → m159 + H+, (0, 1000); `YALI1D28569g or YALI1E19557g` (as R358); group "One carbon pool by folate"; SBO:0000176, BiGG THFAT, EC 2.1.2.1, Rhea 34767, KEGG R02300, MNXR153595 |

**Code and data:**
- Data: [`model/curation/folate/folinic_acid_module.json`](../../curation/folate/folinic_acid_module.json), written by
  [scripts/write_curation.py](scripts/write_curation.py) from [before.json](before.json). before.json was read from
  `e584151f`; its SHA is in [logs/builds_sha256.txt](logs/builds_sha256.txt).
- Step: `platform/scripts/gem_annotate/folinic_acid_module.py` (`apply_folinic_acid_module`). The flag is in `cli.py`;
  `main.py` calls the step after the V-ATPase A1 step and before the biomass and crosswalk steps.
- Test: `platform/tests/test_folinic_acid_module.py` (7 tests, 24 subtests, no LPs).

**How the step works:**
- It pins the full definition of R748 and R990, including name, bounds, rule, stoichiometry, species, annotation,
  notes and groups.
- It pins the identity of m159, m160, H2O and H+ [C_cy], and the reactions of m159 and m160.
- It pins the annotation and reaction links of the three genes, and requires the folate group to exist once.
- The three targets must be all in the starting state or all in the curated state; anything else fails before any
  edit.
- A rerun, also on a saved and reloaded model, is a no-op.
- Inside a `with model:` block the edits are undone on exit, including name, annotation, notes and group membership.

## 2. Mechanism and evidence

| Claim | Tier | Source |
| --- | --- | --- |
| Methenyl-THF + H2O = 5-formyl-THF + H+ is Rhea 34767; its 7 reviewed UniProt enzymes are all SHMTs (E. coli and 6 cytosolic mammalian) | verified 2026-11-12 | [Rhea 34767](https://www.rhea-db.org/rhea/34767), UniProt REST |
| SHMT (cytosolic and mitochondrial rabbit, E. coli) catalyses it; SHMT-null E. coli has no 5-formyl-THF | verified, abstract only | [PMID 2201683](https://pubmed.ncbi.nlm.nih.gov/2201683/) |
| In *S. cerevisiae*, 5-formyl-THF is undetectable without SHMT; FAU1 is the single, non-essential MTHFS | verified, abstract only | [PMID 11923304](https://pubmed.ncbi.nlm.nih.gov/11923304/) |
| KEGG R02300 is filed under EC 2.1.2.10 / K00605 (glycine-cleavage T protein); ENZYME EC 2.1.2.10 describes only the lipoyl-aminomethyl → methylene-THF step | verified | [KEGG R02300](https://www.kegg.jp/entry/R02300), [ENZYME 2.1.2.10](https://enzyme.expasy.org/EC/2.1.2.10) |
| Yeast-GEM 9.0.2 `r_1030` (GCV1/2/3 + LPD1) and iMM904 `THFATm` (GCV1/2/3 + PDA1/PDB1) also carry glycine-cleavage-complex rules. Their subunit sets differ from each other and from iYali26's single T-protein rule. That iYali26 inherited its rule from this lineage is **inferred** | verified (rules); inferred (lineage) | local Yeast-GEM (`9fd2c572…`), [BiGG](https://bigg.ucsd.edu/models/iMM904/reactions/THFATm) |
| `YALI1D28569g` (yli:2910775) and `YALI1E19557g` (yli:2911591) are K00600 SHMTs; `YALI1F04082g` (yli:2907779) is K00605 | verified, orthology only | KEGG REST |
| Which SHMT paralog is mitochondrial | **untested**. KEGG GENES calls yli:2911591 (`E19557g`) "cytosolic", a weak hint only. Both genes go on both reactions, as R358/R359 already have them | — |
| R990 = R64 + R543 | verified (stoichiometry) | [before.json](before.json) |
| No 5-formyl-THF isomerase in KEGG (all reactions of C03479), Rhea (all reactions of ChEBI 57457) or ENZYME; BiGG FTHFI is in iMM904 only, gene-less | verified (searches) | [logs/isomerase_search_20261112.txt](logs/isomerase_search_20261112.txt) |
| R_THFAT charge imbalance −1, the same as R748. Methenyl-THF `m160`/`m803` is stored at charge 0, but C20H20N7O6 is the −1 species | verified | [ChEBI 57455](https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:57455); `check_mass_balance` in the focused test |

Gene names: no verified names. NCBI records all four proteins as "uncharacterized protein".
- `YALI1D28569g` and `YALI1E19557g`: serine hydroxymethyltransferase (KEGG GENES definitions "Serine
  hydroxymethyltransferase" and "…, cytosolic").
- `YALI1F04082g`: putative aminomethyltransferase, the T protein.
- `YALI1D01489g`: predicted 5-formyl-THF cyclo-ligase, the FAU1 orthologue.

All of these are orthology or homology assignments, not experimental evidence in *Y. lipolytica*. The model's UniProt
accessions for these four genes are deleted UniProt entries; this predates the change and is not fixed here.

## 3. Results

Condition: Task A medium `6eca70fd…` + PO1f `35307853…`, via `essentiality_simulation_context.py`. Gurobi at 1e-9,
basis reset before every solve. Inputs SHA-checked in [logs/sweep_inputs_sha256.txt](logs/sweep_inputs_sha256.txt).
Models: the control copy (`e584151f`) and the first build (`d4832b29`), the same optimisation problem as v2.

**Builds** ([logs/builds_sha256.txt](logs/builds_sha256.txt)):
- the control without the flag is byte-identical to `e584151f`;
- the seed-0 and seed-1 builds are byte-identical, both for the first build and for v2.

**Whole-model diff** against `e584151f`: [model_diff_vs_e584151f.json](model_diff_vs_e584151f.json), 11 entries,
[scripts/model_diff.py](scripts/model_diff.py). The script compares cobra-level fields: reactions (all fields plus
objective coefficients), metabolites, genes (name, annotation, links), groups, compartment names, and the model id,
name, annotation and notes. It does not compare gene notes or SBML content that cobra does not load. The byte-level
text diff of the two XML files ([logs/xml_text_diff_v2_vs_e584151f.txt](logs/xml_text_diff_v2_vs_e584151f.txt)) closes
that gap. The audit's independent diff of the first build found hunks only in R748, R990, R_THFAT and group50. The
entries are:
- R748: annotation (`ec-code` only), gene rule, name, notes;
- R990: bounds, notes;
- R_THFAT: added;
- the reaction links of the three genes;
- `group50` ("One carbon pool by folate"): R_THFAT added.

**Targeted solves** ([targeted_control.json](targeted_control.json), [targeted_candidate.json](targeted_candidate.json)):
- WT and the knockouts went through `checked_solve.solve`; all were optimal and passed the targeted check.
- The flux-range LPs were optimal but were not run through the residual and bound check. The candidate's R64 minimum
  at 100 % is −8.7e-10, below its lower bound of 0 but within the 1e-9 tolerance.

| | Control (e584151f) | Candidate |
| --- | --- | --- |
| WT | 0.31248856100769146 | 0.31248856100768774 (Δ 3.7e-15) |
| R64 min/max, no growth floor | 0 / 0 (blocked) | 0 / 55.47 |
| R_THFAT, no growth floor | — | 0 / 55.47 |
| R990, no growth floor | 0 / 0 | 0 / 0 (closed) |
| R748 = R989, no growth floor | 0 / 69.34 | 0 / 69.34 |
| R811 min/max, no growth floor (negative = export) | 0 / 0 | −0.83 / 0 |
| R991 min/max, no growth floor | 0 / 0 | 0 / 0.83 |
| R64, R_THFAT, R748, R989 max, at growth ≥ (1 − 1e-9) × WT | ≤ 7.1e-8 | ≤ 7.1e-8 |
| R542 (mitochondrial cyclohydrolase), any growth | 0 / 0 | 0 / 0 |

| Knockout (KO/WT) | Control | Candidate |
| --- | --- | --- |
| `YALI1F04082g` | 0.98138 (R748 among the 6 closed) | 0.98138 (R748 no longer closed) |
| `YALI1D28569g` or `YALI1E19557g` alone | 1.0 | 1.0 |
| both SHMTs | 0.98527 | 0.98527 (R748 and R_THFAT also closed) |
| `YALI1D01489g` (FAU1-type) | 1.0 | 1.0 |

**Full single-gene sweeps** (`ko_screen_certified.py`, unchanged; [sweep_compare.json](sweep_compare.json),
[sweep_control_summary.json](sweep_control_summary.json), [sweep_candidate_summary.json](sweep_candidate_summary.json);
TSVs in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/folinic_acid_module_20261112/`):
- 1,082 sweep LPs per model, no point failing the check, 20 non-optimal per model, all certified, identical in both.
- 1,081 genes; zero status changes and zero class changes in any rule, scope or cutoff. The largest difference in
  knockout growth is 3.1e-14.
- Rule B with mtDNA and untargeted genes out of scope, both models identical (TP / FN / FP / TN):

  | Cutoff | 1 % | 5 % | 10 % | 15 % (primary) |
  | --- | --- | --- | --- | --- |
  | Both models | 113/213/36/701 | 120/206/38/699 | 124/202/43/694 | 141/185/43/694 |

- **Positive coverage:** 326 of the 1,612 consensus positives are model genes. None of them is out of scope (TP + FN =
  326 in both scopes), so 1,286 positives are outside the model.
- Negatives are model genes not on the consensus list (the user-defined class). The list was used in development, so
  these scores are development references, not independent validation.

**Pytest.** Root `.venv`, Python 3.13, cobra 0.30.0, pytest 9.0.2, `--continue-on-collection-errors`, short-path
archive copies ([logs/pytest_env.txt](logs/pytest_env.txt)):

| Run | Passed | Failed | Errors | Skipped |
| --- | --- | --- | --- | --- |
| Before (HEAD `154eb64c`) | 463 | 3 | 18 | 24 |
| After (first build's data) | 470 | 3 | 18 | 24 |

[logs/pytest_compare.txt](logs/pytest_compare.txt) compares the junit files test by test. The only differences are the
7 new `test_folinic_acid_module` tests, all passed. The 3 failures and 18 errors are the standing ones, identical before
and after. After the note rewording, the focused test passes again (7 / 24 subtests,
[logs/pytest_focused_v2.txt](logs/pytest_focused_v2.txt)); the full suite was not re-run for the note change.

## 4. What it means

- **The dead end is gone.** Cytosolic 5-formyl-THF now has a source, the SHMT side reaction supported by orthology
  and by E. coli, rabbit and *S. cerevisiae* evidence. R64 (FAU1-type) can carry flux.
- **R990 is closed.** Open and gene-less, it would have bypassed a FAU1 knockout once R_THFAT existed (inferred).
  With FAU1 knocked out, 5-formyl-THF can still leave through R811/R991 export.
- **R748 now carries the SHMT orthologues instead of the T protein.** The T-protein knockout closes 5 reactions, not 6.
- **Predictions do not change.** WT, every knockout class and every score are the same. This was expected: no biomass
  term contains 5-formyl-THF.
- **At the optimum the module is idle.** The module reactions carry ≤ 7.1e-8 only within the 1e-9 growth slack:
  R_THFAT + R64 and R748 + R989 are ATP-hydrolysing loops, and R991 export throws folate away. Below the optimum they
  can run.
- **What would make the module carry net flux at maximal growth:** a 5-formyl-THF demand or a folinic-acid uptake.
  Neither was added, and no *Y. lipolytica* data for either was searched in this round.
- **Inactive-reactions rule.** The rule now on `main` (AGENTS.md §5, PR #50; see §7) is respected on two points: no
  exchange, transport or sink was added, and the regression sweep was run.
  - **One tension.** The rule calls a (0, 0) reaction a gap, and closing R990 creates one. The closure follows the
    earlier practice for lumped duplicates (R1725, R2178, and the lumped Kennedy steps R353/R1846/R1771/R2165): keep
    the reaction for traceability.
  - **Why not delete it.** Deletion would need evidence that the reaction does not exist in the organism, and here
    the chemistry exists, as R64 + R543.
  - **Decision for the user:** keep R990 closed, or delete it.

## 5. Could the method be improved, and what was not tested

- **Exact-optimum ranges.** The 100 % ranges use a (1 − 1e-9) × WT floor and unchecked LPs. Exact-optimum FVA through
  the checked solve would separate "zero at the optimum" from "within the slack". The expected answer is zero; this is
  inferred, not run.
- **Untested:**
  - which SHMT paralog is mitochondrial (TargetP / DeepLoc on both would settle the rule);
  - literature beyond abstracts (the Holmes & Appling 2002 full text returned HTTP 403 to WebFetch and curl in this
    session);
  - other media;
  - the full pytest suite on the v2 note text.
- **Out of scope; recorded as follow-ups:**
  - methenyl-THF `m160`/`m803` charge 0 vs ChEBI −1: a charge change that propagates to R542, R543, R549–R551, R64,
    R989, R748, R1907 and R_THFAT;
  - R542 (mitochondrial cyclohydrolase) is blocked at any growth in both models; cause not traced;
  - the 10-formyl-THF metabolites `m532`/`m534` still carry 5-formyl-THF InChI, MetaCyc, SEED and VMH identifiers;
  - R542 carries the R02300 / Rhea 34767 cross-references of the 5-formyl-THF reaction instead of
    [R01655](https://www.kegg.jp/entry/R01655);
  - R748 still lists the cytosolic BiGG ID `THFAT`, now also on R_THFAT;
  - the R989 gene (F3) stays with the lipid-unlump line;
  - the libSBML "Setting notes" error on `R_R304` appears in every build log, the control included, and predates this
    change.

## 6. Verification levels

| Statement | Level |
| --- | --- |
| R64 was blocked in e584151f (while R991's lower bound is 0) | provable from the bounds; observed in the control FVA |
| The control build reproduces e584151f; candidates are seed-independent | observed (byte comparison) |
| Only the intended entities changed | statically verified (cobra-level diff plus byte-level text diff) |
| v2 differs from the first build only in two note texts | statically verified (byte-level diff) |
| WT, the targeted knockouts and the full sweep are unchanged on Task A + PO1f | observed in this run (on the first build; same problem as v2) |
| The module carries no flux at the exact optimum | inferred (≤ 7.1e-8 within the 1e-9 slack) |
| SHMT catalyses R748 / R_THFAT in *Y. lipolytica* | orthology plus other organisms; untested in *Y. lipolytica* |
| Pytest: only the new tests differ | observed (junit, test by test) on the first build's data; focused test re-run on v2 |

## 7. Branch base

**`main` has moved.** It is now at `efd7f780`, three PRs past `154eb64c`: the R2121 and R1172 closures (PRs #50 and
#51) and the AGENTS.md §5 inactive-reactions rule. All three landed before these builds.

**Effect on this candidate.** It is a sibling of those candidates and lacks both closures.

**Effect on the code.** `main.py` and `cli.py` on main add parameters at the same insertion points, so merging this
work needs a textual conflict resolution, then a combined build.

None of that was done in the first round. **Resolved in §8:** the user chose "Merge main, rebuild".

## 8. Rebuild on current main (2026-11-12, user: "Merge main, rebuild")

**Merge.**
- This round's code, data, tests and docs were committed as `c00ec033`. The two `e584151f`-based candidate XMLs were
  not committed.
- `origin/main` (`cf7b8dd9`: PRs #50–#52, the R2121 and R1172 closures and the legacy-biomass plan) was merged as
  `768c187e`.
- The conflicts in `cli.py`, `main.py` and `STATE.md` were resolved by keeping both sides. The folinic step still runs
  after the V-ATPase A1 step and before main's R2121 and R1172 steps; the targets of the three steps are disjoint.
- Not pushed.

**New baseline.** `..._mcnaugh_rxns_v3_r2121_close_r1172_close_20261112.xml`,
`684f3752764fb8f25219008966953863141724f34a451dffb77c235aa3cb2f62` (main's newest candidate). Its argv is
[scripts/args_684f3752.txt](scripts/args_684f3752.txt): the e584151f argv plus `--r2121-malonyl-coa-close
--r1172-hmgcoa-transport-close`.

**Candidate.** `..._mcnaugh_rxns_v3_r2121_close_r1172_close_fthf_20261112.xml`,
`353cefdc16ef1dc510eae49c168fc999254c2652c3994287af45589621a90cd9`. That is the 684f3752 argv plus
`--folinic-acid-module`, built with [scripts/run_build_684f3752.zsh](scripts/run_build_684f3752.zsh). The name uses
`fthf` because `folinic_acid` pushed the builder's temporary file past the macOS 255-byte limit; the first attempt
failed with `Errno 63`.

| Check | Result |
| --- | --- |
| Control build (684f3752 argv, merged code, no flag) | byte-identical to `684f3752` |
| Candidate seed 0 vs seed 1 | byte-identical (`353cefdc`) |
| Whole-model diff vs `684f3752` ([model_diff_vs_684f3752.json](model_diff_vs_684f3752.json)) | the same 11 entries, field for field equal to the v2-vs-e584151f diff (`json ==` True); byte-level hunks identical apart from line numbers ([logs/xml_text_diff_vs_684f3752.txt](logs/xml_text_diff_vs_684f3752.txt)) |
| WT (Task A + PO1f, Gurobi 1e-9) | 0.3111005359259494 → 0.311100535925948 (Δ 1.4e-15) ([targeted684_control.json](targeted684_control.json), [targeted684_candidate.json](targeted684_candidate.json)) |
| R64, no growth floor | 0 / 0 → 0 / 55.47; R_THFAT 0 / 55.47; R990 0 / 0; R748 = R989 0 / 69.34; R542 still 0 / 0 |
| Module at (1 − 1e-9) × WT | ≤ 7.1e-8 (unchecked LPs, as before) |
| Targeted knockouts | `YALI1F04082g` 0.98243 both (R748 no longer closed); both SHMTs 0.98534 both; single SHMTs and `YALI1D01489g` 1.0 |
| Full certified sweep, both models ([sweep684_compare.json](sweep684_compare.json)) | 1,082 LPs each, no check failures, 20 non-optimal, all certified; zero status or class changes; max KO growth difference 2.8e-14 |
| Rule B, mtDNA + untargeted out of scope, TP/FN/FP/TN at 1/5/10/15 % | 117/209/36/701, 124/202/38/699, 128/198/43/694, **145/181/43/694**, identical in both; 326 of 1,612 positives in the model, none out of scope |

**Pytest on the merged code.** Root `.venv`, `--continue-on-collection-errors`, short-path archive copies:

| Run | Passed | Failed | Errors | Skipped |
| --- | --- | --- | --- | --- |
| Before (`origin/main` `cf7b8dd9`) | 477 | 3 | 18 | 24 |
| After (merged HEAD `768c187e`) | 484 | 3 | 18 | 24 |

[logs/pytest684_compare.txt](logs/pytest684_compare.txt) compares the junit files test by test. The only differences
are the 7 new `test_folinic_acid_module` tests, all passed. This run uses the v2 note text, so it also closes the
"full suite not re-run for the note change" gap in §3.

The superseded candidates `d4832b29` (first build) and `85c623c8` (v2 on e584151f) stay untracked in
`model/candidates/` until the user decides.
