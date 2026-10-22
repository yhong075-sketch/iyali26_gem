# TPP requirement: TPP in biomass_C and R1342 closed (2026-10-22)

## Authorization

User, 2026-10-22, after reading the plan in this session:

> Approved: implement the TPP plan as written: --tpp-requirement on top of cf0825b8, TPP in
> biomass_C at the yeast-GEM value 1.2e-6 mmol/gDW, R1342 closed (not deleted), ThiL BLAST check,
> KO comparison included. No commit or push.

The plan answered the user's request "I do the test on claude science, can you fix the gap as
figure showed?". The figure was a Claude Science in-memory test titled "Gap 4 · TPP exists but
nothing requires it". It showed that no reaction drains cytosolic TPP to a sink: its consumers
R753 and R933 end in dead ends. (Wording corrected after the audit; the first version said that no
reaction consumes cytosolic TPP, which is wrong.) With TPP required in
biomass, a THI80 knockout is bypassed through the gene-less TMP kinase R1342, and closing R1342
makes the knockout lethal.

## Starting state

- **Branch:** `claude/science-test-gap-fix-13f440` at `abe7e40` (main after PR #13). Clean worktree.
- **Baseline model:** `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_coproducts_20261021.xml`,
  SHA-256 `cf0825b86c3533d8a12f697a994eba3d740671a68a367143ee32cbfdc8686b7c`.
  - 2,311 reactions, 1,876 metabolites and 1,074 genes.
  - It is `164f235f` plus `--gpi-chain-links --gpi-coproducts`.
  - It is a candidate, not a designated reference.
- **Before-state of every target and its neighbours:** [before.json](before.json), written by
  [scripts/before_state.py](scripts/before_state.py).

## What was reproduced before any change (2026-10-22, in memory, GLPK)

Condition: the Task A medium + PO1f (the benchmark default since 2026-10-21). The figure's WT growth,
0.31243, matches this condition. Script: [scripts/reproduce_figure.py](scripts/reproduce_figure.py);
output: [reproduce_figure_cf0825b8.json](reproduce_figure_cf0825b8.json).

| Test | Growth | Flux on the thiamine map |
| --- | --- | --- |
| As supplied | 0.3124267446 | 0 at the optimum. R752, R753, R933, R1342 and R934 have FVA maximum 0 at fraction 0. |
| THI80 (YALI1E25443g) knocked out | 0.3124267446, KO/WT 1.000 | 0 |
| TPP at 1e-5 in biomass_C | 0.3124267101 | R934 = R752 = 3.124e-6 |
| …and THI80 knocked out | 0.3124266825, KO/WT 1.000 | R934 → R1338 → R408, R413, R645 → R754 → R1342 |
| …and R1342 closed | 4.7e-29 (optimal) | none |
| TPP at 1e-5, R1342 closed, no knockout | 0.3124267101 | R934 = R752 |

Two more in-memory checks:
- With TPP at 1e-5 and thiamine uptake R1340 closed, growth is 0.3124253, so the model can make TPP
  de novo.
  pFBA then routes it R970 → R413/R645 + R1344 → R754 → **R1342**, so the thiamine-free case
  must be rechecked with R1342 closed (through R1337 → R752).
- With R1342 closed and THI80 knocked out, KO growth is about 0 for both coefficients: 1e-5 gives
  4.7e-29 and 1.2e-6 gives −7.4e-30 in the script run (an earlier scratch run gave −5.2e-15 and
  2.8e-16).

## Targets

| Item | Before | After |
| --- | --- | --- |
| `biomass_C` | no TPP term | consumes `m975[C_cy]` (TPP, "TDP_") at 1.2e-6 mmol/gDW; `tpp_requirement_*` notes added |
| `R1342` (thiamine-phosphate kinase, ATP + TMP → ADP + TPP, gene-less) | bounds (0, 1000) | bounds (0, 0); equation, name, annotation and GPR kept; `tpp_requirement_*` notes added |

- No reaction, metabolite or gene is added or removed.
- Curated data: `model/curation/thiamine/tpp_requirement.json`.
- Guarded step: `apply_tpp_requirement()` in `platform/scripts/gem_annotate/tpp_requirement.py`,
  under the opt-in flag `--tpp-requirement`. The flag requires `--gpi-coproducts`, the baseline it was
  reviewed on.
- Focused test: `platform/tests/test_tpp_requirement.py`.

## Evidence

- **No KEGG ThiL (TMP kinase) assignment in either yeast (verified 2026-10-22; heading narrowed
  after the audit, since this is not proof that TMP kinase activity is absent):** the KEGG REST queries
  [link/yli/K00946](https://rest.kegg.jp/link/yli/K00946) and
  [link/sce/K00946](https://rest.kegg.jp/link/sce/K00946) (thiL, EC 2.7.4.16) return nothing.
  KEGG does list thiamine pyrophosphokinase (K00949) in both:
  [yli:2911892](https://rest.kegg.jp/link/yli/K00949) (YALI2_E01340g) and sce:YOR143C (THI80).
  KEGG's *Y. lipolytica* genome is CLIB122, not W29, so the BLAST below checks W29 directly.
- **TPP coefficient (verified 2026-10-22):** yeast-GEM 9.1.1, cofactor pseudoreaction `r_4598`,
  consumes `s_1475` "TDP" at 1.2e-6 mmol/gDW. `s_1475` is thiamine diphosphate (KEGG C00068, BiGG
  thmpp, C12H16N4O7P2S). The value is for *S. cerevisiae*, not *Y. lipolytica*; it is provisional.
  Source: [yeast-GEM.yml](https://github.com/SysBioChalmers/yeast-GEM/blob/main/model/yeast-GEM.yml),
  raw file SHA-256 `2ddfaaa7…`.
- **Yeast TPP route:** TMP is dephosphorylated to thiamine, which THI80 pyrophosphorylates to TPP.
  **Established/uncited**; no source was opened for it in this task.
- **The THI80-like gene:** YALI1E25443g; verified name: none verified; thiamine pyrophosphokinase by
  similarity ("similar to *S. cerevisiae* THI80 (YOR143C)" in the consensus file); model GPR on R752
  and R753. It is on the consensus essential list (at least 2 of 3 screens). That list was used in
  model development, so it is not independent evidence for this change.

## Expected effect

- WT growth falls by about 4e-9 (TPP synthesis cost). R934 and R752 carry flux at the optimum.
- R1342 can no longer carry flux. It was already blocked in the baseline, because TPP had no route
  to a sink (wording corrected after the audit).
- A THI80 knockout becomes lethal.
- Other knockouts are expected to be unchanged: every other route to TPP passes through R752, and the
  medium supplies thiamine. This is checked for all genes.

## Acceptance criteria

1. **ThiL BLAST:** no W29 protein whose reciprocal best Swiss-Prot hits are ThiL (EC 2.7.4.16).
   A credible ThiL hit stops the work.
2. **Tests:** pytest per test equals the baseline run plus the new test file passing.
3. **Builds:** the cf0825b8 argv without the new flag reproduces `cf0825b8…` byte for byte. With the
   flag, `PYTHONHASHSEED=0` and `=1` give the same bytes.
4. **Whole-model diff against `cf0825b8`:**
   - `biomass_C` gains `m975[C_cy]` at −1.2e-6 and the `tpp_requirement_*` notes;
   - `R1342` has bounds (0, 0) and the `tpp_requirement_*` notes;
   - the derived reaction list of `m975[C_cy]` gains `biomass_C`.

   Any other difference stops the work.
5. **Solves (Task A medium + PO1f, GLPK):**
   - WT is optimal;
   - at the optimum TPP comes through R752;
   - THI80 knocked out: KO/WT below 0.01;
   - thiamine uptake closed: growth above 0.
6. **Knockouts:** all 1,074 genes on both models at 1%, 5%, 10% and 15% KO/WT. Every class change is
   listed with its cause. A class change other than YALI1E25443g is reported, never fixed, and stops
   the write-up until it is explained.

Mass balance is not an acceptance criterion here. biomass_C is a pseudo-reaction, and R1342 keeps
its pre-existing imbalance (TMP finding below) while closed.

## Compute budget and stop condition

- **BLAST:** Swiss-Prot ThiL (EC 2.7.4.16) queries against the W29 proteome at E ≤ 1e-5, a
  sensitivity rerun at E ≤ 10, and a reciprocal Swiss-Prot check (top 5 hits) of at most 5 W29 hits
  per query set. Local databases from `orphan_gene_search_20261019`; nothing is downloaded.
- **LPs:** about 30 diagnostic solves plus 2 × (1,074 genes + WT) knockout solves. Same medium,
  strain profile, objective and solver on both models.
- **Stop and report on:**
  - a guard rejection;
  - a ThiL ortholog in W29;
  - a diff outside the targets;
  - seed-dependent bytes;
  - a non-optimal WT;
  - a test regression;
  - a knockout class change other than YALI1E25443g.

## Environment

- **Interpreter:** the main checkout's `.venv` (Python 3.13.5, cobra 0.30.0). BLAST+ 2.17.0 from the
  conda env `seqsearch`.
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- **Condition:** `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (`6eca70fd…`) and
  `state/strain_profiles/po1f_sd_leu.json` (`35307853…`).

## Not in scope (findings only)

- TMP `m977` is annotated as 3,5,6-trichloro-2-methoxypyridine (C6H4Cl3NO), so R754, R1337 and
  R1342 are not mass-balanced.
- TPP `m975[C_cy]` and `m1116[C_mi]` have no formula or cross-reference.
- Mitochondrial TPP has no consumer: R933 stays blocked.
- The loop R754 → R1337 → R1338 → R408/R413/R645 → R754 can carry flux (FVA maximum 21.3). It burns
  ATP and makes no TPP.
- The default build, other media and a benchmark score.
