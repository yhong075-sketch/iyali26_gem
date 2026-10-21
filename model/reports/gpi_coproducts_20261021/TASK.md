# GPI co-products: DAG on R373, Dol-P on R375, R_GPIGNT note wording (2026-10-21)

## Authorization

User, 2026-10-21, after reading the plan in this session:

> Approved: implement the GPI co-product plan as written (new --gpi-coproducts flag requiring
> --gpi-chain-links, R373 + DAG, R375 + Dol-P, R_GPIGNT note fix, KO comparison included; R376 and
> the DAG formula left as findings). No commit or push.

## Starting state

- **Branch:** `claude/gpi-chain-gaps-398472` at `a78a49d`. The previous task's work is uncommitted
  in the same worktree: [gpi_chain_links_20261021](../gpi_chain_links_20261021/REPORT.md).
- **Baseline model:**
  `model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_links_20261021.xml`, SHA-256
  `89559a0f04fa38a58a563d6dc8b2f9279a933b810461e7fa731359ad670db718`. It is `164f235f` plus
  `--gpi-chain-links`. It is a candidate, not a reference.
- **Before-state of the targets and their neighbours:** [before.json](before.json), written by
  [scripts/before_state.py](scripts/before_state.py).

## Targets

| Item | Before | After |
| --- | --- | --- |
| R373 (GPI step 6, EtNP to Man1; GPR YALI1F27343g) | m598 + m599 ⇌ m600 | m598 + m599 ⇌ m600 + **m604[C_er]** (DAG) |
| R375 (GPI step 8, third mannose; GPR YALI1D05270g) | m502 + m601 ⇌ m602 | m502 + m601 ⇌ m602 + **m503[C_er]** (Dol-P) |
| `R_GPIGNT` note `gpi_chain_links_rationale` | "…The 2 H+ reconcile the stored charges (UDP-GlcNAc -2, UDP 0); KEGG does not balance protons…" | Says that the 2 H+ are needed because the model stores m403 at −2 and m11 at 0 (as R778 and R211 do), and that KEGG writes the equation for neutral species |

- Bounds, GPRs, names and annotations of R373 and R375 stay as they are; only `gpi_coproducts_*`
  notes are added.
- The co-products are applied by `apply_gpi_coproducts()` in
  `platform/scripts/gem_annotate/gpi_coproducts.py`, under the flag `--gpi-coproducts`. That flag
  requires `--gpi-chain-links`. Data: `model/curation/chemistry/gpi_coproducts.json`.
- The note is changed in `model/curation/gapfill/gpi_chain_links.json`. As a result, a
  `--gpi-chain-links`-only build no longer reproduces `89559a0f`; it differs only in that note.

## Evidence

All verified, read on 2026-10-20 and 2026-10-21.
- **R373:** KEGG [R12351](https://www.kegg.jp/entry/R12351) (GPI-EtNP transferase 1,
  [K05285](https://www.kegg.jp/entry/K05285)) gives PE + GPI → EtNP-GPI + 1,2-diacylglycerol. All
  five PE-donor EtNP transfers in KEGG map00563 release DAG. R373's gene YALI1F27343g carries
  `kegg.genes` yli:2908691, which KEGG assigns to K05285. That is an ortholog assignment, not
  experimental validation.
- **R375:** KEGG [R05922](https://www.kegg.jp/entry/R05922) (GPI mannosyltransferase 3,
  [K05286](https://www.kegg.jp/entry/K05286)) gives Dol-P-Man + G00148 → Dol-P + G00149. G00148 is
  (Man)2 (EtN)1 and G00149 is (Man)3 (EtN)1, the same compositions as m601 → m602. R375's gene
  YALI1D05270g carries yli:2910629, assigned to K05286. R372 and R374 in the model already return
  Dol-P (m503).

## Expected effect

- WT growth is unchanged: the exit stays open and the chain carries 0 at the optimum.
- With the exit forced, R_DAGt_er_em carries twice the chain flux, R375 returns Dol-P, and the
  growth cost of forced GPI flux falls.
- The R375 balance is checkable only as a difference: m502 − m503 = C6H10O5, one mannosyl residue.
  The GPI intermediates have no formulas.
- R373 is not checkable at all, because the model's DAG formula (C6H8O6) is not a diacylglycerol
  (finding below).

## Acceptance criteria

1. pytest per test equals the last run (217 passed, 11 failed, 17 errors, 31 skipped) plus the new
   tests passing.
2. **Builds** (combined argv):
   - with no GPI flag: reproduces `164f235f`;
   - `--gpi-chain-links` only: differs from `89559a0f` only in R_GPIGNT's `gpi_chain_links_rationale`
     note;
   - both flags: the same bytes under `PYTHONHASHSEED=0` and `1`.
3. Whole-model diff of the new candidate against the links-only build: only R373 and R375
   stoichiometry and their new notes. Against `89559a0f`: the same plus the R_GPIGNT note.
4. On SD-Leu + PO1f with GLPK:
   - WT growth equals 1.430806;
   - the chain FVA maximum is > 0;
   - with the exit forced to ≥ 1e-5 the solve is optimal and R_DAGt_er_em is 2e-5;
   - closing each link in turn makes the forced exit infeasible.
5. Single-gene knockouts compared with `89559a0f`: no class change at 1/5/10/15 %. Any gene whose
   status or growth (by more than 1e-6) differs is solved again in a fresh context and reported.

## Compute budget and stop condition

- **Budget:** 4 builds; about 50 diagnostic LPs including FVA; 2 × (1,074 + 1) knockout solves;
  fresh re-solves of differing genes (at most 10).
- **Stop and report on:**
  - a guard rejection;
  - a diff outside the targets;
  - seed-dependent bytes;
  - a WT change;
  - a test regression;
  - a knockout class change that persists in a fresh solve.

## Environment

As in the links task:
- **Interpreter:** the main checkout's `.venv` (Python 3.13.5, cobra 0.30.0, optlang 1.8.3,
  pytest 9.0.2). pytest runs with `--continue-on-collection-errors`, because `pydantic_ai` is
  absent.
- **Research root:** `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`.
- **Condition:** `sd_leu.csv` (`ed176d26…`) and `po1f_sd_leu.json` (`35307853…`).

## Findings recorded, not changed (user decision)

- **R376 stoichiometry:** its product m603 (CHEBI:53060) carries three EtNP groups; m602 carries
  one. KEGG splits the two additions ([R08107](https://www.kegg.jp/entry/R08107),
  [R12702](https://www.kegg.jp/entry/R12702)), each using one PE.
  - So either R376 lumps two steps and lacks one PE and one DAG, or m603 points to the wrong
    compound.
  - R376's GPR does not settle it:
    - YALI1E04744g is K05287 (GPI11), a subunit shared by both enzymes;
    - YALI1E14925g has no identifiers;
    - KEGG's GPI7 and GPI13 genes (yli:2911071, yli:2912812) are not in the model's gene
      annotations.
- **DAG species (m604, m1649 and the other diglyceride pools):** formula C6H8O6, with an InChI for
  a different small molecule. A diacylglycerol carries two acyl chains. Every DAG reaction checked
  (R1725, R1726, R1728, R1773) is unbalanced.

## Corrections after the independent audit (2026-10-21)

The text above is kept as written. These corrections apply. They were verified this session; see
[REPORT.md §7–8](REPORT.md#7-independent-audit).

- **R376 finding:**
  - YALI1E14925g is GPI13. The NCBI feature table maps YALI1_E14925g to GeneID 2912812 and
    XP_503848.3, which is KEGG yli:2912812, K05288. So KEGG's GPI13 gene is in the model, without
    identifiers.
  - CHEBI:53060's definition describes two EtNP groups, and its KEGG cross-reference is G13044
    ((Man)3 (EtN)2), the product of R08107.
  - So R376 most likely equals R08107 as written. The likelier problems are:
    - m603's three-EtNP name;
    - a missing GPI7 step (R12702; YALI1_D34684g is not in the model);
    - possibly an OR GPR where the GPI11 and GPI13 subunits should be AND.
- **DAG finding:**
  - C6H8O6 cannot be a diacylglycerol; that stands.
  - R1725 and R1726 are imbalanced but closed at (0, 0).
  - R1728 and R1773 contain formula-less species (m1656, m378), so their balance is not computable.
- **Expected effect:** "the chain carries 0 at the optimum" holds for WT, not for every knockout. In
  nine knockouts that share KO/WT 0.1277, the open exit lets the chain turn surplus PE into DAG. The
  gain is +8.477e-7, which is below the 1e-6 tolerance of criterion 5. No class changes.
