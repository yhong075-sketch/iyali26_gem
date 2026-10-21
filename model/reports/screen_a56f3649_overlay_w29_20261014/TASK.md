# TASK: a56f3649 screen under Task A's diagnostic overlay, PO1f + CSM-Leu, glucose 2.414 (2026-10-14)

## Authorization (user, 2026-10-14, this conversation)

The user showed Task A's row for the parent `0ac2addf` (PO1f + CSM-Leu, glucose 2.414; in-memory
diagnostic package with assumed caps R1384 ≤ 8 and R860 ≤ 1; 27 unresolved, all counted as lethal)
and asked: "could you screen at conditions like these for our newest model?" After the plan review
they approved: "screen a56f3649 on PO1f + CSM-Leu glucose 2.414 with the in-memory holistic overlay
(R1384 <= 8, R860 <= 1), including the 4 reaction controls and the NGAM diagnostic, budget at most
1,228 LPs; report the benchmark rule and the lethal view, compare with Task A's 0ac2addf overlay
table, and" (the message was cut there; the approved plan continues "draw the two 15 % panels").

Computation only. The overlay exists in memory only; no model, GPR, medium, label or baseline file
changes; no commit or push.

## Starting state

- Worktree `confusion-matrices-iyali21-7e63ca`, branch `claude/confusion-matrices-iyali21-7e63ca`,
  HEAD `cc278c4` (this session's confusion-matrix commit on top of `main` `4d97b82`), clean.
- Research root `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research` (approved
  earlier in this conversation; `IYALI26_RESEARCH_ROOT` was unset).
- No reference model is designated.

## Inputs

| Input | Path | SHA-256 |
| --- | --- | --- |
| Model | `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` (newest candidate; not a reference) | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` |
| Medium | `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (glucose 2.414) | `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| Strain profile | `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` | `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` |
| Experimental positives | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv` | `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b` |
| mtDNA scope rule | `model/expected/mtdna_genes_out_of_scope.json` | `16b3e9b4418f24821b8775936fd377ae21f77afa8e8af2b8949a246f62c44996` |
| Overlay | `model/reports/r132_mito_adh_20261012/holistic_overlay.py` (verbatim copy of Task A's `apply_overlay`) | `bf6c7ce54d8338c365d98c82ca345ff2bd8c6055b38d402e8144a16b4abf30c4` |
| Task A overlay screen (comparison only) | `$IYALI26_RESEARCH_ROOT/artifacts/tasks/w29_glucose_uptake_20261012/gene_screen_holistic_overlay_diagnostic.tsv` (`0ac2addf`) | `e9f4ca6f9a3dcd74bf142cb80b1c33f9c03eecbe2aa495821611c4230a2e42ca` |

**Overlay contents** (diagnostic, in memory): R1864/R1866 rewritten from 2 NAD(H) to 1; R2069
rewritten to KEGG R01752; R1882, R1883, R1884, R2255, R2259, R2263, R_MDHy and R_ALCDH_nadp_hi
closed; R1113 upper bound 0; R171 GPR reduced to its 14-gene AND core; **R1384 ≤ 8 and R860 ≤ 1,
ASSUMPTION (diagnostic)**. Task A audit F3: the R1864/R1866 rewrite leaves +1 H and +1 charge,
and m1949 has no formula, so overlay results are diagnostic only.

**Pre-check (this session, no solve):** the overlay's structural assertions pass on a56f3649
(R132 is not referenced), and after the overlay the only forced bound is
`xMAINTENANCE` ≥ 7.8625.

## Method

[run_screen.py](run_screen.py), adapted from
[../confusion_matrices_csmleu_w29_20261014/run_screen.py](../confusion_matrices_csmleu_w29_20261014/run_screen.py):
same solver settings, static FBA, objective `biomass_C`, strict unrounded KO/WT < cutoff,
failed solves unresolved. Added:
1. the overlay, applied after the effective context is built, with its file SHA checked;
2. Task A's four reaction controls (R2062, R305, R304, R171 closed);
3. Task A's NGAM diagnostic for every knockout with no solution (genes including mtDNA, and
   controls): re-confirm infeasible with the NGAM floor; set only the `xMAINTENANCE` lower
   bound to 0 and re-solve growth; maximize `xMAINTENANCE`. Label "cannot meet maintenance" if
   the relaxed problem is feasible and the maximum is below 7.8625; otherwise "unresolved, other
   cause".

[summarize.py](summarize.py) (no solves) scores three views, compares every gene with Task A's
`0ac2addf` overlay table and draws the figure.

## Budget and stop conditions

- WT + 1,073 gene KOs + 4 controls = 1,078 LPs; NGAM diagnostic at most 150 LPs (3 per knockout
  with no solution). Total at most 1,228 LPs. One thread, 60 s per LP, 900 s wall. No retries.
- Stop on any input, overlay or code hash change, an overlay assertion failure, an invalid WT, or
  the budget being reached. A failed or non-finite solve is unresolved, not lethal.

## Scoring (benchmark.md)

- **Benchmark rule (primary):** a knockout with no solution is unresolved.
- **Lethal view (labelled second view, user decision 2026-10-13 for Task A):** every knockout
  with no solution is counted as lethal (KO/WT = 0).
- **Maintenance view (labelled, Task A §5 proposal, not adopted):** only "cannot meet
  maintenance" knockouts are counted as lethal; others stay unresolved.
- Cutoffs 1, 5, 10 and 15 % (15 % primary). Negatives are model genes not on the list
  (user-defined class). The 16 mtDNA genes are out of scope. Development reference, not
  independent validation.

## Outputs

- This folder: TASK.md, scripts, `screen.manifest.json`, `summary.json`, the figure (PNG and
  PDF), REPORT.md.
- Raw tables: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/screen_a56f3649_overlay_w29_20261014/`.
- A `model/STATE.md` entry.
