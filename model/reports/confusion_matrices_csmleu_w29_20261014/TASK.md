# TASK: confusion matrices for iYali21 and the newest candidate, PO1f + CSM-Leu, glucose 2.414 (2026-10-14)

## Authorization (user, 2026-10-14, this conversation)

Request: "plz using the iyali21 and our newest model to draw two confusion matrix", with the
condition "PO1f + CSM-Leu, glucose uptake 2.414 (W29)". After the plan review the user approved:
"use the research root you found, rerun both models fresh on PO1f + CSM-Leu glucose 2.414
(2,133 LPs), and draw the two 15% matrices as planned."

Computation only. No model, GPR, medium, label or baseline change; no commit or push.

## Starting state

- Worktree `confusion-matrices-iyali21-7e63ca`, branch `claude/confusion-matrices-iyali21-7e63ca`,
  HEAD `4d97b828ac2b0dc95d592060605af75069705c09` (= `main` after PR #5), clean.
- `IYALI26_RESEARCH_ROOT` was unset. The user approved
  `/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research`; the experimental list and
  the PO1f profile match the SHA-256 values in [benchmark.md](../../expected/benchmark.md)
  (checked 2026-10-14).
- No reference model is designated ([candidates/README.md](../../candidates/README.md)).
- Related work: a sibling session drew the same two matrices on SD-Leu (glucose 10) today, commit
  `01beb03` on branch `claude/confusion-matrices-iyali21-fc6ecb` (not merged here). Its
  `run_screen.py` and `summarize.py` are the starting point for this folder.

## Models

| Label | File | SHA-256 | Why |
| --- | --- | --- | --- |
| `iyali21` | `model/source/iyli21.xml` | `6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631` | Original iYali21 |
| `a56f3649` | `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml` | `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515` | Newest built candidate (Task B, R132 removed, 2026-10-13); newest is not the same as reference |

iYali21 had never been screened on this medium. Task B screened `a56f3649` on it
([r132_mito_adh_20261012/run.json](../r132_mito_adh_20261012/run.json)).

## Inputs

| Input | Path | SHA-256 |
| --- | --- | --- |
| Medium | `model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv` (glucose R1070 = 2.414) | `6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45` |
| Strain profile | `$IYALI26_RESEARCH_ROOT/state/strain_profiles/po1f_sd_leu.json` | `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383` |
| Experimental positives | `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv` | `1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b` |
| mtDNA scope rule | `model/expected/mtdna_genes_out_of_scope.json` | `16b3e9b4418f24821b8775936fd377ae21f77afa8e8af2b8949a246f62c44996` |
| Task B raw screen (comparison only) | `$IYALI26_RESEARCH_ROOT/artifacts/tasks/r132_mito_adh_20261012/runs/gene_screen_candidate.tsv` | `c9185ecdb26fcc99620b7522b091640b1938397ead865c2c46958d4427c30024` |

## Method

- [run_screen.py](run_screen.py): the sibling SD-Leu runner with only the medium path and SHA
  changed. Same Gurobi parameters, static FBA, objective `biomass_C`, strict unrounded KO/WT <
  cutoff, failed solves unresolved.
- [summarize.py](summarize.py) scores both screens from the raw rows and draws the figure. Its
  SD-Leu comparisons are replaced by one no-solve check: the `a56f3649` rerun against Task B's raw
  table on the same medium. It also reports a second, labelled view that counts "no solution"
  as lethal, as Task B did.

## Budget and stop conditions

- `a56f3649`: 1 WT + 1,073 KOs = 1,074 LPs. `iyali21`: 1 WT + 1,058 KOs = 1,059 LPs.
  Total 2,133 LPs. One thread, 60 s per LP, 900 s wall per model. No retries.
- Stop on any input or code hash change, an invalid WT, or the budget being reached. A
  failed or non-finite solve is unresolved, not lethal. An infeasible WT is reported, never
  fixed by tuning the medium.

## Scoring (benchmark.md)

- Primary cutoff 15 % (drawn); 1, 5 and 10 % in the table.
- Positive: model gene on the consensus list. Negative: model gene not on it (user-defined
  class). Listed genes with no model counterpart are out of scope, not FN.
- The 16 mtDNA genes are out of scope.
- The list was used during development: this is a development reference, not independent
  validation.

## Outputs

- This folder: TASK.md, scripts, `screen_<label>.manifest.json`, `summary.json`, the figure
  (PNG and PDF), REPORT.md.
- Raw per-KO tables: `$IYALI26_RESEARCH_ROOT/artifacts/tasks/confusion_matrices_csmleu_w29_20261014/`.
- A `model/STATE.md` entry.
