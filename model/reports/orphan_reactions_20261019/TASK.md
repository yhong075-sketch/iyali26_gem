# Orphan reactions in the newest candidate (2026-10-19)

## Authorization (user, this session)

1. "find the orphan reaction in our newest model" — read-only scan.
2. "Run the flux-evidence step (pFBA, FVA at 0% and 99.9%, single-reaction KO) on the orphans
   under SD-Leu + PO1f, as a read-only computation."
3. "yes. save the table. and analysis what should we do for these orphan reactions." — save the
   tables in this folder, add a `model/STATE.md` entry, and recommend next actions.

Out of scope: any model, GPR, chemistry, direction or medium change; BLAST or other sequence
searches; literature curation; commits and pushes.

## Starting state

- Worktree branch `claude/orphan-reaction-newest-model-b8ddc8` at `4d97b828ac2b0dc95d592060605af75069705c09`
  (= local `main`; `origin/main` `285eaee` adds Task D reports only, no model XML).
- Model executed: `model/candidates/E5_coq9_resp_pkg_r132_20261012.xml`,
  SHA-256 `a56f3649af2163ddedd36d1f53b830f2f550aed514aa4a28f6fb9ea577b58515`. Chosen as the newest
  candidate by last-commit date over all tracked candidate XMLs (2026-10-13, `e985dcd`); no other
  branch adds a newer XML. No reference model is designated (`model/candidates/README.md`), so
  this is the newest candidate, not a reference model.
- Condition: `model/reports/reference_pipeline_restore_20260909/research/state/media/sd_leu.csv`
  (`ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85`) with
  `.../strain_profiles/po1f_sd_leu.json` (`35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383`):
  whole medium replaced, R612 closed (ura3-302), uracil uptake override applied.
- `IYALI26_RESEARCH_ROOT` is unset; all outputs are small and kept here.

## Code identity

- Classification: `classify_gene_less` and `flag_passive_diffusion` from the
  `gem-orphan-reaction-curation` skill `kernel.py` (SHA-256 `19fe8cfa…`, not copied here).
- `scripts/flux_evidence.py` and `scripts/common.py`: copied unchanged from that skill.
- `scripts/orphan_evidence.py`, `scripts/loop_check.py`, `scripts/triage.py`: written this session.
- cobrapy 0.30.0, optlang 1.8.3, GLPK with presolve (Gurobi refused the model size on its
  restricted licence and was not used).

## Commands (from the worktree root; `$S` = report `scripts/`, `$D` = the state folder above)

```
python $S/flux_evidence.py --model model/candidates/E5_coq9_resp_pkg_r132_20261012.xml \
  --medium $D/media/sd_leu.csv --profile $D/strain_profiles/po1f_sd_leu.json --ids ids.csv --out orph
python $S/orphan_evidence.py --model <same> --flux flux_evidence.csv --ids ids.csv \
  --selection model/curation/selection/metadata_reaction_selection.json --out orphan_evidence.csv
python $S/loop_check.py $S <same model> $D/media/sd_leu.csv $D/strain_profiles/po1f_sd_leu.json > loop_check.txt
python $S/triage.py orphan_evidence.csv triage.csv
```

`ids.csv` = the 86 rows of `gene_less_classification.csv` with `category == metabolic`
(columns `id, name, spontaneous_note`).
