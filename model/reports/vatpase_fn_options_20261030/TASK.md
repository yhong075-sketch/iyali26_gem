# TASK: options for the V-ATPase false negatives (investigation only)

## Authorization

- 2026-10-30, the user: "how can we solve the VTPase problem?" This was a question, so the deliverable is an
  assessment. Ultracode was on, so the work ran as a read-only multi-agent workflow.
- 2026-10-30, the user: "Save as report". This covers saving the assessment as this folder plus a `model/STATE.md`
  entry. It does not cover a model change.

## Scope

- **In scope:** why the 11 V-ATPase genes on R794/R795 are false negatives, which fixes could work, and which are
  backed by evidence. Static model reading, Boolean GPR evaluation, report and literature reading.
- **Out of scope:** any model, GPR, bound, medium, biomass or label change; any LP, FBA or FVA; builds; commits.

## Starting state

| Item | Value |
| --- | --- |
| Branch / HEAD | `claude/essential-gene-prediction-e76be3` @ `ff7fc62` (the analysed model is on `origin/main`, not on this branch) |
| Model analysed | `E5_…_r171_atp_synthase_gpr_r673_prolyl_split_notes_v3_20261029.xml` from `origin/main` (`059316ef34417c88aadba61cfd5606134d27879a842b3835d2cc7df917a9248d`), the newest candidate on main; not a designated reference |
| Reference model | none designated (`model/candidates/README.md`) |
| Experimental positives | `consensus_essential_genes.csv` (`1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b`), 1,612 genes |
| Raw screen calls | `cas9_cas12a_fitness.csv` (`97ef559651cd99ce63144ffe08ae64e983ecaec60d568f859da8926c14c7d9ff`) |
| FN labels | `ko_6af03c40` sweep (Task A + PO1f, Gurobi 1e-9). `6af03c40` and `059316ef` differ only in notes (recalled from the 2026-10-29 STATE entries of PR #37; not re-solved here) |

## Compute

No LP was solved. Static and Boolean checks only, re-run by the lead on `059316ef` with
[scripts/static_checks.py](scripts/static_checks.py).

## Stop condition

The report is written and the STATE entry added. Nothing is committed.
