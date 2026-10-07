---
name: curation-change
description: Carry one user-authorized iYali26 model correction (GPR, reaction direction, chemistry/formula/charge, reaction selection) from scope to a new candidate model with whole-model verification, independent audit, REPORT.md and a model/STATE.md entry. Use when the user authorizes a specific model change. Not for exploratory review or for publishing (publishing needs a separate explicit instruction).
---

# Curation change

Stop and ask whenever the scope, the inputs or the baseline model is unclear. Make one change
per run. Run commands from `platform/`.

## 1. Scope (no writes yet)

- Quote the user's authorization.
- Name:
  - the baseline model actually used (path and SHA-256);
  - the target reactions, genes and metabolites;
  - the expected effect;
  - the acceptance criteria;
  - the compute budget and stop condition (root AGENTS.md §2).
- Create `model/reports/<topic>_<yyyymmdd>/TASK.md` with the above, and `before.json`
  recording the baseline state of every target.

## 2. Implement through the pipeline

1. Write the decision as a data file in `model/curation/<topic>/` with evidence fields: source
   URL, evidence tier and reviewer. The file name must be unique within `model/curation/`.
2. Add or extend one guarded `apply_*` function in `platform/scripts/gem_annotate/`. It must
   fail closed on stale rows or unexpected IDs and must never touch anything outside the
   targets. Load the file with `MODEL.curation_file(name)`.
3. Add a focused test in `platform/tests/` that runs on an in-memory model.
4. Call the function from `main.py` at the right stage.
5. Run `python -m pytest` and compare with the previous run test by test.

## 3. Build to a new file

```bash
python -m scripts.gem_annotate --research-root "$IYALI26_RESEARCH_ROOT" \
  --offline --no-solve --coq9-curation metadata \
  --output-model ../model/candidates/<new_name>.xml
```

Record the command, the output SHA-256 and the `.build.json`.

## 4. Verify the whole model

- Compare the baseline and the new model across every reaction, metabolite, gene, GPR, bound,
  formula, charge, annotation and note.
- **Only the intended targets may differ.** Any other difference means: stop and report it.
- Check mass and charge balance for every touched reaction.
- Run only the authorized solves (WT growth, target controls) and record their status. A
  non-optimal solve is not biological evidence.

## 5. Independent audit

- Give a fresh read-only reviewer only TASK.md, the diff, the logs and the deliverables.
- It returns one verdict per claim: supported, partially supported, unsupported or unverified.
- Fix the cause or report the disagreement; never argue a verdict away.

## 6. Write up

- **`REPORT.md`:** explain
  - what changed and why;
  - the mechanism and the evidence, each tagged verified / inferred / recalled;
  - results, limits, what was not tested, and next gates.

  Put large outputs in `$IYALI26_RESEARCH_ROOT/artifacts/tasks/<same folder name>/` and link
  to them.
- **`model/STATE.md`:** add a newest-first dated entry with the authorization, outcome,
  verification level, links and limits.
- **Publishing:** stop here. Publishing needs a separate explicit instruction.
