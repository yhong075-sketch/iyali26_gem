# Essentiality false-negative review workflow

This repository uses a human-gated, literature-driven workflow for SD-Leu
essentiality false negatives. These instructions apply whenever the user's
request matches one of the four Chinese commands below.

## `审查下一批 essentiality FN`

1. Resolve the external workspace from `IYALI26_RESEARCH_ROOT`; fail closed if
   it is unset or incomplete. Run a fresh validation from the current
   `model.xml` with the positive-only experimental list,
   `$IYALI26_RESEARCH_ROOT/state/media/sd_leu.csv`, cutoffs 1%, 5%, 10%, 15%,
   `--diagnose --prepare-agent-cases --batch-size 3`. Never reuse an old result
   whose model, experimental, or medium SHA-256 differs.
2. Read
   `$IYALI26_RESEARCH_ROOT/artifacts/results/essentiality/essentiality_agent_batch.json`.
   Delegate each of
   its three independent packets to a separate
   `yarrowia-essentiality-literature-reviewer` with no conversation history.
   Before delegation, move only those three ledger rows from `queued` to
   `researching` using the guarded state helper. Do not give one reviewer
   multiple cases. Wait for all three results.
3. Send all fresh packets and reviewer results together to one
   `essentiality-evidence-skeptic`. It runs after, not in parallel with, the
   three reviewers.
4. Validate the returned structures with
   `scripts.gem_annotate.essentiality_evidence`, then update the durable dossier
   and move `researching -> reviewed`. A skeptic pass may then move a supported
   candidate to `awaiting_human`; it must never move it to `accepted`.
5. Write the diagnosis and screen-test update to
   `/Users/david/Desktop/Lab/Ian wheeldon/code/Genome-wide/weekly_reports`.
   Do not modify the model during this command.

All research agents are read-only. Their remit is direct experimental evidence
in *Yarrowia lipolytica*. UniProt/KEGG are identity cross-checks; related yeasts
and databases alone cannot make a patch acceptable. Recall is never evidence.

## `接受 EGC-xxxxxxxxxxxx`

Treat acceptance as valid only when the current user message explicitly names
that exact case. Run the human-decision recorder, verify dossier, skeptic pass,
input SHAs and current target fingerprint, and then delegate that one case to
`essentiality-patch-builder`. The builder may change only curated patch data,
pipeline functions and tests. It must never edit `model.xml`; the final model is
rebuilt from `data/iyali26.xml` by the normal pipeline.

After rebuild, report model changes and final SHA in
`$IYALI26_RESEARCH_ROOT/artifacts/weekly_briefing`, and
write essentiality diagnosis/screen regression to the external
`Genome-wide/weekly_reports` directory. Report 1%, 5%, 10%, and 15%, with 10% as
the primary cutoff. Do not report FP, TN, precision, accuracy or MCC for the
positive-only reference.

## `拒绝 EGC-xxxxxxxxxxxx`

Record `human_decision=rejected`, the decision time and `status=rejected` in
the durable ledger/dossier. Do not call the patch builder and do not change the
model.

## `延后 EGC-xxxxxxxxxxxx`

Record `human_decision=deferred` and `status=needs_more_evidence`. Do not call
the patch builder and do not change the model.

## Invariants

- State order is `detected -> queued -> researching -> reviewed ->
  awaiting_human/needs_more_evidence/rejected -> accepted -> implemented ->
  regression_passed`.
- Only an explicit human `接受 EGC-...` command may create `status=accepted`,
  `approved_by=human_user`, and `approved_at`.
- A changed target fingerprint invalidates simulation evidence. Literature may
  be reused, but diagnosis and adversarial review must be rerun.
- Do not tune SD-Leu uptake, close a bypass, alter a GPR, or add biomass demand
  merely to improve recall.
- A new `supported_patch_candidate` must have a current-SHA
  `chemistry_review` showing balanced reaction chemistry and a current-SHA
  `identity_review` with `status=verified`. If a connected-component
  microspecies audit is present, it must report `ready_for_activation=true`.
  Otherwise the case remains `needs_more_evidence` and cannot enter human
  approval.
- The legacy `EG-GPR-001` remains active under schema v1 while its evidence is
  backfilled. Every new patch must use schema v2 gates.

# M00-M11 metabolic module agent workflow

These instructions apply when the user explicitly asks to use the metabolic
module workflow, requests a whole-model module audit, or assigns work to one or
more modules `M00` through `M11`. The machine-readable source of truth is
`config/module-agent-workflow.json`; validate it before delegation with:

```bash
python -m scripts.gem_annotate.module_workflow
```

## Orchestration

1. The primary agent is the global orchestrator. Record the validated workflow
   SHA-256 and current `model.xml` SHA-256 in every module packet. Fail closed
   when either fingerprint changes during a run.
2. Derive execution waves from the configured dependency DAG. Do not maintain
   or follow a separate handwritten wave list. With the five-agent runtime
   budget, run only one complete module team at a time.
3. Delegate one bounded module to `metabolic-module-lead`. The lead owns
   scientific reasoning, module interfaces, acceptance criteria, and
   cross-module escalation. It does not edit files.
4. The lead delegates to `metabolic-workflow-coordinator`. The coordinator
   resolves explicit target reaction/metabolite/gene IDs, assigns exclusive
   file ownership, sequences work, and consolidates results. It does not edit
   files.
5. For an authorized change, the coordinator delegates the bounded work packet
   to `metabolic-module-implementer`. Only after implementation and targeted
   tests finish may it delegate the resulting diff and evidence packet to
   `metabolic-module-reviewer`. The implementer never reviews its own work.
6. The coordinator reports to the lead, and the lead reports to the primary
   agent using every field in the configured report contract. Cross-module
   findings use every field in the configured handoff contract and pass through
   the primary agent; do not create untracked all-to-all agent conversations.

## Module workflow invariants

- The `M00`-`M11` classification is currently conceptual. SBML subsystem groups
  are not authoritative module membership. Default to a read-only inventory
  until an explicit work packet identifies the target model entities and files.
- A request to review, explain, or diagnose does not authorize implementation.
  Activate the write-capable implementer only for an explicit user change
  request with a bounded scope, baseline SHA, file ownership, and acceptance
  gates.
- Never edit `model.xml` directly. Change curated inputs, pipeline functions,
  and tests, then rebuild from `data/iyali26.xml` through
  `scripts.gem_annotate` when a model rebuild is authorized.
- At most one agent owns a writable file at a time. Preserve unrelated changes
  and stop if ownership cannot be made exclusive.
- Treat module arrows as dependency and handoff edges, not as proof of reaction
  irreversibility. In particular, `M09 -> M04` is a real prosthetic-group
  dependency whose coupling may be absent from model stoichiometry.
- Any named gene must be reported with its model identifier, established name
  or symbol, concise protein function, and evidence status. Database annotation
  alone is not experimental proof.
- The essentiality false-negative commands above always take precedence. The
  generic module workflow may assist diagnosis but may not accept a case,
  bypass its dedicated reviewers/skeptic, or modify the model on its behalf.
