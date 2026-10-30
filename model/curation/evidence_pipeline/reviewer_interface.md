# The existing AI reviewer: interface as found (M0, 2026-10-29)

This document answers spec §6 ("locate it, read it, and write `docs/reviewer_interface.md`") for the
evidence pipeline. It describes what exists. It wires nothing and changes no prompt; the user confirms this
document before anything is connected (spec §6).

Everything below was read this session from `origin/main` @ `f212874` (2026-10-29), from branch
`codex/module-agent-workflow` @ `35c959b` (2026-07-23), and from the research workspace
(`$IYALI26_RESEARCH_ROOT`). Statements are **verified** (file read this session) unless marked
**inferred**.

## 1. Which reviewer

The repo has three review mechanisms. Only the first is the "adversarial, multi-agent review workflow" the
spec describes.

| # | Reviewer | What it reviews | Where | Runs on |
| --- | --- | --- | --- | --- |
| A | **Essentiality FN loop**: one literature reviewer per case, then one skeptic per batch of three, then a human gate | Biological claims behind an essentiality false negative (GPR, bypass, inactive reaction) | Agents: `.codex/agents/*.toml` on `codex/module-agent-workflow` (not on `main`). Ledger code: `platform/scripts/gem_annotate/essentiality_evidence.py` (on `main`). Case preparation: `validate_essential_genes --prepare-agent-cases` (only on `codex/module-agent-workflow`). Data: `$IYALI26_RESEARCH_ROOT/state/essentiality/repository/` | Codex |
| B | `code-reviewer` subagent | Code diffs against the `platform/AGENTS.md` rules | `.claude/agents/code-reviewer.md` (on `main`; not yet in this worktree, which is 63 commits behind) | Claude Code |
| C | Independent audit in the `curation-change` skill | One model change: a verdict per claim in TASK.md/REPORT.md (supported / partially / unsupported / unverified) | `.agents/skills/curation-change/SKILL.md` §5; 134 `AUDIT.md` files under `model/reports/` | Fresh read-only agent briefed by hand |

`platform/gem_review/lit` (Europe PMC → quote-verified claims → logic chains) is an extraction layer, not a
reviewer: one agent, no adversarial pass. It is covered in the M0 report as the base for the spec's
L0–L2.

The rest of this document is about **A**. B and C are noted at the end.

## 2. Reviewer A: entry point and invocation

The procedure is written in the `AGENTS.md` of `codex/module-agent-workflow`. The current root `AGENTS.md`
§6 keeps the rules and says the loop is paused on this line. It is triggered by four user commands in
Chinese:

| Command | Effect |
| --- | --- |
| `审查下一批 essentiality FN` ("review the next batch of FN") | Fresh screen → three case packets → three literature reviewers → one skeptic → validated import → `reviewed` (a supported case may go to `awaiting_human`) |
| `接受 EGC-…` ("accept") | Human decision recorded; that one case goes to `essentiality-patch-builder` |
| `拒绝 EGC-…` ("reject") | `rejected`; no builder |
| `延后 EGC-…` ("defer") | `deferred` → `needs_more_evidence`; no builder |

Steps of `审查下一批 essentiality FN` (as written on the Codex branch):

1. Fresh screen with case preparation (old layout, old contract):
   `python -m scripts.gem_annotate.validate_essential_genes --model model.xml --experimental …/consensus_essential_genes.csv --media …/sd_leu.csv --positive-only --primary-cutoff 0.10 --growth-cutoffs 0.01,0.05,0.10,0.15 --diagnose --prepare-agent-cases --batch-size 3 [--case-category …] --solver gurobi`.
   `build_agent_cases()` groups FN diagnostics by (category, reaction IDs, GPR) into deterministic packets.
   `prepare_agent_case_files()` writes `essentiality_agent_cases.jsonl` and
   `essentiality_agent_batch.json` (three cases).
2. `python -m scripts.gem_annotate.essentiality_evidence start-batch --batch-file <batch.json>` moves the
   three ledger rows `queued → researching`.
3. Each packet goes to its own `yarrowia-essentiality-literature-reviewer`, with no conversation history and
   one case per reviewer.
4. All packets and all three results go together to one `essentiality-evidence-skeptic`, after the
   reviewers finish.
5. `essentiality_evidence import-batch-review --batch-file … --reviewer-results … --skeptic-result …`
   validates every structure, writes the dossier and moves the cases `researching → reviewed`.
6. Optional: `import-identity-review` and `import-chemistry-review` per case; `accept | reject | defer
   <case_id>` records the human decision.

How the agents are spawned: by the Codex orchestrator as Codex custom agents from `.codex/agents/` (inferred
from the file location and the `sandbox_mode` keys; no spawn code exists in the repo). No Claude Code agent
definition exists for the literature reviewer or the skeptic. A Claude session could only run them as
`general-purpose` agents briefed with the TOML's `developer_instructions`, which would be a different
runtime from the one that produced the stored results.

## 3. Input format: the case packet

Schema version `2.0` (`build_agent_cases`, Codex branch):

```
case_id ("EGC-" + 12 hex, stable over category, genes, reactions), category, priority, ranking_reason,
gene_ids, reaction_ids, model_sha256, experimental_sha256, media_sha256,
target_fingerprint, chemistry_fingerprint,
claim_under_review, falsifiable_questions        # fixed text per category (CASE_CLAIMS / CASE_QUESTIONS)
threshold_results {gene: {cutoff_label: bool}},
model_context {model_id, wt_growth, primary_cutoff_fraction_of_wt, medium{path, sha256, leucine…},
               genes[{gene_id, name, annotations, gpr_partners, reaction_ids}],
               reactions[{id, equation, gpr, gpr_gene_ids, compartments, ec_codes, annotations, …}],
               diagnostics[ko_growth, ko_growth_ratio, verified_single_bypasses, isozyme_interpretation, …],
               verified_bypasses[reaction contexts]},
case_packet_sha256
```

Categories in the ledger: `inactive_reaction` 136, `metabolic_bypass` 66, `isozyme_redundancy` 35,
`nutrient_bypass` 3.

The packet's unit is **one FN case (genes + reactions)**. It is not a (subject, predicate) claim set.

## 4. Output formats

**Literature reviewer** (one JSON object per case):

```
case_id, claim_under_review, primary_sources[], identity_crosschecks[], contradictions[],
verdict ∈ {supported_patch_candidate, retain_current_model, experimental_conflict,
           outside_metabolic_scope, needs_more_evidence},
proposed_operation {operation, target_id, value, rationale}, confidence, unresolved_questions[]
```

The stored dossier also carries `search_audit` (databases, queries, inclusion/exclusion criteria, screened
and excluded sources, `direct_evidence_found`), `reasoning`, and the binding hashes
(`case_packet_sha256`, `reviewer_result_sha256`, model, experimental and medium SHA, target and chemistry
fingerprints). Each `primary_sources` item has:

```
source_id, url, doi, title, year, species, strain, culture_conditions, source_type, evidence_type,
stance ∈ {supports, contradicts, contextual}, claim, location (figure/table/section), evidence_tags,
genes, reactions, pathways, methods, result,
condition_match ∈ {exact, partial, mismatch, unknown}, condition_mismatch_reason, relevance, confidence
```

`evidence_type` is free text: 15 different values in the 196 stored sources, for example
`direct_experiment` 151 and `primary_experiment` 19.

**Skeptic** (stored under the dossier key `adversarial_review`):

```
case_id, status="complete", verdict ∈ {pass, block, downgrade}, findings[], unresolved_contradictions[],
corrected_verdict, confidence
```

**Validation and gates** (`essentiality_evidence.py`, on `main`): the validator checks the required fields,
the stance and condition enums, provenance binding and the state transitions
(`detected → queued → researching → reviewed → awaiting_human | needs_more_evidence | rejected → accepted →
implemented → regression_passed`). Only `accept` by a human creates `accepted` with
`approved_by=human_user`. A changed target fingerprint invalidates simulation evidence.

**Patch builder** (`essentiality-patch-builder`): a writer, not a reviewer. It changes curated patch data,
pipeline code and tests for one accepted case.

## 5. Current state of the stored output

Ledger `curation_cases.csv` (240 rows) and 240 dossiers, read this session:

- Status: queued 153, reviewed 85, researching 1, implemented 1. Human decisions: 1 `legacy_accepted`, none
  recorded under schema v2.
- Literature verdicts: needs_more_evidence 86, experimental_conflict 2, supported_patch_candidate 1.
- Skeptic (verdict → corrected verdict): block → needs_more_evidence 77, pass → needs_more_evidence 4,
  downgrade → needs_more_evidence 4.
- Every case is bound to an old model SHA (`b3f60933` 205, `39f4cae1` 34, `0eb99eaf` 1) and to the SD-Leu
  medium. None is bound to any current candidate (newest on `main`: `d8d1bac7`). By the reviewer's own rule
  their simulation evidence is stale; the literature parts may be reused.

## 6. What it is good at

- **Adversarial structure that already exists:** isolated per-case reviewers, then a separate skeptic that
  must try to falsify; the skeptic blocked 77 of 85 reviewed cases.
- **Strict evidence scope:** only direct *Y. lipolytica* experiments can support a patch. Related yeasts and
  databases are background; UniProt/KEGG are identity checks only; recall is never evidence.
- **Checks the spec also wants:** YALI0/YALI1 normalization and paralog confusion; circular support (a
  database repeating the model's own annotation); strain, temperature, medium, carbon source, compartment and
  assay endpoint; and whether the cited figure or table supports the exact operation.
- **Hard provenance gates in code:** packet and result hashes, input-SHA binding, fingerprint invalidation,
  a guarded state machine and a human-only acceptance gate.

## 7. What it is not asked to do (gaps against spec §6)

| Spec role | Closest existing piece | Gap |
| --- | --- | --- |
| **R1 Adjudicator** (`ConflictPacket` → `ReviewVerdict`) | Literature reviewer on one case packet | It is not given claims, passages, pair relations or a rule-engine trace, so it cannot agree or disagree with a rule. Its verdict vocabulary differs from `side_positive / side_negative / unresolved / scope_error / extraction_error`, and it has no `disputed_rule` field. It searches for evidence itself instead of judging a fixed evidence set. |
| **R2 Red team of extraction** (`ClaimChallenge`) | Skeptic's "cited figure/table supports the exact operation" check | Nothing re-reads a stored passage to challenge polarity, assertion strength, claim level or scope, and there is no re-extraction trigger. |
| **R3 Decision gate** (`DecisionReview`) | Skeptic pass/block/downgrade, before `awaiting_human` | It gives a verdict and findings, not "the strongest case for and against plus a recommendation". It is scoped to essentiality FN cases, not to arbitrary L6 actions (relocate, change bounds, reassign subsystem). |
| `reviewer_disagreements` per rule | none | No rule IDs exist to disagree with. |

Other limits:
- **Scope:** essentiality false negatives only. It is never asked about reactions with no FN gene, about
  localization, directionality or cofactors outside an FN case, or about subsystems.
- **Old contract baked into the prompts:** SD-Leu medium, 10 % primary cutoff, positive-only labels,
  `data/iyali26.xml` / `model.xml` paths. The current contract (`model/expected/benchmark.md`) is the Task A
  medium + PO1f, a 15 % primary cutoff, user-defined negatives and the `model/` layout. Changing the prompts
  needs the user's permission (spec §6).
- **Not runnable from `main` as it stands:** the case-preparation half and the agent TOMLs are only on
  `codex/module-agent-workflow` (July layout). `essentiality_evidence.py` on `main` points its defaults at
  the Task A medium and `MODEL.canonical_model` (`candidates/legacy/model.xml`), which is not a designated
  reference.
- **Free-text `evidence_type`:** it would need a deterministic mapping to the spec's `claim_level`, or a
  re-extraction, before it could feed L2.
- **Runtime:** Codex agents. Their stored outputs were produced by Codex and recorded with no model ID or
  prompt hash (the spec requires both for every LLM call).

## 8. Reviewers B and C, briefly

- **B `code-reviewer`:** input is the task statement, a diff range and the tests run. Output is ranked
  findings `file:line · severity · what · failure scenario · verified/inferred`. Read-only; no builds,
  solves or network. Useful for reviewing the pipeline's own code at every milestone, not for scientific
  adjudication.
- **C curation-change audit:** input is TASK.md, the diff, the logs and the deliverables. Output is one
  verdict per claim (supported / partially / unsupported / unverified). This is the existing gate after a
  model change and the natural place for the spec's "apply approved patch, then validate" step.

## 9. Questions for the user before anything is wired

1. Is reviewer A the one you mean, and which of its outputs may gate decisions: the skeptic's
   `pass/block/downgrade`, the literature reviewer's verdict, or neither (advisory only)?
2. Should the R1–R3 contracts reuse A's agents with their prompts unchanged (adapting the packets around
   them), get new prompts (your permission needed), or be new Claude agents that leave A untouched?
3. Should A's 85 stored dossiers (196 sources) be imported as L0/L2 seed evidence? This needs a mapping of
   the free-text `evidence_type` to `claim_level` and re-verification of every quote. The dossiers store
   paraphrased `claim`/`result` text, not verified verbatim passages, so most would enter L2 only after
   re-extraction.
4. Which runtime should the reviewer run on, Codex or Claude? The spec requires the model ID and prompt hash
   to be recorded for every call.
