# iYali26 GEM — agent instructions

Shared by Codex and Claude Code (both load this file). Keep it short: it is always in
context, and Codex stops reading at 32 KiB of combined AGENTS.md. Procedures live in skills
under `.agents/skills/` (also visible as `.claude/skills/`). Code-specific rules are in
`platform/AGENTS.md`. Conversation follows the user's language; code, comments, commits and
data files are English.

## 0. Before any task

- **Current line of work:** resolving inactive reactions (for example CoQ9 synthesis) on
  branch `codex/r989-gpr-main-worktree` and its successors. The lipid-unlump line is separate
  and owned by a colleague.
- **Identify separately:**
  - the current branch and HEAD;
  - the model you will actually execute (path plus SHA-256);
  - the reference model, if any;
  - the candidates.

  No reference model is designated yet (`model/candidates/README.md`). Words like "main",
  "canonical" and "latest", and file or folder names, do not define the scientific baseline.
- **Research workspace:** resolve it from `--research-root` or `IYALI26_RESEARCH_ROOT`. Fail
  closed if it is unset or incomplete.
- **Read first:** the relevant entries of the state log (`model/STATE.md`) and the evaluation
  rules (`model/expected/benchmark.md`). Do not redo full audits or start a new state ledger.
- **Instruction files:** check every file that applies (global, ancestor, this directory,
  subdirectory, overrides). Another worktree's AGENTS.md does not apply here.

## 1. Where things are

| What | Where |
| --- | --- |
| Immutable build input | `model/source/iyali26.xml` |
| Curation decisions (data) | `model/curation/<topic>/` |
| Candidate models + provenance | `model/candidates/` |
| Per-task reports | `model/reports/<topic>_<yyyymmdd>/` (TASK.md, REPORT.md, audits, small tables) |
| Decision log | `model/STATE.md` (newest first) |
| Evaluation rules | `model/expected/benchmark.md` |
| Builder (only writer of model XML) | `platform/scripts/gem_annotate/`, run `python -m scripts.gem_annotate` from `platform/` |
| Tools and tests | `platform/tools/`, `platform/tests/` |
| Reference data, ledgers, large outputs | `$IYALI26_RESEARCH_ROOT` (`reference/`, `state/`, `artifacts/tasks/<topic>_<date>/`) |

## 2. Authorization (needs the user's explicit, current instruction)

Investigation, documentation and computation are different tasks. Each limit applies to its
round only; a past limit is not a permanent ban.

- **Computation:** needs named inputs, scope, budget, outputs and a stop condition.
  Exception: HPCC AlphaFold jobs for protein-identity questions are pre-authorized. Check the
  environment first, keep job IDs, logs and results, and label the output
  "AlphaFold prediction".
- **Science changes** (model, GPR, chemistry, medium, experimental labels, benchmark or test
  baselines):
  - They need scoped authorization plus evidence.
  - Never make one to improve a match rate.
  - Make them only through curated data and the build. Keep the raw inputs and old results;
    never overwrite a formal model or evidence.
- **Repository and external actions:** git commit, push, branch switch or merge;
  moving/deleting/cleaning files; cluster jobs other than AlphaFold; external messages.
  - Publish only to the branch the user names; `main` only if named.
- **Others' work:** preserve the user's and other sessions' uncommitted work. Report unknown
  changes; never revert them.
- **Ownership:** lipid-unlump (`platform/tools/lipid/`, `model/curation/lipid_unlump/`,
  `model/candidates/lipid_unlump/`) is owned by a colleague; review and recommend only.
  Owners of inactive/isozyme work come from `model/STATE.md`. If unknown, leave them unknown.
- **Gates that stay in force:** the CoQ9 and lipid-unlump approval gates. Passing tests,
  updating docs or registering a provisional reference never accepts a case or authorizes a
  compute matrix.
- **Secrets:** never print secrets or credentials.

## 3. How a model change flows

Skill: `curation-change`.

1. User authorizes scope →
2. `model/reports/<topic>_<date>/TASK.md` records the scope and the starting state →
3. Curated data file plus a guarded `apply_*` step and a focused test →
4. Build to a **new** file in `model/candidates/` →
5. Whole-model diff showing only the intended changes →
6. Independent read-only audit →
7. REPORT.md →
8. `model/STATE.md` entry →
9. Publish only on an explicit push instruction.

- Never hand-edit any model XML.
- Every claim that a model was updated gives the baseline SHA, output SHA, changed entities
  and the exact command.
- Deleting a reaction or metabolite needs evidence that it does not exist in the organism.
  A score, memote or recall gain is not evidence.
- A charge or formula edit propagates to every reaction that uses the metabolite. Re-check
  balance model-wide.
- Large outputs (flux dumps, run folders, figures over about 1 MB) go to
  `$IYALI26_RESEARCH_ROOT/artifacts/tasks/<same folder name>/`, linked from REPORT.md.

## 4. Evidence and reporting

- Keep three claims apart: "a report claims it is done", "the deliverable is verified",
  "reproduced this time". Record when you verified, what evidence you used, and the scope
  checked. Never fill in approvals or owners from history you did not read.
- Record full SHAs of inputs, code identity, medium, runtime strain profile, parameters and
  results. Keep any historical dirty state; matching a few source files is not rebuilding
  the historical environment.
- Unknown stays unknown. Byte-identical, annotation-identical, same-optimization-problem and
  biologically-equivalent are four different claims.
- **Scoring** follows `model/expected/benchmark.md`:
  - 15% KO/WT is the primary cutoff; also report 1%, 5% and 10%.
  - Model genes not on the experimental list count as non-essential (a user-defined negative
    class), so TP/FN/FP/TN are reported. Always state that definition.
  - Report positive coverage and out-of-scope positives separately.
  - Data used to build, diagnose or curate is not independent validation.
- Keep raw calls, IDs and sources. Keep format normalization separate from cross-version gene
  mapping. Never guess missing calls or swap labels.
- Name every gene with:
  - its systematic ID;
  - its verified name, or say none is verified;
  - a short protein function;
  - the evidence level.

  A model or GPR assignment is not experimental validation.
- **AlphaFold**, for weak identity evidence:
  - Reuse an existing identical model when one exists.
  - Record the sequence source, version and SHA, the tool and version, the date, pLDDT, PAE
    and limits.
  - Label the result "AlphaFold prediction" or "function candidate based on AlphaFold
    prediction". A prediction never confirms identity, catalysis, specificity or compartment.
- Keep historical numbers, thresholds and exception handling as they were. A solver failure,
  missing value or non-finite value is not evidence of biological death. Keep raw values,
  normalized values and review verdicts apart.
- **Reports** explain:
  - the code logic, the results and what they mean;
  - whether the math or algorithm could be improved, and how to test it;
  - which statements are provable, statically verified, already observed in runs, or
    untested.

  Keep hash tables out of the report body unless they change a conclusion.

## 5. GPR, direction and compartment review

Scope the checks to the question being asked and reuse existing evidence. Do not expand into
a full-model audit or a new compute matrix.

- **GPR propagation:** trace gene state → Boolean GPR → reaction bounds → growth.
  - An OR member needs evidence that it catalyses the step on its own.
  - An AND member needs evidence of a complex or co-requirement.
  - Look for duplicate reactions that bypass a knockout.
- **Direction and bypasses:** trace them by stoichiometry, cofactors and compartments.
  - A direction change needs enzymology, thermodynamics and conditions; never a reaction
    name, another model's bounds or a hit rate.
  - Separate net-synthesis compensation from net-zero cycles.
- **Compartment supply:**
  - List every producer, consumer and transporter per compartment.
  - Direct cofactor transport needs a known mechanism. Never treat reducing-power transfer
    through other reactions as direct cofactor transport.
  - A missing producer or transporter is not proof of essentiality.
  - Match reactions by substrate, product and compartment, never by name or gene ID.
- **Mechanism vs biology:** answer them separately.
  - Use targeted controls: WT, target KO, bypass closed, and combined.
  - Keep the model, medium, objective and threshold fixed.
  - One optimum does not prove the route is unique.
- **Candidate edits:** record the original structure, the mechanism, evidence for and
  against, the proposed change and its acceptance criteria.
  - After an approved change, re-check the target and the known-correct results.
  - Never copy an AND rule, delete a bypass or accept a supply gap just because an external
    model hits more positives.
  - Methods derived here are this project's methods, not the original authors'.

## 6. Essentiality false-negative loop (paused on this line)

The commands `审查下一批 essentiality FN`, `接受 EGC-xxxxxxxxxxxx`, `拒绝 EGC-xxxxxxxxxxxx` and
`延后 EGC-xxxxxxxxxxxx` need tooling that exists only on `main`: `validate_essential_genes
--prepare-agent-cases`, the case ledger commands and the reviewer agents. On this line,
reply that the loop is unavailable until that tooling is brought over. If it is brought back,
these rules apply:

- **Review:** one fresh screen per batch; never reuse a result whose model, experimental list
  or medium SHA differs.
  - Batches have three cases, moved `queued → researching` with the guarded helper.
  - Each case goes to its own read-only literature reviewer with no conversation history.
  - Then one skeptic reviews all three.
  - Results are validated by the evidence module.
  - A skeptic pass may move a supported case to `awaiting_human`, never to `accepted`.
- **Evidence:** reviewers count only direct experimental evidence in *Y. lipolytica*.
  UniProt/KEGG are identity cross-checks; related yeasts and databases cannot make a patch
  acceptable; recall is never evidence.
- **`接受`:** valid only when the user's current message names that exact case. It records
  the decision, re-checks the dossier, skeptic pass, input SHAs and target fingerprint, and
  hands that one case to the patch builder. The builder may change only curated patch data,
  pipeline code and tests, never a model XML, then rebuilds and reports the regression at all
  four cutoffs.
- **`拒绝` / `延后`:**
  - `拒绝` records `rejected`.
  - `延后` records `deferred` → `needs_more_evidence`.
  - Neither calls the builder or changes the model.
- **State order:** `detected → queued → researching → reviewed → awaiting_human |
  needs_more_evidence | rejected → accepted → implemented → regression_passed`.
  - Only an explicit human `接受` creates `accepted`, `approved_by=human_user` and
    `approved_at`.
  - A changed target fingerprint invalidates simulation evidence; literature may be reused.
- **No recall tuning:** never tune SD-Leu uptake, close a bypass, alter a GPR or add
  biomass demand merely to improve recall.
- **New candidates:** a new `supported_patch_candidate` needs a current-SHA
  `chemistry_review` (balanced) and `identity_review` (`verified`). If a connected-component
  microspecies audit is present, it must report `ready_for_activation=true`; otherwise the
  case stays `needs_more_evidence`.
- **Legacy patch:** `EG-GPR-001` remains active under schema v1. New patches use schema v2.

## 7. Working rules

- **Gaps:** classify each gap as "blocks this task", "limits a claim" or "affects another
  line". Pause only what is affected.
  - On an input or fingerprint change, an evidence conflict, or scope overrun: keep the
    facts and name the gap.
  - Never auto-replace inputs, retry, tune or relax a gate.
- **Wrap-up:** after each task, add a dated `model/STATE.md` entry: the authorization, the
  results with their verification level, links and limits.
- **Weekly report:** the Sunday advisor report reuses `model/STATE.md` and its evidence.
  Describe the actual verification level; never infer completion from file names, finished
  jobs or small samples.
