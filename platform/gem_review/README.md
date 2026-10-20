# gem_review — literature review for iYali26 gene–reaction pairs

Read-only review tooling: it searches Europe PMC, extracts verbatim, validated quotes with an
LLM, and writes `review.xlsx`. It never edits a model file; anything it suggests still goes
through the normal curation flow (`AGENTS.md` §3).

All commands run from `platform/`:

```
uv run python -m gem_review.lit <command> [options]
```

Every command needs the research workspace: pass `--research-root <path>` or export
`IYALI26_RESEARCH_ROOT`.

## 0. Choose the AI model (once)

```
uv run python -m gem_review.lit setup
```

A short terminal wizard: pick a provider (Anthropic, OpenAI, Google Gemini, Groq, Mistral,
DeepSeek, OpenRouter, xAI, or a local Ollama), accept or edit the model name, paste the API key
(input hidden), optionally test it with one tiny call, and it is saved to `platform/.env`
(git-ignored, owner-only permissions) as `GEM_REVIEW_LLM=provider:model` plus the provider's key
variable. Every later command reads that file; `--llm provider:model` overrides it for one run.
Providers marked "not installed" need one command first, which the wizard prints, e.g.
`uv add 'pydantic-ai-slim[groq]'`.

The suggested model names are defaults to edit, not recommendations: check your provider's
current model list.

## 1. Retrieve literature (no LLM, free)

```
uv run python -m gem_review.lit retrieve --model ../model/candidates/<model>.xml --tier 1 --task <name>
uv run python -m gem_review.lit retrieve --model ../model/candidates/<model>.xml --reactions R348,R1771 --task <name>
```

Selects reactions by tier (tier 1 = lipid subsystems in
`model/curation/gem_review/tier1_lipid_subsystems.csv`) or by ID — never the whole model. For
each gene–reaction pair it queries Europe PMC with the gene's names (YALI1/YALI0 tags, UniProt
symbols) and the reaction's EC/name, fetches open-access full texts, and stores abstracts and
matching paragraphs. Output: `model/reports/<task>/hit_counts.csv` and per-pair JSON under
`$IYALI26_RESEARCH_ROOT/artifacts/tasks/<task>/out/lit/retrieval/`. All HTTP responses are
cached, so re-runs are free.

## 2. Extract validated claims (LLM, costs tokens)

```
uv run python -m gem_review.lit extract --task <name> [--pairs GENE:RXN,...] [--limit N]
```

Only pairs with hits are sent. Each claim must cite a retrieved PMID and a quote (≤ 25 words)
that appears verbatim in the retrieved text; otherwise the model is asked again (twice), then
the pair is recorded `EXTRACTION_FAILED`. The evidence level (0 none … 4 biochemical) is
computed from the validated claims; evidence from other organisms counts at most as homology.
Token usage per pair is logged to `out/lit/usage.jsonl`.

## 3. Build the spreadsheet

```
uv run python -m gem_review.lit report --task <name>
```

`model/reports/<task>/review.xlsx` with sheets `literature` (one row per pair),
`review_queue` (pairs needing a human: `contradicted`, `claims_disagree`, `low_evidence`,
`no_hits`) and `logic_chains`.

## 4. Ask a question as a logic chain

```
uv run python -m gem_review.lit chain --task <name> --pair YALI1D32754g:R100
uv run python -m gem_review.lit chain --task <name> --question "Does Dga1p acylate diacylglycerol in the ER?" \
    --gene YALI1E38810g --aliases DGA1,YALI0E32769g --terms "2.3.1.20,diacylglycerol acyltransferase" [--retrieve-only]
```

The answer is an ordered chain of verbatim quote steps and inference steps ending in a verdict
(yes / no / unclear) with an evidence level that can never exceed the quotes it rests on.
`--retrieve-only` shows the hit count without any model call.

## Reading the results

Quotes, PMIDs and chunk membership are verified by code. The labels `supports`,
`evidence_type`, `organism_is_yarrowia` and `compartment` are the model's judgment: use the
levels as triage and read the quotes before treating anything as evidence.
