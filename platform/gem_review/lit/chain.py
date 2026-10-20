"""Logic-chain layer: one question answered as quote → inference → verdict steps.

The agent sees the retrieved chunks as JSON in its instructions (as in extract.py) and returns a
`LogicChain`. Structure is checked by Pydantic (every inference traces to a quote, the verdict's
level never exceeds the quotes it rests on); every quote is checked against the chunks by the
same rules as `validate_claims`. After the retry budget the run is recorded as CHAIN_FAILED.
"""

import hashlib
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, Field, field_validator, model_validator
from pydantic_ai import Agent, ModelRetry, RunContext
from pydantic_ai.exceptions import UnexpectedModelBehavior
from pydantic_ai.models import Model
from pydantic_ai.models.anthropic import AnthropicModelSettings

from .client import CachedClient
from .extract import DEFAULT_MODEL, check_quote, chunks_json, pair_context
from .models import MAX_QUOTE_WORDS, EvidenceType, short_quote
from .pairs import Pair
from .retrieve import retrieve_pair
from .terms import PairTerms, locus_spellings

logger = logging.getLogger(__name__)

INSTRUCTIONS = f"""\
You answer ONE question for the review of a genome-scale metabolic model of the yeast Yarrowia
lipolytica, as a logic chain built only from the retrieved text chunks given below as JSON.

The chain is an ordered list of steps, numbered from 0:
- quote: one pmid and one chunk_id from the chunks and at most {MAX_QUOTE_WORDS} words copied
  verbatim from that chunk. evidence_type: 4 biochemical (purified enzyme or in vitro activity with
  a defined substrate), 3 genetic (knockout/overexpression phenotype in the organism),
  2 physiological (growth or secretion phenotype), 1 homology (ortholog of a characterised
  enzyme), 0 none. organism_is_yarrowia is false when the quoted evidence was obtained in another
  organism.
- inference: one short reasoning step; rests_on lists the indices of the earlier steps it uses.
  Every inference must rest, directly or through earlier inferences, on at least one quote.
- verdict: the last step and the only verdict. answer is yes, no or unclear; evidence_level is the
  strength of the evidence for the verdict as stated and may not exceed the strongest quote it
  rests on (a quote from another organism counts at most as 1); rests_on lists the steps it uses;
  verdict is one sentence.

Rules:
- Answer ONLY from the chunks. Never use background knowledge, even if you are sure of it.
- Quote text that contradicts your answer as well, when the chunks contain it.
- A sentence that names no gene is not evidence about a specific gene; say so in an inference
  instead of treating it as gene-specific.
- If the chunks do not bear on the question, return only a verdict with answer unclear,
  evidence_level 0 and an empty rests_on.
"""


class QuoteStep(BaseModel):
    kind: Literal["quote"]
    pmid: str = Field(description="PMID of the cited article; must be one of the retrieved PMIDs")
    chunk_id: str = Field(description="chunk_id of the chunk the quote is copied from")
    quote: str = Field(description=f"At most {MAX_QUOTE_WORDS} words copied verbatim from the cited chunk")
    evidence_type: EvidenceType
    organism_is_yarrowia: bool = Field(description="False if the evidence comes from another organism")

    @field_validator("quote")
    @classmethod
    def _short_quote(cls, value: str) -> str:
        return short_quote(value)

    @property
    def effective_level(self) -> EvidenceType:
        """Evidence from another organism counts at most as homology (as for Claim)."""
        if not self.organism_is_yarrowia:
            return min(self.evidence_type, EvidenceType.homology)
        return self.evidence_type


class InferenceStep(BaseModel):
    kind: Literal["inference"]
    text: str = Field(min_length=1)
    rests_on: list[int] = Field(min_length=1, description="indices of the earlier steps this inference uses")


class VerdictStep(BaseModel):
    kind: Literal["verdict"]
    answer: Literal["yes", "no", "unclear"]
    evidence_level: EvidenceType = Field(description="strength of the evidence for the verdict as stated")
    verdict: str = Field(min_length=1, description="one sentence")
    rests_on: list[int] = Field(default_factory=list, description="indices of the steps the verdict uses")


Step = Annotated[QuoteStep | InferenceStep | VerdictStep, Field(discriminator="kind")]


class LogicChain(BaseModel):
    steps: list[Step] = Field(min_length=1, description="quote and inference steps, then one verdict")

    @model_validator(mode="after")
    def _check_structure(self) -> "LogicChain":
        problems = chain_problems(self.steps)
        if problems:
            raise ValueError("; ".join(problems))
        return self

    @property
    def verdict(self) -> VerdictStep:
        return self.steps[-1]

    @property
    def quotes(self) -> list[QuoteStep]:
        return [s for s in self.steps if isinstance(s, QuoteStep)]


def _supporting_quotes(steps: list, roots: list[int]) -> list[QuoteStep]:
    """Quotes reachable from `roots` through rests_on (indices already range-checked)."""
    seen, stack, found = set(), list(roots), []
    while stack:
        index = stack.pop()
        if index in seen:
            continue
        seen.add(index)
        step = steps[index]
        if isinstance(step, QuoteStep):
            found.append(step)
        elif isinstance(step, InferenceStep):
            stack.extend(step.rests_on)
    return found


def chain_problems(steps: list) -> list[str]:
    """Structural problems. rests_on may only point backwards, so tracing step by step is transitive."""
    problems = []
    traced: set[int] = set()  # quotes, and inferences resting on a traced step
    for index, step in enumerate(steps):
        last = index == len(steps) - 1
        if isinstance(step, VerdictStep) and not last:
            problems.append(f"step {index}: the verdict must be the last step and appear once")
        if last and not isinstance(step, VerdictStep):
            problems.append("the chain must end with a verdict step")
        if isinstance(step, QuoteStep):
            traced.add(index)
            continue
        bad = [j for j in step.rests_on if not 0 <= j < index]
        if bad:
            problems.append(f"step {index}: rests_on {bad} must point to earlier steps")
            continue
        if isinstance(step, InferenceStep):
            if any(j in traced for j in step.rests_on):
                traced.add(index)
            else:
                problems.append(f"step {index}: inference does not trace to any quote")
        elif last:
            support = _supporting_quotes(steps, step.rests_on)
            cap = max((q.effective_level for q in support), default=EvidenceType.none)
            if step.evidence_level > cap:
                problems.append(
                    f"verdict evidence_level {int(step.evidence_level)} exceeds {int(cap)}, the strongest "
                    "quote it rests on (quotes from another organism count at most as homology = 1; "
                    "a verdict resting on no quote must be level 0)")
    return problems


@dataclass(slots=True)
class ChainContext:
    """Run dependencies: the prompt and what the quote validator is allowed to accept."""

    prompt: str
    pmids: set[str]
    chunks: dict[str, dict]


def validate_chain(ctx: ChainContext, output: LogicChain) -> LogicChain:
    problems = []
    for index, step in enumerate(output.steps):
        if isinstance(step, QuoteStep):
            problems += [f"step {index}: {p}" for p in
                         check_quote(ctx.pmids, ctx.chunks, step.pmid, step.chunk_id, step.quote)]
    if problems:
        raise ModelRetry("Fix these quote steps or drop them (and the steps resting on them):\n"
                         + "\n".join(problems))
    return output


def build_chain_agent(model: str | Model = DEFAULT_MODEL, *, retries: int = 2) -> Agent[ChainContext, LogicChain]:
    settings = AnthropicModelSettings(max_tokens=16000, anthropic_effort="medium") \
        if isinstance(model, str) and model.startswith("anthropic:") else None
    agent: Agent[ChainContext, LogicChain] = Agent(
        model, output_type=LogicChain, deps_type=ChainContext, instructions=INSTRUCTIONS,
        retries=retries, model_settings=settings, name="gem_review_lit_chain")

    @agent.instructions
    def chunk_context(ctx: RunContext[ChainContext]) -> str:
        return "Retrieved chunks (JSON):\n" + chunks_json(ctx.deps.chunks)

    @agent.output_validator
    def _validate(ctx: RunContext[ChainContext], output: LogicChain) -> LogicChain:
        return validate_chain(ctx.deps, output)

    return agent


@dataclass(slots=True)
class ChainResult:
    mode: str                        # pair | question
    label: str                       # "<gene>__<reaction>" or the question slug
    question: str
    prompt: str
    status: str                      # pending | ok | no_hits | CHAIN_FAILED | MODEL_ERROR
    gene: str | None = None
    reaction_id: str | None = None
    retrieval_file: str | None = None
    chain: LogicChain | None = None
    error: str | None = None
    usage: dict = field(default_factory=dict)
    model: str | None = None

    def to_dict(self) -> dict:
        return {"mode": self.mode, "label": self.label, "question": self.question, "prompt": self.prompt,
                "gene": self.gene, "reaction_id": self.reaction_id, "status": self.status,
                "retrieval_file": self.retrieval_file,
                "chain": self.chain.model_dump(mode="json") if self.chain else None,
                "error": self.error, "usage": self.usage, "model": self.model,
                "chained_at": datetime.now(timezone.utc).isoformat()}


async def run_chain(agent: Agent[ChainContext, LogicChain], ctx: ChainContext, result: ChainResult, *,
                    model_name: str) -> ChainResult:
    """Fill `result` from one agent run; failures are recorded, never dropped."""
    result.model = model_name
    try:
        run = await agent.run(ctx.prompt, deps=ctx)
    except UnexpectedModelBehavior as exc:
        logger.warning("%s: CHAIN_FAILED (%s)", result.label, exc)
        result.status, result.error = "CHAIN_FAILED", str(exc)
        return result
    except Exception as exc:  # provider/network errors: recorded, never silently dropped
        logger.error("%s: MODEL_ERROR %s: %s", result.label, type(exc).__name__, exc)
        result.status, result.error = "MODEL_ERROR", f"{type(exc).__name__}: {exc}"
        return result
    result.status, result.chain = "ok", run.output
    result.usage = {"input_tokens": run.usage.input_tokens, "output_tokens": run.usage.output_tokens,
                    "cache_read_tokens": run.usage.cache_read_tokens, "requests": run.usage.requests}
    return result


def pair_chain_inputs(record: dict, retrieval_file: Path | None = None) -> tuple[ChainContext, ChainResult]:
    """From a step-1 per-pair retrieval JSON; the question lists the aliases papers use."""
    pctx = pair_context(record)
    question = pctx.question()
    ctx = ChainContext(question, pctx.pmids, pctx.chunks)
    result = ChainResult("pair", f"{pctx.gene}__{pctx.reaction_id}", question, question, "pending",
                         gene=pctx.gene, reaction_id=pctx.reaction_id,
                         retrieval_file=str(retrieval_file) if retrieval_file else None)
    return ctx, result


def question_slug(question: str, gene: str | None, aliases: list[str], terms: list[str]) -> str:
    payload = json.dumps([question, gene, aliases, terms])
    return "q_" + hashlib.sha256(payload.encode()).hexdigest()[:12]


def question_terms(gene: str | None, aliases: list[str], terms: list[str]) -> PairTerms:
    gene_terms = []
    for name in ([gene] if gene else []) + aliases:
        gene_terms += locus_spellings(name)
    return PairTerms(tuple(dict.fromkeys(gene_terms)), tuple(dict.fromkeys(terms)))


async def retrieve_question(client: CachedClient, question: str, *, gene: str | None, aliases: list[str],
                            terms: list[str], page_size: int = 25, max_fulltext: int = 10) -> dict:
    """Step-1 retrieval for an ad-hoc question; returns a record shaped like the per-pair JSON."""
    slug = question_slug(question, gene, aliases, terms)
    pair = Pair(gene or aliases[0], slug, question, "", (), None)
    found = await retrieve_pair(client, pair, question_terms(gene, aliases, terms), yali0=None,
                                symbols=aliases, page_size=page_size, max_fulltext=max_fulltext)
    return {**found.to_dict(), "gene": gene, "question": question, "aliases": aliases, "terms": terms,
            "slug": slug}


def question_chain_inputs(record: dict, retrieval_file: Path) -> tuple[ChainContext, ChainResult]:
    names = [n for n in [record["gene"], *record["aliases"]] if n]
    prompt = f"{record['question']} (names this gene goes by in papers: {', '.join(dict.fromkeys(names))})"
    ctx = ChainContext(prompt, set(record["pmids"]), {c["chunk_id"]: c for c in record["chunks"]})
    result = ChainResult("question", record["slug"], record["question"], prompt, "pending",
                         gene=record["gene"], retrieval_file=str(retrieval_file))
    return ctx, result


def hit_line(record: dict) -> str:
    sources = [c["source"] for c in record["chunks"]]
    return (f"gene_hits={record['gene_hits']} pair_hits={record['pair_hits']} PMIDs={len(record['pmids'])} "
            f"chunks={len(sources)} (abstract {sources.count('abstract')}, fulltext {sources.count('fulltext')})\n"
            f"pair query: {record['pair_query']}")


def chain_path(out_dir: Path, result: ChainResult) -> Path:
    if result.mode == "question":
        return out_dir / "questions" / f"{result.label}.chain.json"
    return out_dir / f"{result.label}.chain.json"


def write_chain(out_dir: Path, result: ChainResult) -> Path:
    path = chain_path(out_dir, result)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(result.to_dict(), handle, indent=1)
    return path


def load_chains(out_dir: Path) -> list[dict]:
    rows = []
    for path in sorted([*out_dir.glob("*.chain.json"), *(out_dir / "questions").glob("*.chain.json")]):
        with path.open() as handle:
            rows.append(json.load(handle))
    return rows


__all__ = ["ChainContext", "ChainResult", "InferenceStep", "LogicChain", "QuoteStep", "VerdictStep",
           "build_chain_agent", "chain_problems", "hit_line", "load_chains", "pair_chain_inputs",
           "question_chain_inputs", "retrieve_question", "run_chain", "validate_chain", "write_chain"]
