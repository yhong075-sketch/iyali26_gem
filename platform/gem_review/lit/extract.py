"""Step 2: Pydantic AI extraction of validated claims from the retrieved chunks.

Per pair the agent sees the chunks as structured JSON in its instructions (never pasted into the
user message), answers with a `ClaimSet`, and an output validator rejects any claim whose PMID
is not in the retrieved set or whose quote is not verbatim in the cited chunk. After the retry
budget the pair is recorded as EXTRACTION_FAILED; no unvalidated answer is ever kept.
"""

import asyncio
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from pydantic_ai import Agent, ModelRetry, RunContext
from pydantic_ai.exceptions import UnexpectedModelBehavior
from pydantic_ai.models import Model

from .models import Claim, ClaimSet, PairReview, normalise_text
from .providers import DEFAULT_LLM, settings_for
from .terms import build_terms

logger = logging.getLogger(__name__)

DEFAULT_MODEL = DEFAULT_LLM
MAX_CONCURRENT = 5

INSTRUCTIONS = """\
You review the evidence for ONE gene-reaction pair in a genome-scale metabolic model of the yeast
Yarrowia lipolytica. The retrieved text chunks are given below as JSON.

Rules:
- Answer ONLY from the chunks. Never use background knowledge, even if you are sure of it.
- If the chunks do not establish (or contradict) that this gene product catalyses this reaction,
  return an empty claims list.
- Every claim cites one pmid and one chunk_id from the chunks, and a quote of at most 25 words
  copied verbatim from that chunk. Do not paraphrase.
- evidence_type: 4 biochemical (purified enzyme or in vitro activity with a defined substrate),
  3 genetic (knockout/overexpression phenotype in the organism), 2 physiological (growth or
  secretion phenotype consistent with the reaction), 1 homology (ortholog of a characterised
  enzyme), 0 none.
- organism_is_yarrowia is false when the quoted evidence was obtained in another organism.
- compartment only as stated in the paper, otherwise null.
"""


def chunks_json(chunks: dict[str, dict]) -> str:
    rows = [{"chunk_id": cid, "pmid": c["pmid"], "source": c["source"], "text": c["text"]}
            for cid, c in sorted(chunks.items())]
    return json.dumps(rows, ensure_ascii=False)


@dataclass(slots=True)
class PairContext:
    """Run dependencies: what the validator is allowed to accept."""

    gene: str
    reaction_id: str
    reaction_name: str
    ec_codes: list[str]
    pmids: set[str]
    chunks: dict[str, dict]  # chunk_id -> retrieval chunk record
    aliases: list[str] = field(default_factory=list)  # other spellings, YALI0 tag, gene symbols

    def context_json(self) -> str:
        return chunks_json(self.chunks)

    def question(self) -> str:
        """Papers name the gene by symbol or CLIB122 tag, never by the model's YALI1 ID."""
        ec = ", ".join(self.ec_codes) if self.ec_codes else "no EC"
        alias = f" (the same gene is also written {', '.join(self.aliases)})" if self.aliases else ""
        return (f"Gene {self.gene}{alias}; reaction {self.reaction_id} \"{self.reaction_name}\" ({ec}). "
                "Does the retrieved text establish that this gene product catalyses this reaction?")

    def share_key(self) -> tuple:
        """Only pairs asking the identical question share one model call.

        The full reaction name is part of the key on purpose: "acyl-CoA oxidase (C16)" and
        "(C4)" retrieve the same chunks, but the same POX3 quotes contradict the first and
        support the second, so chain-length variants must each be judged.
        """
        return (self.gene, self.reaction_name, tuple(self.ec_codes), tuple(sorted(self.chunks)))


def check_quote(pmids: set[str], chunks: dict[str, dict], pmid: str, chunk_id: str, quote: str) -> list[str]:
    """Problems with one cited quote: PMID retrieved, chunk retrieved and of that PMID, quote verbatim."""
    if pmid not in pmids:
        return [f"pmid {pmid} is not among the retrieved PMIDs"]
    chunk = chunks.get(chunk_id)
    if chunk is None:
        return [f"chunk_id {chunk_id} is not a retrieved chunk"]
    problems = []
    if chunk["pmid"] != pmid:
        problems.append(f"chunk {chunk_id} belongs to pmid {chunk['pmid']}, not {pmid}")
    if normalise_text(quote) not in normalise_text(chunk["text"]):
        problems.append(f"quote is not verbatim in chunk {chunk_id}")
    return problems


def validate_claims(ctx: PairContext, output: ClaimSet) -> ClaimSet:
    """Reject claims outside the retrieved set; raises ModelRetry with every problem listed."""
    problems = []
    for index, claim in enumerate(output.claims):
        problems += [f"claim {index}: {p}" for p in
                     check_quote(ctx.pmids, ctx.chunks, claim.pmid, claim.chunk_id, claim.quote)]
    if problems:
        raise ModelRetry("Fix these claims or drop them:\n" + "\n".join(problems))
    return output


def build_agent(model: str | Model = DEFAULT_MODEL, *, retries: int = 2) -> Agent[PairContext, ClaimSet]:
    settings = settings_for(model) if isinstance(model, str) else None
    agent: Agent[PairContext, ClaimSet] = Agent(
        model, output_type=ClaimSet, deps_type=PairContext, instructions=INSTRUCTIONS,
        retries=retries, model_settings=settings, name="gem_review_lit_extract")

    @agent.instructions
    def chunk_context(ctx: RunContext[PairContext]) -> str:
        return "Retrieved chunks (JSON):\n" + ctx.deps.context_json()

    @agent.output_validator
    def _validate(ctx: RunContext[PairContext], output: ClaimSet) -> ClaimSet:
        return validate_claims(ctx.deps, output)

    return agent


@dataclass(slots=True)
class PairResult:
    gene: str
    reaction_id: str
    status: str                      # ok | no_hits | EXTRACTION_FAILED | MODEL_ERROR
    review: PairReview | None = None
    error: str | None = None
    usage: dict = field(default_factory=dict)
    shared_with: str | None = None   # "<gene>__<reaction>" whose identical question was reused
    model: str | None = None

    def to_dict(self) -> dict:
        return {"gene": self.gene, "reaction_id": self.reaction_id, "status": self.status,
                "review": self.review.model_dump(mode="json") if self.review else None,
                "error": self.error, "usage": self.usage, "shared_with": self.shared_with,
                "model": self.model, "extracted_at": datetime.now(timezone.utc).isoformat()}


def pair_context(record: dict) -> PairContext:
    """From a step-1 per-pair retrieval JSON; aliases are the gene search terms minus the model ID."""
    terms = build_terms(record["gene"], record.get("yali0"), record.get("symbols", []), "", ())
    aliases = [t for t in terms.gene_terms if t != record["gene"]]
    return PairContext(record["gene"], record["reaction_id"], record["reaction_name"],
                       list(record["ec_codes"]), set(record["pmids"]),
                       {c["chunk_id"]: c for c in record["chunks"]}, aliases)


async def extract_pair(agent: Agent[PairContext, ClaimSet], ctx: PairContext, *,
                       semaphore: asyncio.Semaphore, model_name: str) -> PairResult:
    async with semaphore:
        try:
            run = await agent.run(ctx.question(), deps=ctx)
        except UnexpectedModelBehavior as exc:
            logger.warning("%s %s: EXTRACTION_FAILED (%s)", ctx.gene, ctx.reaction_id, exc)
            return PairResult(ctx.gene, ctx.reaction_id, "EXTRACTION_FAILED", error=str(exc), model=model_name)
        except Exception as exc:  # provider/network errors: recorded, never silently dropped
            logger.error("%s %s: MODEL_ERROR %s: %s", ctx.gene, ctx.reaction_id, type(exc).__name__, exc)
            return PairResult(ctx.gene, ctx.reaction_id, "MODEL_ERROR",
                              error=f"{type(exc).__name__}: {exc}", model=model_name)
    usage = {"input_tokens": run.usage.input_tokens, "output_tokens": run.usage.output_tokens,
             "cache_read_tokens": run.usage.cache_read_tokens, "requests": run.usage.requests}
    logger.info("%s %s: %d claims, tokens in=%d out=%d requests=%d", ctx.gene, ctx.reaction_id,
                len(run.output.claims), usage["input_tokens"], usage["output_tokens"], usage["requests"])
    review = PairReview.from_claims(ctx.gene, ctx.reaction_id, list(run.output.claims))
    return PairResult(ctx.gene, ctx.reaction_id, "ok", review=review, usage=usage, model=model_name)


async def run_extraction(records: list[dict], *, agent: Agent[PairContext, ClaimSet], model_name: str,
                         out_dir: Path, max_concurrent: int = MAX_CONCURRENT) -> list[PairResult]:
    """Pairs with no hits get a `no_hits` result without a model call. Identical questions share one call."""
    if max_concurrent > MAX_CONCURRENT:
        raise ValueError(f"max_concurrent is capped at {MAX_CONCURRENT}")
    out_dir.mkdir(parents=True, exist_ok=True)
    semaphore = asyncio.Semaphore(max_concurrent)
    results: list[PairResult] = []
    pending: dict[tuple, asyncio.Task] = {}
    followers: list[tuple[dict, tuple]] = []
    for record in records:
        if int(record.get("pair_hits", 0)) == 0 or not record["chunks"]:
            results.append(PairResult(record["gene"], record["reaction_id"], "no_hits",
                                      review=PairReview.from_claims(record["gene"], record["reaction_id"], [])))
            continue
        ctx = pair_context(record)
        key = ctx.share_key()
        if key in pending:
            followers.append((record, key))
            continue
        pending[key] = asyncio.create_task(extract_pair(agent, ctx, semaphore=semaphore, model_name=model_name))
    leaders = {key: await task for key, task in pending.items()}
    results.extend(leaders.values())
    for record, key in followers:
        lead = leaders[key]
        review = None if lead.review is None else PairReview.from_claims(
            record["gene"], record["reaction_id"], list(lead.review.claims))
        results.append(PairResult(record["gene"], record["reaction_id"], lead.status, review=review,
                                  error=lead.error, usage={}, model=lead.model,
                                  shared_with=f"{lead.gene}__{lead.reaction_id}"))
    with (out_dir / "usage.jsonl").open("a") as usage_log:
        for result in results:
            with (out_dir / f"{result.gene}__{result.reaction_id}.review.json").open("w") as handle:
                json.dump(result.to_dict(), handle, indent=1)
            if result.usage:
                usage_log.write(json.dumps({"gene": result.gene, "reaction_id": result.reaction_id,
                                            "model": result.model, **result.usage}) + "\n")
    return results


def load_retrieval(retrieval_dir: Path, *, only: set[tuple[str, str]] | None = None) -> list[dict]:
    records = []
    for path in sorted(retrieval_dir.glob("*.json")):
        with path.open() as handle:
            record = json.load(handle)
        if only is None or (record["gene"], record["reaction_id"]) in only:
            records.append(record)
    return records


def load_results(out_dir: Path) -> list[dict]:
    rows = []
    for path in sorted(out_dir.glob("*.review.json")):
        with path.open() as handle:
            rows.append(json.load(handle))
    return rows


__all__ = ["Claim", "ClaimSet", "PairContext", "PairResult", "PairReview", "build_agent",
           "check_quote", "extract_pair", "load_results", "load_retrieval", "run_extraction", "validate_claims"]
