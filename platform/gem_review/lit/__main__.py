"""CLI for the literature layer.

  python -m gem_review.lit retrieve --model <xml> (--tier 1 | --reactions R1,R2) [--task T]
  python -m gem_review.lit extract  --task T [--llm anthropic:claude-opus-5-5] [--pairs g:r,...]
  python -m gem_review.lit report   --task T
"""

import argparse
import asyncio
import logging
import os
import sys
from collections import Counter
from datetime import date
from pathlib import Path

from scripts.gem_annotate.config import load_project_paths
from scripts.gem_annotate.locus_resolver import _DEFAULT_IDENTITY_EXCLUSIONS, LocusCrosswalk
from scripts.gem_annotate.model_layout import MODEL

from .client import CachedClient, DiskCache
from .extract import DEFAULT_MODEL, build_agent, load_results, load_retrieval, run_extraction
from .pairs import geneless, list_pairs, select_reactions
from .report import literature_rows, review_queue_rows, write_review_xlsx
from .retrieve import retrieve_all, sha256_of, write_outputs


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--research-root", help="overrides IYALI26_RESEARCH_ROOT")
    parser.add_argument("--task", default=f"gem_review_lit_{date.today():%Y%m%d}",
                        help="task folder name under model/reports and artifacts/tasks")


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="python -m gem_review.lit")
    sub = parser.add_subparsers(dest="command", required=True)
    retrieve = sub.add_parser("retrieve", help="Europe PMC retrieval and hit-count table (no LLM)")
    retrieve.add_argument("--model", required=True, type=Path, help="SBML model to read pairs from")
    retrieve.add_argument("--tier", type=int, help="priority tier (only 1 is defined)")
    retrieve.add_argument("--reactions", help="comma-separated reaction IDs")
    _common(retrieve)
    retrieve.add_argument("--page-size", type=int, default=25)
    retrieve.add_argument("--max-fulltext", type=int, default=10)
    retrieve.add_argument("--max-concurrent", type=int, default=5)

    extract = sub.add_parser("extract", help="Pydantic AI claim extraction on pairs with hits")
    _common(extract)
    extract.add_argument("--llm", default=DEFAULT_MODEL,
                         help="pydantic-ai model string, e.g. anthropic:claude-opus-5-5 (key from env)")
    extract.add_argument("--pairs", help="comma-separated gene:reaction pairs to restrict to")
    extract.add_argument("--limit", type=int, help="stop after this many records (smoke runs)")
    extract.add_argument("--max-concurrent", type=int, default=5)
    extract.add_argument("--retries", type=int, default=2)

    report = sub.add_parser("report", help="write the literature and review_queue sheets to review.xlsx")
    _common(report)

    args = parser.parse_args(argv)
    if args.command == "retrieve" and not args.reactions and args.tier is None:
        parser.error("pass --tier and/or --reactions; the tool never runs on the whole model")
    return args


def _task_dirs(args) -> tuple[Path, Path]:
    paths = load_project_paths(args.research_root, required=True)
    task_dir = paths.task_outputs / args.task
    return task_dir / "out" / "lit" / "retrieval", task_dir / "out" / "lit"


def extract_main(args) -> int:
    if args.llm.startswith("anthropic:") and not os.environ.get("ANTHROPIC_API_KEY"):
        print("ANTHROPIC_API_KEY is not set; export it before running extraction.", file=sys.stderr)
        return 2
    retrieval_dir, out_dir = _task_dirs(args)
    only = {tuple(p.split(":", 1)) for p in args.pairs.split(",")} if args.pairs else None
    records = load_retrieval(retrieval_dir, only=only)
    if args.limit:
        records = records[:args.limit]
    with_hits = sum(int(r["pair_hits"]) > 0 for r in records)
    logging.info("%d records, %d with hits -> model calls; LLM %s", len(records), with_hits, args.llm)
    agent = build_agent(args.llm, retries=args.retries)
    results = asyncio.run(run_extraction(records, agent=agent, model_name=args.llm, out_dir=out_dir,
                                         max_concurrent=args.max_concurrent))
    status = Counter(r.status for r in results)
    levels = Counter(r.review.best_evidence_level.name for r in results if r.review)
    tokens_in = sum(r.usage.get("input_tokens", 0) for r in results)
    tokens_out = sum(r.usage.get("output_tokens", 0) for r in results)
    print(f"status: {dict(status)}")
    print(f"best_evidence_level: {dict(levels)}")
    print(f"tokens: input={tokens_in} output={tokens_out}; per-pair log: {out_dir / 'usage.jsonl'}")
    return 0


def report_main(args) -> int:
    retrieval_dir, out_dir = _task_dirs(args)
    rows = literature_rows(load_retrieval(retrieval_dir), load_results(out_dir))
    queue = review_queue_rows(rows)
    path = write_review_xlsx(MODEL.reports / args.task / "review.xlsx", rows, queue)
    print(f"{len(rows)} literature rows, {len(queue)} review_queue rows -> {path}")
    return 0


def _summary(rows: list[dict]) -> str:
    bucket = Counter()
    for row in rows:
        hits = row["pair_hits"]
        bucket["0" if hits == 0 else "1-4" if hits < 5 else "5-19" if hits < 20 else "20+"] += 1
    with_ft = sum(row["fulltext_available"] == "Y" for row in rows)
    lines = [f"pairs: {len(rows)}",
             f"pair_hits  0: {bucket['0']}   1-4: {bucket['1-4']}   5-19: {bucket['5-19']}   20+: {bucket['20+']}",
             f"pairs with at least one open-access full text fetched: {with_ft}"]
    return "\n".join(lines)


def main(argv=None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    logging.getLogger("httpx").setLevel(logging.WARNING)
    os.environ.setdefault("PYDANTIC_AI_NO_BANNER", "1")
    if args.command == "extract":
        return extract_main(args)
    if args.command == "report":
        return report_main(args)
    import cobra

    paths = load_project_paths(args.research_root, required=True)
    model = cobra.io.read_sbml_model(str(args.model))
    reaction_ids = [r.strip() for r in args.reactions.split(",")] if args.reactions else None
    reactions = select_reactions(model, reaction_ids=reaction_ids, tier=args.tier)
    pairs = list_pairs(reactions, args.tier)
    no_gpr = geneless(reactions)
    logging.info("%d reactions, %d gene-reaction pairs, %d reactions without GPR",
                 len(reactions), len(pairs), len(no_gpr))

    crosswalk = LocusCrosswalk.from_csvs(paths.locus_map / "iyli21_genes_vs_S2.csv",
                                         paths.locus_map / "s2_metabolic_genes.csv",
                                         identity_exclusions_path=_DEFAULT_IDENTITY_EXCLUSIONS)

    def accession(gene_id: str) -> str | None:
        value = model.genes.get_by_id(gene_id).annotation.get("uniprot")
        return value[0] if isinstance(value, list) else value

    task_dir = paths.task_outputs / args.task
    client = CachedClient(DiskCache(task_dir / "cache"), max_concurrent=args.max_concurrent)

    async def run():
        try:
            return await retrieve_all(client, pairs, counterpart=crosswalk.counterpart,
                                      uniprot_accession=accession, page_size=args.page_size,
                                      max_fulltext=args.max_fulltext)
        finally:
            await client.aclose()

    results = asyncio.run(run())
    manifest = {
        "model": str(args.model), "model_sha256": sha256_of(args.model),
        "tier": args.tier, "reactions_arg": reaction_ids, "n_reactions": len(reactions),
        "reactions_without_gpr": no_gpr, "crosswalk_fingerprint": crosswalk.fingerprint,
        "page_size": args.page_size, "max_fulltext": args.max_fulltext,
        "network_requests": client.network_requests, "cache_hits": client.cache_hits,
    }
    table = write_outputs(results, out_dir=task_dir / "out" / "lit" / "retrieval",
                          report_dir=MODEL.reports / args.task, manifest=manifest)
    print(_summary([r.hit_row() for r in results]))
    print(f"hit table: {table}")
    print(f"network requests: {client.network_requests}, cache hits: {client.cache_hits}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
