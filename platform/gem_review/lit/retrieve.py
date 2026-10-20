"""Step 1: retrieve Europe PMC evidence per gene–reaction pair and write the hit-count table."""

import asyncio
import csv
import hashlib
import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from .chunks import Chunk, abstract_chunk, fulltext_chunks
from .client import CachedClient, epmc_full_text, epmc_hit_count, epmc_search, uniprot_gene_symbols
from .pairs import Pair
from .terms import PairTerms, build_terms, locus_spellings

logger = logging.getLogger(__name__)

HIT_TABLE_COLUMNS = [
    "reaction_id", "reaction_name", "ec", "subsystem", "gene", "yali0", "symbols",
    "gene_hits", "pair_hits", "pair_oa_hits", "fulltext_fetched", "fulltext_available",
    "n_abstract_chunks", "n_fulltext_chunks", "pmids",
]


@dataclass(slots=True)
class PairRetrieval:
    gene: str
    reaction_id: str
    reaction_name: str
    ec_codes: list[str]
    subsystem: str
    tier: int | None
    yali0: str | None
    symbols: list[str]
    gene_query: str
    pair_query: str | None
    gene_hits: int = 0
    pair_hits: int = 0
    pair_oa_hits: int = 0
    pmids: list[str] = field(default_factory=list)
    pmcids: list[str] = field(default_factory=list)
    fulltext_fetched: list[str] = field(default_factory=list)
    chunks: list[Chunk] = field(default_factory=list)

    @property
    def fulltext_available(self) -> bool:
        return bool(self.fulltext_fetched)

    def to_dict(self) -> dict:
        data = asdict(self)
        data["fulltext_available"] = self.fulltext_available
        return data

    def hit_row(self) -> dict:
        return {
            "reaction_id": self.reaction_id, "reaction_name": self.reaction_name,
            "ec": ";".join(self.ec_codes), "subsystem": self.subsystem, "gene": self.gene,
            "yali0": self.yali0 or "", "symbols": ";".join(self.symbols),
            "gene_hits": self.gene_hits, "pair_hits": self.pair_hits,
            "pair_oa_hits": self.pair_oa_hits, "fulltext_fetched": len(self.fulltext_fetched),
            "fulltext_available": "Y" if self.fulltext_available else "N",
            "n_abstract_chunks": sum(c.source == "abstract" for c in self.chunks),
            "n_fulltext_chunks": sum(c.source == "fulltext" for c in self.chunks),
            "pmids": ";".join(self.pmids),
        }


async def retrieve_pair(client: CachedClient, pair: Pair, terms: PairTerms, *, yali0, symbols,
                        page_size: int = 25, max_fulltext: int = 10) -> PairRetrieval:
    out = PairRetrieval(pair.gene, pair.reaction_id, pair.reaction_name, list(pair.ec_codes),
                        pair.subsystem, pair.tier, yali0, list(symbols),
                        terms.gene_query(), terms.pair_query())
    out.gene_hits = await epmc_hit_count(client, out.gene_query)
    # Without a reaction term the pair query would equal the gene query: record it as no pair hits
    # and take no abstracts, so the table shows where the reaction is unsearchable.
    if out.pair_query is None:
        return out
    search = await epmc_search(client, out.pair_query, result_type="core", page_size=page_size)
    out.pair_hits = int(search["hitCount"])
    results = search.get("resultList", {}).get("result", [])
    all_terms = terms.gene_terms + terms.reaction_terms
    for result in results:
        pmid = result.get("pmid")
        if pmid:
            out.pmids.append(pmid)
        if result.get("pmcid"):
            out.pmcids.append(result["pmcid"])
        chunk = abstract_chunk(result)
        if chunk:
            out.chunks.append(Chunk(chunk.chunk_id, chunk.pmid, chunk.pmcid, chunk.source, chunk.text,
                                    tuple(t for t in all_terms if t.lower() in chunk.text.lower())))
    open_access = [r for r in results if r.get("isOpenAccess") == "Y" and r.get("pmcid") and r.get("pmid")]
    out.pair_oa_hits = len(open_access)
    for result in open_access[:max_fulltext]:
        xml_text = await epmc_full_text(client, result["pmcid"])
        if xml_text is None:
            continue
        out.fulltext_fetched.append(result["pmcid"])
        out.chunks.extend(fulltext_chunks(xml_text, pmid=result["pmid"], pmcid=result["pmcid"],
                                          terms=all_terms))
    return out


async def retrieve_all(client: CachedClient, pairs: list[Pair], *, counterpart, uniprot_accession,
                       page_size: int = 25, max_fulltext: int = 10) -> list[PairRetrieval]:
    """`counterpart(gene) -> yali0|None`; `uniprot_accession(gene) -> accession|None`."""
    symbols_by_gene: dict[str, list[str]] = {}

    async def symbols_for(gene: str) -> list[str]:
        if gene not in symbols_by_gene:
            accession = uniprot_accession(gene)
            symbols_by_gene[gene] = await uniprot_gene_symbols(client, accession) if accession else []
        return symbols_by_gene[gene]

    async def one(pair: Pair) -> PairRetrieval:
        yali0 = counterpart(pair.gene)
        if yali0:
            yali0 = locus_spellings(yali0)[0]
        symbols = await symbols_for(pair.gene)
        terms = build_terms(pair.gene, yali0, symbols, pair.reaction_name, pair.ec_codes)
        return await retrieve_pair(client, pair, terms, yali0=yali0, symbols=symbols,
                                   page_size=page_size, max_fulltext=max_fulltext)

    return list(await asyncio.gather(*(one(pair) for pair in pairs)))


def write_outputs(results: list[PairRetrieval], *, out_dir: Path, report_dir: Path,
                  manifest: dict) -> Path:
    """Per-pair JSON under out_dir, hit_counts.csv plus manifest.json under report_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    for result in results:
        with (out_dir / f"{result.gene}__{result.reaction_id}.json").open("w") as handle:
            json.dump(result.to_dict(), handle, indent=1)
    table = report_dir / "hit_counts.csv"
    rows = sorted((r.hit_row() for r in results), key=lambda row: (row["reaction_id"], row["gene"]))
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_TABLE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    with (report_dir / "manifest.json").open("w") as handle:
        json.dump({**manifest, "written_at": datetime.now(timezone.utc).isoformat(),
                   "n_pairs": len(results), "out_dir": str(out_dir)}, handle, indent=1)
    return table


def sha256_of(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
