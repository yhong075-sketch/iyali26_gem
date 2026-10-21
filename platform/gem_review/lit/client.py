"""Cached, rate-capped HTTP client for Europe PMC and UniProt (no LLM involved).

Every GET is keyed by SHA-256 of (url, params) and stored on disk, so a re-run makes no
network request for anything already fetched. Only 200 and 404 responses are cached: a 404
on fullTextXML means "no full text", which is a stable answer.
"""

import asyncio
import hashlib
import json
import logging
import time
from datetime import datetime, timezone
from pathlib import Path

import httpx

logger = logging.getLogger(__name__)

EPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EPMC_FULLTEXT_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

_RETRYABLE = {429, 500, 502, 503, 504}
_CACHEABLE = {200, 404}


class DiskCache:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def key(url: str, params: dict | None) -> str:
        payload = json.dumps([url, sorted((params or {}).items())], sort_keys=True)
        return hashlib.sha256(payload.encode()).hexdigest()

    def path(self, key: str) -> Path:
        return self.root / key[:2] / f"{key}.json"

    def get(self, key: str) -> dict | None:
        path = self.path(key)
        if not path.exists():
            return None
        with path.open() as handle:
            return json.load(handle)

    def put(self, key: str, record: dict) -> None:
        path = self.path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        with tmp.open("w") as handle:
            json.dump(record, handle)
        tmp.replace(path)


class CachedClient:
    """GET with disk cache, at most `max_concurrent` in flight and one start per `min_interval`."""

    def __init__(self, cache: DiskCache, *, transport=None, max_concurrent: int = 5,
                 min_interval: float = 0.34, retries: int = 3, timeout: float = 60.0):
        if max_concurrent > 5:
            raise ValueError("max_concurrent is capped at 5")
        self.cache = cache
        self._client = httpx.AsyncClient(transport=transport, timeout=timeout,
                                         headers={"User-Agent": "iyali26-gem-review/0.1"})
        self._semaphore = asyncio.Semaphore(max_concurrent)
        self._lock = asyncio.Lock()
        self._next_start = 0.0
        self._min_interval = min_interval
        self._retries = retries
        self.network_requests = 0
        self.cache_hits = 0

    async def aclose(self) -> None:
        await self._client.aclose()

    async def _throttle(self) -> None:
        async with self._lock:
            now = time.monotonic()
            wait = self._next_start - now
            if wait > 0:
                await asyncio.sleep(wait)
            self._next_start = max(now, self._next_start) + self._min_interval

    async def get(self, url: str, params: dict | None = None) -> tuple[int, str]:
        """Return (status, body). Cached 200/404 answers never touch the network."""
        key = self.cache.key(url, params)
        cached = self.cache.get(key)
        if cached is not None:
            self.cache_hits += 1
            return cached["status"], cached["body"]
        async with self._semaphore:
            response = None
            for attempt in range(self._retries):
                await self._throttle()
                response = await self._client.get(url, params=params)
                self.network_requests += 1
                if response.status_code not in _RETRYABLE:
                    break
                delay = 2 ** attempt
                if response.status_code == 429:
                    try:
                        delay = int(response.headers.get("Retry-After", delay))
                    except (TypeError, ValueError):
                        pass
                logger.debug("HTTP %s on %s; retry %d/%d in %ss", response.status_code, url,
                             attempt + 1, self._retries, delay)
                await asyncio.sleep(delay)
        status, body = response.status_code, response.text
        if status in _CACHEABLE:
            self.cache.put(key, {"url": url, "params": params or {}, "status": status,
                                 "fetched_at": datetime.now(timezone.utc).isoformat(), "body": body})
        return status, body


async def epmc_search(client: CachedClient, query: str, *, result_type: str = "core",
                      page_size: int = 25) -> dict:
    params = {"query": query, "format": "json", "resultType": result_type, "pageSize": page_size}
    status, body = await client.get(EPMC_SEARCH_URL, params)
    if status != 200:
        raise RuntimeError(f"Europe PMC search failed ({status}) for query: {query}")
    return json.loads(body)


async def epmc_hit_count(client: CachedClient, query: str) -> int:
    data = await epmc_search(client, query, result_type="lite", page_size=1)
    return int(data["hitCount"])


async def epmc_full_text(client: CachedClient, pmcid: str) -> str | None:
    """JATS XML for an open-access article, or None when Europe PMC has none."""
    status, body = await client.get(EPMC_FULLTEXT_URL.format(pmcid=pmcid))
    if status == 404:
        return None
    if status != 200:
        raise RuntimeError(f"Europe PMC fullTextXML failed ({status}) for {pmcid}")
    return body


async def uniprot_gene_symbols(client: CachedClient, accession: str) -> list[str]:
    """Primary gene name and synonyms (not ORF names, which are locus tags)."""
    status, body = await client.get(UNIPROT_ENTRY_URL.format(accession=accession),
                                    {"fields": "gene_names"})
    if status == 404:
        return []
    if status != 200:
        raise RuntimeError(f"UniProt lookup failed ({status}) for {accession}")
    symbols = []
    for gene in json.loads(body).get("genes", []):
        if "geneName" in gene:
            symbols.append(gene["geneName"]["value"])
        symbols.extend(s["value"] for s in gene.get("synonyms", []))
    return list(dict.fromkeys(symbols))
