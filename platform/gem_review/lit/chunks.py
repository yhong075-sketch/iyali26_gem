"""Text chunks handed to the extraction step: abstracts and matching full-text paragraphs."""

import re
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass

_WS = re.compile(r"\s+")


@dataclass(frozen=True, slots=True)
class Chunk:
    chunk_id: str       # "<pmid>:abstract" or "<pmcid>:p<n>"
    pmid: str
    pmcid: str | None
    source: str         # "abstract" | "fulltext"
    text: str
    matched_terms: tuple[str, ...]

    def to_dict(self) -> dict:
        return asdict(self)


def _normalise(text: str) -> str:
    return _WS.sub(" ", text).strip()


def matched(text: str, terms) -> tuple[str, ...]:
    lowered = text.lower()
    return tuple(term for term in terms if term.lower() in lowered)


def abstract_chunk(result: dict) -> Chunk | None:
    text = _normalise(result.get("abstractText") or "")
    pmid = result.get("pmid")
    if not text or not pmid:
        return None
    return Chunk(f"{pmid}:abstract", pmid, result.get("pmcid"), "abstract", text, ())


def fulltext_chunks(xml_text: str, *, pmid: str, pmcid: str, terms) -> list[Chunk]:
    """Body paragraphs (in document order) that mention at least one gene or reaction term."""
    root = ET.fromstring(xml_text)
    body = root.find(".//body")
    if body is None:
        return []
    chunks = []
    for index, paragraph in enumerate(body.iter("p")):
        text = _normalise("".join(paragraph.itertext()))
        hits = matched(text, terms)
        if hits:
            chunks.append(Chunk(f"{pmcid}:p{index}", pmid, pmcid, "fulltext", text, hits))
    return chunks
