"""Gene identifier map: YALI0 (CLIB122) ↔ YALI1 (W29) ↔ UniProt ↔ NCBI (spec §1.1; D6).

Built on the repo's S2-based `LocusCrosswalk` sources. Every cross-assembly link carries a
cardinality (1:1, 1:n, n:1, n:m). Only 1:1 links resolve; anything else is UNRESOLVED and an
identifier with no link at all is UNMAPPED. A YALI1 locus needs no link and resolves to itself.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Literal

from scripts.gem_annotate.locus_resolver import _parsed_yali_locus, canonical_locus_key

Status = Literal["resolved", "unresolved", "unmapped"]


def display_locus(raw: str) -> str | None:
    """SBML-style spelling (YALI1A00309g) or None for a non-YALI identifier."""
    parsed = _parsed_yali_locus(raw)
    return None if parsed is None else parsed[0] + parsed[1]


@dataclass(frozen=True)
class IdMapRow:
    model_gene: str                 # as written in the model
    canonical: str | None           # YALI1 ID used for comparison
    status: Status
    cardinality: str                # 1:1, 1:n, n:1, n:m, none (no cross-assembly link), n/a
    yali0: tuple[str, ...]
    uniprot: tuple[str, ...]
    ncbigene: tuple[str, ...]
    note: str


@dataclass(frozen=True)
class Normalized:
    canonical: str | None
    status: Status
    cardinality: str
    candidates: tuple[str, ...] = ()


@dataclass
class IdMap:
    """Cross-assembly graph (casefolded keys) plus model gene annotations."""

    graph: dict[str, set[str]]
    rows: dict[str, IdMapRow] = field(default_factory=dict)
    fingerprint: str = ""
    aliases: dict[str, set[str]] = field(default_factory=dict)   # UniProt / NCBI ID -> canonical genes

    @classmethod
    def from_pairs(cls, pairs: Iterable[tuple[str, str]], excluded: Iterable[tuple[str, str]] = (),
                   fingerprint: str = "") -> "IdMap":
        excluded = {(canonical_locus_key(a), canonical_locus_key(b)) for a, b in excluded}
        graph: dict[str, set[str]] = defaultdict(set)
        for yali1, yali0 in pairs:
            a, b = canonical_locus_key(yali1), canonical_locus_key(yali0)
            if (a, b) in excluded or a == b:
                continue
            graph[a].add(b)
            graph[b].add(a)
        return cls(graph=dict(graph), fingerprint=fingerprint)

    def cardinality(self, key: str) -> str:
        """Cardinality of the link seen from `key` (its own side first)."""
        partners = self.graph.get(key, set())
        if not partners:
            return "none"
        back = set().union(*(self.graph.get(p, set()) for p in partners))
        if len(partners) == 1 and back == {key}:
            return "1:1"
        if len(partners) == 1:
            return "n:1"
        if back == {key}:
            return "1:n"
        return "n:m"

    def normalize(self, raw: str) -> Normalized:
        display = display_locus(raw)
        if display is None:
            targets = sorted(self.aliases.get(raw.strip(), set()))
            if len(targets) == 1:
                return Normalized(targets[0], "resolved", "1:1")
            if targets:
                return Normalized(None, "unresolved", "1:n", tuple(targets))
            return Normalized(None, "unmapped", "n/a")
        key = canonical_locus_key(display)
        cardinality = self.cardinality(key)
        if display.startswith("YALI1"):
            return Normalized(display, "resolved", cardinality)
        partners = tuple(sorted(display_locus(p) or p for p in self.graph.get(key, set())))
        if not partners:
            return Normalized(None, "unmapped", cardinality)
        if cardinality != "1:1":
            return Normalized(None, "unresolved", cardinality, partners)
        return Normalized(partners[0], "resolved", cardinality)

    def add_model_genes(self, genes: dict[str, dict]) -> None:
        """Rows for the model's genes; `genes` maps gene ID -> cobra-style annotation dict."""
        for gene_id, annotation in sorted(genes.items()):
            norm = self.normalize(gene_id)
            key = canonical_locus_key(norm.canonical or gene_id)
            yali0 = tuple(sorted(display_locus(p) or p for p in self.graph.get(key, set())
                                 if p.startswith("yali0")))
            self.rows[gene_id] = IdMapRow(
                model_gene=gene_id, canonical=norm.canonical, status=norm.status, cardinality=norm.cardinality,
                yali0=yali0, uniprot=_as_tuple(annotation.get("uniprot")),
                ncbigene=_as_tuple(annotation.get("ncbigene")), note=_note(gene_id, norm))
        for row in self.rows.values():
            if row.status != "resolved":        # an unmapped gene's accession must not resolve (D6)
                continue
            for accession in row.uniprot + row.ncbigene:
                self.aliases.setdefault(accession, set()).add(row.canonical)

    def model_gene_canonical(self, gene_id: str) -> str:
        """Canonical ID for a model gene, or the model ID itself when it cannot be mapped."""
        row = self.rows.get(gene_id)
        return row.canonical if row is not None and row.canonical else gene_id


def _as_tuple(value) -> tuple[str, ...]:
    if value is None:
        return ()
    return tuple(sorted(value)) if isinstance(value, (list, tuple, set)) else (str(value),)


def _note(gene_id: str, norm: Normalized) -> str:
    if norm.status == "resolved" and not gene_id.startswith("YALI1"):
        return f"model uses a non-YALI1 ID; 1:1 crosswalk to {norm.canonical}"
    if norm.status == "unresolved":
        return f"cross-assembly link is {norm.cardinality}: {', '.join(norm.candidates)}"
    if norm.status == "unmapped":
        return "no YALI1 identifier and no crosswalk link"
    return ""


def load_id_map(locus_map_dir: Path, exclusions_path: Path) -> IdMap:
    """The same sources and exclusions as `LocusCrosswalk.from_csvs`; checked against it."""
    from scripts.gem_annotate.locus_resolver import LocusCrosswalk, _read_pairs

    iyli21 = locus_map_dir / "iyli21_genes_vs_S2.csv"
    s2 = locus_map_dir / "s2_metabolic_genes.csv"
    crosswalk = LocusCrosswalk.from_csvs(iyli21, s2, exclusions_path)
    pairs = _read_pairs(iyli21, "yali1_s2", "yali0") + _read_pairs(s2, "yali1", "yali0")
    id_map = IdMap.from_pairs(pairs, crosswalk.excluded_pairs, fingerprint=crosswalk.fingerprint)
    one_to_one = {k for k in id_map.graph if id_map.cardinality(k) == "1:1"}
    if one_to_one != set(crosswalk.aliases):
        raise ValueError("id_map 1:1 links disagree with LocusCrosswalk aliases")
    return id_map
