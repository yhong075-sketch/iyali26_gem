"""Deterministic L0–L2 builders (no LLM). Every claim gets a passage verified against its source text.

Converters: parent-model inclusion (iYli21), the current model's evidence notes, W29 13C flux rows,
the essentiality screens and consensus list, the legacy essentiality cases (sources only), and
hand-written literature claims checked against the cached gem_review chunks.
"""

import csv
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .id_map import IdMap
from .schemas import (Claim, ClaimLevel, Condition, FlagValue, GeneEncodesValue, Passage, Predicate, RejectedClaim,
                      Source, SourceKind, SourceTier, SubjectRef, SubsystemValue, ClaimValue)
from .taxonomy import taxon_distance
from .weights import Weights

MAX_WORDS = 25
YL = "Yarrowia lipolytica"


@dataclass
class Output:
    sources: dict[str, Source] = field(default_factory=dict)
    passages: dict[str, Passage] = field(default_factory=dict)
    claims: list[Claim] = field(default_factory=list)
    rejected: list[RejectedClaim] = field(default_factory=list)

    def merge(self, other: "Output") -> "Output":
        for sid, s in other.sources.items():
            if sid in self.sources and self.sources[sid] != s:
                raise ValueError(f"source {sid} defined twice with different content")
            self.sources[sid] = s
        self.passages.update(other.passages)
        self.claims += other.claims
        self.rejected += other.rejected
        return self


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _digest(*parts) -> str:
    return hashlib.sha256("\0".join(str(p) for p in parts).encode()).hexdigest()[:12]


def find_quote(text: str, quote: str) -> tuple[int, int] | None:
    """Offsets of a whitespace-normalised exact match of `quote` in `text`."""
    tokens = quote.split()
    if not tokens:
        return None
    match = re.search(r"\s+".join(re.escape(t) for t in tokens), text)
    return (match.start(), match.end()) if match else None


def first_words(text: str, n: int = MAX_WORDS) -> str:
    return " ".join(text.split()[:n])


class Builder:
    """Shared helpers: passages, subject normalisation and claim construction."""

    def __init__(self, weights: Weights, id_map: IdMap, run_id: str):
        self.weights, self.id_map, self.run_id = weights, id_map, run_id
        self.out = Output()

    def reject(self, source_id: str | None, reason: str, detail: str) -> None:
        """One row per rejection; the run ID and position keep identical details apart."""
        self.out.rejected.append(RejectedClaim(
            reject_id="X-" + _digest(self.run_id, len(self.out.rejected), source_id, reason, detail),
            source_id=source_id, reason=reason, detail=detail))

    def passage(self, source: Source, text: str, quote: str, section: str | None,
                offset_hint: int = 0) -> Passage | None:
        """A verified passage, or None with a rejected row. Over-long quotes are rejected, never cut."""
        if len(quote.split()) > MAX_WORDS:
            self.reject(source.source_id, "quote_too_long", f"{section}: {len(quote.split())} words: {quote}")
            return None
        span = find_quote(text[offset_hint:], quote)
        if span is None:
            self.reject(source.source_id, "quote_not_found", f"{section}: {quote}")
            return None
        start, end = span[0] + offset_hint, span[1] + offset_hint
        p = Passage(passage_id="L1-" + _digest(source.source_id, section, start, end), source_id=source.source_id,
                    quote=quote, section=section, char_start=start, char_end=end,
                    quote_sha256=sha256_text(" ".join(quote.split())), verified=True)
        self.out.passages[p.passage_id] = p
        return p

    def gene_subject(self, raw_gene: str, reaction_id: str | None) -> SubjectRef:
        norm = self.id_map.normalize(raw_gene)
        kind = "gene_reaction" if reaction_id else "gene"
        return SubjectRef(kind=kind, reaction_id=reaction_id, gene_id=norm.canonical or raw_gene,
                          status=norm.status, raw=raw_gene, candidate_subjects=norm.candidates)

    def claim(self, source: Source, passage: Passage, subject: SubjectRef, predicate: Predicate,
              value: ClaimValue, polarity: str, level: ClaimLevel, *, organism: str = YL,
              strain: str | None = None, strength: str = "assertive", method: str,
              original: bool = True, condition: Condition | None = None, cited_support: tuple = (),
              rebuts: tuple = (), transfer_from: ClaimLevel | None = None, justification: str | None = None,
              strain_specific: bool = False, year: int | None = None) -> Claim:
        c = Claim(
            claim_id="pending", source_id=source.source_id, passage_id=passage.passage_id, subject=subject, predicate=predicate,
            value=value, polarity=polarity, assertion_strength=strength, claim_level=level,
            level_justification=justification, transfer_from_level=transfer_from, organism=organism, strain=strain,
            taxon_distance=taxon_distance(organism, strain, self.weights), condition=condition,
            strain_specific=strain_specific, method=method, is_original_finding=original,
            cited_support=cited_support, rebuts=rebuts, eco=self.weights.eco_for(level),
            extraction_run_id=self.run_id, authors=source.authors, corresponding_author=source.corresponding_author,
            year=year if year is not None else (source.year or 0))
        # The ID digests every field, so two claims differing in any field never collide.
        content = json.dumps(c.model_dump(mode="json", exclude={"claim_id"}), sort_keys=True)
        c = c.model_copy(update={"claim_id": "C-" + hashlib.sha256(content.encode()).hexdigest()[:16]})
        self.out.claims.append(c)
        return c


# ── Models (parent inclusion and the current model's notes) ──────────────────────────────────

def _reaction_span(xml: str, reaction_id: str) -> tuple[int, int] | None:
    for raw_id in (reaction_id, f"R_{reaction_id}"):
        match = re.search(r'<reaction [^>]*\bid="%s"[^>]*>.*?</reaction>' % re.escape(raw_id), xml, re.S)
        if match:
            return match.start(), match.end()
    return None


def _species_set(element: str) -> set[str]:
    """Species of a reaction element as model IDs (iYli21 escapes m570[C_lp] as m570__91__C_lp__93__)."""
    ids = re.findall(r'<speciesReference [^>]*species="([^"]+)"', element)
    return {re.sub(r"^M_", "", i).replace("__91__", "[").replace("__93__", "]") for i in ids}


def parent_model_claims(xml_path: Path, reaction_ids: list[str], weights: Weights, id_map: IdMap,
                        lineage: dict, model_species: dict[str, set[str]],
                        protons: frozenset[str] = frozenset()) -> Output:
    """iYli21 inclusion: reaction present, its GPR genes, its groups (MODEL_INCLUSION, D8).
    A reaction counts as the same only when its species (substrates, products, compartments)
    match the current model's; an ID match alone is not enough (AGENTS.md §5). Protons are left
    out of the match because iYali26 rebalanced protonation after iYli21."""
    raw = xml_path.read_bytes()
    xml = raw.decode("utf-8")
    sha = hashlib.sha256(raw).hexdigest()
    info = lineage["iYli21"]
    b = Builder(weights, id_map, "det-parent-model-v1")
    source = Source(source_id=f"model:iYli21:{sha[:12]}", kind=SourceKind.PRIOR_MODEL, title=info["title"],
                    year=info["year"], tier=SourceTier.D, peer_reviewed=True, organisms=(YL,),
                    parent_model=info["parent"], text_ref=xml_path.name, text_sha256=sha,
                    properties={"model_name": "iYli21", "role": "parent", "publication": info["publication"],
                                "origin_tags": info["origin_tags"]})
    b.out.sources[source.source_id] = source
    groups = [(m.start(), m.group(1), xml[m.start():xml.find("</groups:group>", m.start())])
              for m in re.finditer(r'<groups:group [^>]*groups:name="([^"]*)"[^>]*>', xml)]
    for rid in reaction_ids:
        span = _reaction_span(xml, rid)
        if span is None:
            b.reject(source.source_id, "not_in_parent", f"{rid} absent from iYli21")
            continue
        element = xml[span[0]:span[1]]
        if _species_set(element) - protons != model_species[rid] - protons:
            differ = sorted((_species_set(element) ^ model_species[rid]) - protons)
            b.reject(source.source_id, "parent_reaction_differs", f"{rid}: species differ {differ}")
            continue
        name = re.search(r'\bname="([^"]*)"', element)
        p = b.passage(source, xml, f'name="{name.group(1)}"', f"reaction {rid}", span[0])
        if p:
            b.claim(source, p, SubjectRef(kind="reaction", reaction_id=rid), Predicate.REACTION_EXISTS,
                    FlagValue(predicate=Predicate.REACTION_EXISTS), "positive", ClaimLevel.MODEL_INCLUSION,
                    strain="W29", method="present in parent GEM iYli21")
        for gm in re.finditer(r'fbc:geneProduct="G_([^"]+)"', element):
            gp = b.passage(source, xml, gm.group(0), f"reaction {rid} GPR", span[0] + gm.start())
            if gp:
                subject = b.gene_subject(gm.group(1), rid)
                b.claim(source, gp, subject, Predicate.GENE_ENCODES,
                        GeneEncodesValue(predicate=Predicate.GENE_ENCODES, gene_id=subject.gene_id, reaction_id=rid),
                        "positive", ClaimLevel.MODEL_INCLUSION, strain="W29", method="GPR in parent GEM iYli21")
        for start, gname, body in groups:
            if f'groups:idRef="{rid}"' in body:
                gp = b.passage(source, xml, f'groups:name="{gname}"', f"group {gname}", start)
                if gp:
                    b.claim(source, gp, SubjectRef(kind="reaction", reaction_id=rid), Predicate.SUBSYSTEM_MEMBER,
                            SubsystemValue(predicate=Predicate.SUBSYSTEM_MEMBER, pathways=frozenset({gname}),
                                           ontology="model_groups"),
                            "positive", ClaimLevel.MODEL_INCLUSION, strain="W29", method="group in parent GEM iYli21")
    return b.out


def lineage_sources(lineage: dict) -> Output:
    """The ancestors above iYli21, as L0 sources without claims (D8 parent_model chain)."""
    out = Output()
    for name, info in lineage.items():
        if name in ("iYli21", "iYali26"):
            continue
        out.sources[f"model:{name}"] = Source(
            source_id=f"model:{name}", kind=SourceKind.PRIOR_MODEL, title=info["title"], year=info.get("year"),
            tier=SourceTier.D, peer_reviewed=True, parent_model=info.get("parent"),
            properties={"model_name": name, "role": "ancestor", "publication": info.get("publication", "")})
    return out


NOTE_RE = re.compile(r"<p>([^<:]+):\s*([^<]*)</p>")
NOTE_KEY_RE = re.compile(r"(^|_)evidence(_status|_tier)?$")   # evidence-status keys, not limits or counterevidence


def model_note_claims(xml_path: Path, model_commit: str, reaction_ids: list[str], weights: Weights,
                      id_map: IdMap, lineage: dict) -> Output:
    """The model's own evidence-status notes as MODEL_INCLUSION (kind prior_model, role self; D8)."""
    raw = xml_path.read_bytes()
    xml = raw.decode("utf-8")
    info = lineage["iYali26"]
    b = Builder(weights, id_map, "det-model-notes-v1")
    source = Source(source_id=f"model:iYali26:{model_commit[:12]}", kind=SourceKind.PRIOR_MODEL,
                    title=f"iYali26 {xml_path.name}", year=info["year"], tier=SourceTier.D, peer_reviewed=False,
                    organisms=(YL,), parent_model=info["parent"], text_ref=xml_path.name,
                    text_sha256=hashlib.sha256(raw).hexdigest(),
                    properties={"model_name": "iYali26", "role": "self", "commit": model_commit})
    b.out.sources[source.source_id] = source
    for rid in reaction_ids:
        span = _reaction_span(xml, rid)
        if span is None:
            b.reject(source.source_id, "not_in_model", f"{rid} absent from {xml_path.name}")
            continue
        element = xml[span[0]:span[1]]
        for m in NOTE_RE.finditer(element):
            key = m.group(1).strip()
            if not NOTE_KEY_RE.search(key) or "counter" in key or "conflicting" in key:
                continue
            excerpt = first_words(f"{key}: {m.group(2)}")       # deliberate excerpt of a long note
            p = b.passage(source, xml, excerpt, f"reaction {rid} note {key}", span[0] + m.start())
            if p:
                b.claim(source, p, SubjectRef(kind="reaction", reaction_id=rid), Predicate.REACTION_EXISTS,
                        FlagValue(predicate=Predicate.REACTION_EXISTS), "positive", ClaimLevel.MODEL_INCLUSION,
                        strain="W29", method=f"model note {key}")
    return b.out


# ── Datasets ─────────────────────────────────────────────────────────────────────────────────

CARBON_SOURCES = {"glucose": "glucose", "glycerol": "glycerol", "oleic_acid": "oleic acid"}


def _num(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def w29_13c_claims(csv_path: Path, map_path: Path, weights: Weights, id_map: IdMap) -> Output:
    """FLUX_ACTIVE_UNDER from W29 13C-MFA intervals: positive when [LB, UB] excludes 0, negative when
    both bounds are 0, otherwise no claim (logged). Unpublished internal data (D9)."""
    raw = csv_path.read_bytes()
    text = raw.decode("utf-8")
    b = Builder(weights, id_map, "det-w29-13c-v1")
    source = Source(source_id=f"dataset:w29_13c:{hashlib.sha256(raw).hexdigest()[:12]}", kind=SourceKind.DATASET,
                    title="W29 13C-MFA central fluxes on glucose, glycerol and oleic acid (internal)", year=2024,
                    tier=SourceTier.E, peer_reviewed=False, unpublished=True, organisms=(YL,),
                    text_ref=csv_path.name, text_sha256=hashlib.sha256(raw).hexdigest(),
                    properties={"interval_rule": "positive if LB > 0 or UB < 0; negative if LB == UB == 0",
                                "attribution": "README says unpublished; benchmark.md cites Worland et al. 2024 (inferred)"})
    b.out.sources[source.source_id] = source
    with map_path.open(newline="") as handle:
        mapping = {row["row_id"]: row["model_reaction_id"] for row in csv.DictReader(handle)}
    lines = text.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1] + len(line))
    reader = list(csv.DictReader(lines))
    for missing in sorted(set(mapping) - {row["ID"] for row in reader}):
        b.reject(source.source_id, "map_row_missing", f"{missing} ({mapping[missing]}) not in {csv_path.name}")
    for index, row in enumerate(reader, start=1):
        rid = mapping.get(row["ID"])
        if rid is None:
            continue
        for column, carbon in CARBON_SOURCES.items():
            raw_lb, raw_ub = row.get(f"{column}_LB"), row.get(f"{column}_UB")
            lb, ub = _num(raw_lb), _num(raw_ub)
            if lb is None or ub is None:
                if raw_lb or raw_ub:
                    b.reject(source.source_id, "bounds_not_numeric", f"{row['ID']} {carbon}: [{raw_lb}, {raw_ub}]")
                continue
            if lb > 0 or ub < 0:
                polarity = "positive"
            elif lb == 0 and ub == 0:
                polarity = "negative"
            else:
                b.reject(source.source_id, "interval_includes_zero", f"{row['ID']} ({rid}) {carbon}: [{lb}, {ub}]")
                continue
            p = b.passage(source, text, f"{row['ID']},{row['Equation']}", f"row {index}", offsets[index])
            if p:
                b.claim(source, p, SubjectRef(kind="reaction", reaction_id=rid), Predicate.FLUX_ACTIVE_UNDER,
                        FlagValue(predicate=Predicate.FLUX_ACTIVE_UNDER), polarity, ClaimLevel.FLUX_13C, strain="W29",
                        method=f"13C-MFA interval [{lb}, {ub}]", condition=Condition(carbon_source=carbon))
    return b.out


def essentiality_claims(fitness_csv: Path, consensus_csv: Path, genes: set[str], weights: Weights,
                        id_map: IdMap) -> Output:
    """ESSENTIAL_UNDER per screen row (Cas9, Cas12a) and the consensus list as a derivative source."""
    b = Builder(weights, id_map, "det-essentiality-v1")
    raw = fitness_csv.read_bytes()
    text = raw.decode("utf-8")
    rows = list(csv.DictReader(text.splitlines()))
    cutoffs = {}
    for assay in sorted({r["assay"] for r in rows}):
        ess = [float(r["fitness_score"]) for r in rows if r["assay"] == assay and r["experimental_call"] == "essential"]
        non = [float(r["fitness_score"]) for r in rows if r["assay"] == assay and r["experimental_call"] != "essential"]
        cutoffs[f"observed_boundary_{assay}"] = f"essential FS <= {max(ess):.3f}; non-essential FS >= {min(non):.3f}"
    screens = Source(source_id="dataset:ramesh2023_screens", kind=SourceKind.DATASET,
                     title="Cas9 and Cas12a genome-wide screens (Ramesh et al. 2023, Commun Biol)", year=2023,
                     doi="10.1038/s42003-023-04996-8", tier=SourceTier.E, peer_reviewed=True, organisms=(YL,),
                     text_ref=fitness_csv.name, text_sha256=hashlib.sha256(raw).hexdigest(),
                     properties={**cutoffs, "experimental_cutoff_note": "boundaries observed in the normalised file; "
                                 "the paper's own calling rule not re-read", "medium": "not recorded in the file"})
    b.out.sources[screens.source_id] = screens
    for row in rows:
        gene = id_map.normalize(row["gene_id"])
        if (gene.canonical or row["gene_id"]) not in genes:
            continue
        p = b.passage(screens, text, f"{row['gene_id']},{row['source_gene_id']},{row['assay']}", f"{row['assay']} row")
        if p:
            b.claim(screens, p, b.gene_subject(row["gene_id"], None), Predicate.ESSENTIAL_UNDER,
                    FlagValue(predicate=Predicate.ESSENTIAL_UNDER),
                    "positive" if row["experimental_call"] == "essential" else "negative", ClaimLevel.PHENOTYPE_KO,
                    method=f"CRISPR-{row['assay']} pooled screen, fitness score {row['fitness_score']}")
    raw_c = consensus_csv.read_bytes()
    text_c = raw_c.decode("utf-8")
    consensus = Source(source_id="dataset:ramesh2023_consensus", kind=SourceKind.DATASET,
                       title="Consensus essential genes, essential in at least 2 of 3 screens (Ramesh et al. 2023)",
                       year=2023, doi="10.1038/s42003-023-04996-8", tier=SourceTier.E, peer_reviewed=True,
                       organisms=(YL,), text_ref=consensus_csv.name, text_sha256=hashlib.sha256(raw_c).hexdigest(),
                       properties={"experimental_cutoff": "essential in >= 2 of 3 screens; positive-only list"})
    b.out.sources[consensus.source_id] = consensus
    for row in csv.DictReader(text_c.splitlines()):
        gene = id_map.normalize(row["gene_id"])
        if (gene.canonical or row["gene_id"]) not in genes:
            continue
        p = b.passage(consensus, text_c, f"{row['gene_id']},{row['source_gene_id']}", "consensus row")
        if p:
            b.claim(consensus, p, b.gene_subject(row["gene_id"], None), Predicate.ESSENTIAL_UNDER,
                    FlagValue(predicate=Predicate.ESSENTIAL_UNDER), "positive", ClaimLevel.PHENOTYPE_KO,
                    method="consensus of three screens", original=False, cited_support=(screens.source_id,))
    return b.out


def legacy_case_sources(ledger_csv: Path, evidence_dir: Path) -> Output:
    """The 240 essentiality FN cases as L0 sources only (legacy_scoring=True; D7)."""
    out = Output()
    with ledger_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            dossier = evidence_dir / f"{row['case_id']}.json"
            data = json.loads(dossier.read_text()) if dossier.exists() else {}
            lit = data.get("literature_review") or {}
            skeptic = data.get("adversarial_review") or {}
            out.sources[row["case_id"]] = Source(
                source_id=row["case_id"], kind=SourceKind.LEGACY_CASE, legacy_scoring=True, peer_reviewed=False,
                title=data.get("claim_under_review") or row["category"], year=2026,
                subjects=tuple(x for x in (row["reaction_ids"] + ";" + row["gene_ids"]).split(";") if x),
                text_ref=dossier.name if dossier.exists() else None,
                text_sha256=hashlib.sha256(dossier.read_bytes()).hexdigest() if dossier.exists() else None,
                properties={"status": row["status"], "category": row["category"],
                            "model_sha256": row["model_sha256"], "literature_verdict": str(lit.get("verdict") or ""),
                            "skeptic_verdict": str(skeptic.get("verdict") or ""),
                            "skeptic_corrected_verdict": str(skeptic.get("corrected_verdict") or "")})
    return out


# ── Hand-written literature claims (M1) ──────────────────────────────────────────────────────

def europe_pmc_index(cache_dir: Path) -> dict[str, dict]:
    """PMID -> Europe PMC search metadata from the gem_review disk cache (offline)."""
    index: dict[str, dict] = {}
    for path in sorted(cache_dir.rglob("*.json")):
        body = json.loads(path.read_text()).get("body")
        if isinstance(body, str):
            try:
                body = json.loads(body)
            except json.JSONDecodeError:
                continue
        if isinstance(body, dict) and "resultList" in body:
            for r in body["resultList"].get("result", []):
                if r.get("pmid") and r["pmid"] not in index:
                    index[r["pmid"]] = r
    return index


def handwritten_claims(spec_path: Path, retrieval_dir: Path, cache_dir: Path, weights: Weights,
                       id_map: IdMap) -> Output:
    spec = yaml.safe_load(spec_path.read_text())
    b = Builder(weights, id_map, spec["extraction_run_id"])
    meta = europe_pmc_index(cache_dir)
    for entry in spec["claims"]:
        record = json.loads((retrieval_dir / entry["retrieval"]).read_text())
        chunk = next((c for c in record["chunks"] if c["chunk_id"] == entry["chunk_id"]), None)
        m = meta.get(str(entry["pmid"]), {})
        source = Source(source_id=f"pmid:{entry['pmid']}", kind=SourceKind.ARTICLE, title=m.get("title", ""),
                        year=int(m["pubYear"]) if m.get("pubYear") else None, doi=m.get("doi"), pmid=str(entry["pmid"]),
                        tier=SourceTier.A, peer_reviewed=True,
                        organisms=(YL,), authors=tuple(a.strip() for a in (m.get("authorString") or "").rstrip(".").split(",")
                                                        if a.strip()),
                        text_ref=f"gem_review retrieval {entry['retrieval']}")
        b.out.sources.setdefault(source.source_id, source)
        source = b.out.sources[source.source_id]
        if chunk is None:
            b.reject(source.source_id, "chunk_not_found", entry["chunk_id"])
            continue
        p = b.passage(source, chunk["text"], entry["quote"], f"{chunk['source']} {chunk['chunk_id']}")
        if p is None:
            continue
        predicate = Predicate(entry["predicate"])
        rid = entry["reaction_id"]
        if predicate == Predicate.GENE_ENCODES:
            subject = b.gene_subject(entry["gene"], rid)
            value = GeneEncodesValue(predicate=predicate, gene_id=subject.gene_id, reaction_id=rid,
                                     exclusive=entry.get("exclusive", False))
        else:
            subject = SubjectRef(kind="reaction", reaction_id=rid)
            value = FlagValue(predicate=predicate)
        b.claim(source, p, subject, predicate, value, entry["polarity"], ClaimLevel(entry["claim_level"]),
                organism=entry.get("organism", YL), strain=entry.get("strain"), strength=entry["assertion_strength"],
                method=entry["method"], original=entry.get("is_original_finding", True),
                transfer_from=ClaimLevel(entry["transfer_from_level"]) if entry.get("transfer_from_level") else None,
                justification=entry.get("level_justification"))
    return b.out
