"""Shared builders for the evidence-pipeline tests (no LLM, no research workspace)."""

import hashlib

from evidence.schemas import (Claim, ClaimLevel, CofactorValue, Condition, DirectionalityValue, FlagValue,
                              GeneEncodesValue, LocalizedValue, Passage, Predicate, Source, SourceKind, SubjectRef)
from evidence.taxonomy import taxon_distance
from evidence.weights import load_weights

W = load_weights()
YL = "Yarrowia lipolytica"


def make_claim(cid: str, predicate: Predicate, polarity: str = "positive",
               level: ClaimLevel = ClaimLevel.ACTIVITY_PURIFIED, *, reaction: str = "R1", gene: str | None = None,
               subject: SubjectRef | None = None, strength: str = "assertive", organism: str = YL,
               strain: str | None = "W29", condition: Condition | None = None, source: str | None = None,
               year: int = 2020, original: bool = True, cited: tuple = (), rebuts: tuple = (),
               authors: tuple = (), exclusive: bool = False, compartments: set | None = None,
               direction: str = "forward", dlevel: str = "biochemical", roles: dict | None = None,
               method: str = "assay", strain_specific: bool = False) -> Claim:
    if predicate == Predicate.GENE_ENCODES:
        subject = subject or SubjectRef(kind="gene_reaction", reaction_id=reaction, gene_id=gene)
        value = GeneEncodesValue(predicate=predicate, gene_id=subject.gene_id, reaction_id=reaction,
                                 exclusive=exclusive)
    elif predicate == Predicate.LOCALIZED_IN:
        value = LocalizedValue(predicate=predicate, compartments=frozenset(compartments or {"C_cy"}))
    elif predicate == Predicate.DIRECTIONALITY:
        value = DirectionalityValue(predicate=predicate, direction=direction, level=dlevel)
    elif predicate == Predicate.COFACTOR:
        value = CofactorValue(predicate=predicate, roles={k: frozenset(v) for k, v in (roles or {}).items()})
    else:
        value = FlagValue(predicate=predicate)
    if subject is None:
        subject = (SubjectRef(kind="gene", gene_id=gene) if predicate == Predicate.ESSENTIAL_UNDER and gene
                   else SubjectRef(kind="reaction", reaction_id=reaction))
    return Claim(claim_id=cid, source_id=source or f"S-{cid}", passage_id=f"P-{cid}", subject=subject,
                 predicate=predicate, value=value, polarity=polarity, assertion_strength=strength, claim_level=level,
                 organism=organism, strain=strain, taxon_distance=taxon_distance(organism, strain, W),
                 condition=condition, strain_specific=strain_specific, method=method, is_original_finding=original,
                 cited_support=cited, rebuts=rebuts, eco=W.eco_for(level), extraction_run_id="test",
                 authors=authors, year=year)


def support_for(claims: list[Claim], kind: SourceKind = SourceKind.ARTICLE) -> tuple[dict, dict]:
    """Sources and verified passages for test claims."""
    sources, passages = {}, {}
    for c in claims:
        sources.setdefault(c.source_id, Source(source_id=c.source_id, kind=kind, title=c.source_id, year=c.year,
                                               peer_reviewed=True))
        quote = f"quote of {c.claim_id}"
        passages[c.passage_id] = Passage(passage_id=c.passage_id, source_id=c.source_id, quote=quote, section="t",
                                         char_start=0, char_end=len(quote),
                                         quote_sha256=hashlib.sha256(quote.encode()).hexdigest(), verified=True)
    return sources, passages
