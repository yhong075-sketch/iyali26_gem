"""Contradiction engine (spec §1.5–1.7, §5). Pure Python, no LLM.

`classify_pair` returns the L3 relation between two claims, or None when they are not
comparable (different bucket or predicate, speculative, or a transfer at taxon distance ≥ 2).
"""

import hashlib
from collections import defaultdict
from itertools import combinations
from typing import Iterable

from .schemas import (ACTIVITY_LEVELS, Claim, ClaimLevel, PairClass, PairRelation, Predicate,
                      FLAG_PREDICATES)
from .weights import Weights

CONDITION_KEYED = (Predicate.FLUX_ACTIVE_UNDER, Predicate.ESSENTIAL_UNDER)
DIRECT_LOCALIZATION = (ClaimLevel.LOCALIZATION_DIRECT,)

COMPATIBLE, CONFLICT, SCOPED = "compatible", "conflict", "scoped"


def _fractionated(claim: Claim) -> bool:
    return claim.claim_level == ClaimLevel.ACTIVITY_PURIFIED and "fraction" in claim.method.casefold()


def _direct_localization(claim: Claim) -> bool:
    return claim.claim_level in DIRECT_LOCALIZATION or _fractionated(claim)


def value_conflict(a: Claim, b: Claim) -> tuple[str, str | None]:
    """Per-predicate incompatibility (spec §1.6): compatible, conflict, or scoped with a reason."""
    opposite = a.polarity != b.polarity
    p = a.predicate
    if p in FLAG_PREDICATES:
        return (CONFLICT, None) if opposite else (COMPATIBLE, None)

    if p == Predicate.GENE_ENCODES:
        if a.value.gene_id == b.value.gene_id:
            if not opposite:
                return COMPATIBLE, None
            if a.claim_level == b.claim_level:
                return CONFLICT, None
            negative = a if a.polarity == "negative" else b
            if negative.claim_level == ClaimLevel.PHENOTYPE_KO:
                return SCOPED, "redundancy_possible"
            return SCOPED, "level_mismatch"
        both_positive = a.polarity == b.polarity == "positive"
        if both_positive and (a.value.exclusive or b.value.exclusive):
            return CONFLICT, "exclusive"
        return COMPATIBLE, None

    if p == Predicate.LOCALIZED_IN:
        overlap = bool(a.value.compartments & b.value.compartments)
        if opposite and overlap:
            clash = True
        elif not opposite and a.polarity == "positive" and not overlap:
            clash = True
        else:
            return COMPATIBLE, None
        if clash and _direct_localization(a) and _direct_localization(b):
            return CONFLICT, None
        return SCOPED, "level_mismatch"

    if p == Predicate.DIRECTIONALITY:
        da, db = a.value.direction, b.value.direction
        if opposite:
            incompatible = da == db
        elif {da, db} == {"forward", "reverse"}:
            incompatible = True
        elif da != db:   # reversible vs an irreversible value
            incompatible = a.assertion_strength == b.assertion_strength == "assertive"
        else:
            incompatible = False
        if not incompatible:
            return COMPATIBLE, None
        if a.value.level != b.value.level:
            return SCOPED, "direction_level"
        return CONFLICT, None

    if p == Predicate.COFACTOR:
        shared_roles = set(a.value.roles) & set(b.value.roles)
        if opposite:
            disjoint = not any(a.value.roles[r] & b.value.roles[r] for r in shared_roles)
            if disjoint:
                return COMPATIBLE, None
        elif not any(not (a.value.roles[r] & b.value.roles[r]) for r in shared_roles):
            return COMPATIBLE, None
        if a.claim_level in ACTIVITY_LEVELS and b.claim_level in ACTIVITY_LEVELS:
            return CONFLICT, None
        return SCOPED, "level_mismatch"

    return COMPATIBLE, None   # SUBSYSTEM_MEMBER: never incompatible


def _strain_key(strain: str | None) -> str | None:
    return None if strain is None else "".join(strain.split()).casefold()


def shares_independence_class(a: Claim, b: Claim, weights: Weights) -> bool:
    """Same source, shared corresponding author, or ≥ N shared authors within the window."""
    if a.source_id == b.source_id:
        return True
    if abs(a.year - b.year) > weights.independence_window_years:
        return False
    if a.corresponding_author and a.corresponding_author == b.corresponding_author:
        return True
    return len(set(a.authors) & set(b.authors)) >= weights.independence_min_shared_authors


def is_derivative(a: Claim, b: Claim, earliest_original: dict[tuple, str]) -> bool:
    """b repeats a (or a repeats b): citation, shared root source, or restated earliest finding."""
    if a.source_id == b.source_id:
        return True
    if a.source_id in b.cited_support or b.source_id in a.cited_support:
        return True
    for x, y in ((a, b), (b, a)):
        if not y.is_original_finding and earliest_original.get(_finding_key(y)) == x.claim_id:
            return True
    return False


def _finding_key(claim: Claim) -> tuple:
    return claim.subject.key, claim.predicate, claim.polarity


def earliest_originals(claims: Iterable[Claim]) -> dict[tuple, str]:
    best: dict[tuple, Claim] = {}
    for c in claims:
        if not c.is_original_finding:
            continue
        key = _finding_key(c)
        if key not in best or (c.year, c.claim_id) < (best[key].year, best[key].claim_id):
            best[key] = c
    return {k: v.claim_id for k, v in best.items()}


def pair_id(a: str, b: str) -> str:
    first, second = sorted((a, b))
    return "P-" + hashlib.sha256(f"{first}\0{second}".encode()).hexdigest()[:12]


def classify_pair(a: Claim, b: Claim, weights: Weights,
                  earliest_original: dict[tuple, str] | None = None) -> PairRelation | None:
    if a.subject.bucket != b.subject.bucket or a.predicate != b.predicate:
        return None
    if "speculative" in (a.assertion_strength, b.assertion_strength):
        return None
    if a.taxon_distance > 1 or b.taxon_distance > 1:
        return None
    if a.subject.status != "resolved" or b.subject.status != "resolved":
        return None
    earliest_original = earliest_original if earliest_original is not None else earliest_originals((a, b))
    shared = shares_independence_class(a, b, weights)

    def relation(cls: PairClass, reason: str | None = None) -> PairRelation:
        return PairRelation(pair_id=pair_id(a.claim_id, b.claim_id), claim_a=a.claim_id, claim_b=b.claim_id,
                            pair_class=cls, reason=reason, independence_shared=shared)

    kind, reason = value_conflict(a, b)
    if kind == COMPATIBLE:
        if is_derivative(a, b, earliest_original):
            return relation(PairClass.AGREE_DERIVATIVE)
        return relation(PairClass.AGREE_INDEPENDENT)

    if (a.strain_specific or b.strain_specific) and _strain_key(a.strain) != _strain_key(b.strain):
        return relation(PairClass.SCOPED_APART, "strain")
    if a.predicate in CONDITION_KEYED:
        cond_a, cond_b = a.condition, b.condition
        if cond_a is not None and not cond_a.overlaps(cond_b):
            return relation(PairClass.SCOPED_APART, "condition")
    band_a, band_b = weights.band(a.claim_level), weights.band(b.claim_level)
    if a.polarity != b.polarity and band_a != band_b:
        negative, positive = (a, b) if a.polarity == "negative" else (b, a)
        if weights.band(negative.claim_level) == "parts" and weights.band(positive.claim_level) == "functional":
            return relation(PairClass.ANOMALY_CROSS_LEVEL)
        return relation(PairClass.SCOPED_APART, "level_mismatch")
    if kind == SCOPED:
        return relation(PairClass.SCOPED_APART, reason)

    for earlier, later in ((a, b), (b, a)):
        if earlier.source_id in later.rebuts:
            return relation(PairClass.SUPERSEDED, f"{later.claim_id} rebuts {earlier.source_id}")
    if a.assertion_strength == b.assertion_strength == "assertive":
        return relation(PairClass.CONTRADICTION_STRONG, reason)
    return relation(PairClass.CONTRADICTION_WEAK, reason)


def classify_all(claims: list[Claim], weights: Weights) -> list[PairRelation]:
    """All L3 relations, bucket by bucket (buckets are small; O(n²) per bucket)."""
    earliest = earliest_originals(claims)
    buckets: dict[tuple, list[Claim]] = defaultdict(list)
    for c in claims:
        buckets[(c.subject.bucket, c.predicate)].append(c)
    relations = []
    for key in sorted(buckets, key=lambda k: (k[0] or "", k[1].value)):
        members = sorted(buckets[key], key=lambda c: c.claim_id)
        for a, b in combinations(members, 2):
            rel = classify_pair(a, b, weights, earliest)
            if rel is not None:
                relations.append(rel)
    return relations


def derive_existence_claims(claims: list[Claim]) -> list[Claim]:
    """A positive FLUX_ACTIVE_UNDER claim implies REACTION_EXISTS at the same level (spec §1.2)."""
    from .schemas import FlagValue, SubjectRef

    derived = []
    for c in claims:
        if c.predicate != Predicate.FLUX_ACTIVE_UNDER or c.polarity != "positive" or c.derived_from:
            continue
        derived.append(c.model_copy(update={
            "claim_id": f"{c.claim_id}~exists", "predicate": Predicate.REACTION_EXISTS,
            "value": FlagValue(predicate=Predicate.REACTION_EXISTS), "condition": None,
            "derived_from": c.claim_id,
            "subject": SubjectRef(kind="reaction", reaction_id=c.subject.reaction_id, status=c.subject.status,
                                  raw=c.subject.raw)}))
    return derived
