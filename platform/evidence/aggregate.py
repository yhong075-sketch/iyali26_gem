"""L4 forward pass (spec §2.1): one Assertion per (subject, predicate), weights frozen in the YAML.

Sides are (polarity, value) groups of admissible claims, each with its own score. Contradictions come
from L3. When sides contradict, resolution applies R5 (assertive rebuttal) first, because it wins
"regardless of R2–R4", then R2 directness, R3 taxon distance, R4 independence and R6 recency; the
winner and every side outside the contradiction stand. Sides that are only scoped apart all stand.
R7 network consistency is advisory and never decides.
"""

import hashlib
import json
from collections import defaultdict

from .rules import shares_independence_class
from .schemas import (Assertion, Claim, ClaimLevel, CONFLICT_CLASSES, PairClass, PairRelation, Predicate, Side,
                      Source, SourceKind)
from .weights import Weights

FUNCTIONAL = "functional"


def effective_level(claim: Claim) -> tuple[ClaimLevel, ClaimLevel | None]:
    """Claims at distance ≥ 2 enter as transfers: ORTHOLOGY with the original level as transfer_from."""
    if claim.taxon_distance >= 2 and claim.claim_level != ClaimLevel.ORTHOLOGY:
        return ClaimLevel.ORTHOLOGY, claim.claim_level
    return claim.claim_level, claim.transfer_from_level


def admissible(claim: Claim, weights: Weights) -> bool:
    if claim.assertion_strength == "speculative" or claim.subject.status != "resolved":
        return False
    level, _ = effective_level(claim)
    if not weights.admissible(claim.predicate_key, level):
        return False
    if claim.predicate == Predicate.LOCALIZED_IN and level == ClaimLevel.ACTIVITY_PURIFIED:
        return "fraction" in claim.method.casefold()
    return True


def value_key(claim: Claim) -> str:
    v = claim.value
    if claim.predicate == Predicate.LOCALIZED_IN:
        return ",".join(sorted(v.compartments))
    if claim.predicate == Predicate.DIRECTIONALITY:
        return f"{v.direction}@{v.level}"
    if claim.predicate == Predicate.COFACTOR:
        return json.dumps({role: sorted(ids) for role, ids in sorted(v.roles.items())})
    if claim.predicate == Predicate.GENE_ENCODES:
        return "exclusive" if v.exclusive else ""
    return ""


def assertion_id(subject_key: str, predicate: Predicate) -> str:
    return "A-" + hashlib.sha256(f"{subject_key}\0{predicate.value}".encode()).hexdigest()[:12]


def _n_independent(claims: list[Claim], weights: Weights) -> int:
    """Distinct independence classes among original findings; restatements are excluded (§1.7).
    A derived existence claim counts with its flux source, so 13C support is one class, not zero."""
    originals = [c for c in claims if c.is_original_finding]
    parent = {c.claim_id: c.claim_id for c in originals}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(originals):
        for b in originals[i + 1:]:
            if shares_independence_class(a, b, weights) or a.source_id in b.cited_support \
                    or b.source_id in a.cited_support:
                parent[find(a.claim_id)] = find(b.claim_id)
    return len({find(c.claim_id) for c in originals})


def _rank(claim: Claim, weights: Weights) -> int:
    return weights.rank(claim.predicate_key, effective_level(claim)[0])


def _score(claims: list[Claim], weights: Weights) -> int:
    scores = []
    for c in claims:
        level, transfer = effective_level(c)
        scores.append(weights.claim_score(c.predicate_key, level, c.taxon_distance, transfer))
    return max(scores, default=0)


def _side(claims: list[Claim], weights: Weights) -> Side:
    best = min(claims, key=lambda c: (_rank(c, weights), c.claim_id))
    return Side(polarity=best.polarity, value_key=value_key(best), claim_ids=tuple(sorted(c.claim_id for c in claims)),
                best_level=effective_level(best)[0], best_rank=_rank(best, weights), score=_score(claims, weights),
                n_independent=_n_independent(claims, weights),
                min_taxon_distance=min(c.taxon_distance for c in claims),
                latest_year=max(c.year for c in claims))


def _strength(rel: PairRelation, claims_by_id: dict[str, Claim]) -> str:
    """SUPERSEDED pairs are strong only when both claims are assertive (§1.7)."""
    if rel.pair_class == PairClass.SUPERSEDED:
        both = all(claims_by_id[c].assertion_strength == "assertive" for c in (rel.claim_a, rel.claim_b))
        return "CONTRADICTION_STRONG" if both else "CONTRADICTION_WEAK"
    return {PairClass.CONTRADICTION_STRONG: "CONTRADICTION_STRONG", PairClass.CONTRADICTION_WEAK:
            "CONTRADICTION_WEAK", PairClass.ANOMALY_CROSS_LEVEL: "ANOMALY"}[rel.pair_class]


def _conflict_status(rels: list[PairRelation], claims_by_id: dict[str, Claim]) -> str:
    found = {_strength(r, claims_by_id) for r in rels}
    for status in ("CONTRADICTION_STRONG", "CONTRADICTION_WEAK", "ANOMALY"):
        if status in found:
            return status
    return "NONE"


def _resolve(contenders: list[int], sides: list[Side], superseded: list[PairRelation],
             claims_by_id: dict[str, Claim], side_of: dict[str, int]) -> tuple[int | None, str | None, list]:
    trace: list[tuple[str, str]] = []
    rebut_sides = set()
    for rel in superseded:
        rebutter = (rel.reason or "").split(" rebuts ")[0]
        if rebutter in side_of and claims_by_id[rebutter].assertion_strength == "assertive":
            rebut_sides.add(side_of[rebutter])
    if len(rebut_sides) == 1:
        winner = rebut_sides.pop()
        trace.append(("R5", f"assertive rebuttal resolves toward side {winner}"))
        return winner, "R5", trace
    trace.append(("R5", "no single assertive rebuttal" if not rebut_sides else "rebuttals point both ways"))

    remaining = list(contenders)
    steps = (("R2", lambda s: -s.best_rank, "best directness"),
             ("R3", lambda s: -s.min_taxon_distance, "lowest taxon distance"),
             ("R4", lambda s: s.n_independent, "most independent classes"),
             ("R6", lambda s: s.latest_year, "most recent"))
    for rule, key, label in steps:
        best = max(key(sides[i]) for i in remaining)
        remaining = [i for i in remaining if key(sides[i]) == best]
        if len(remaining) == 1:
            trace.append((rule, f"{label}: side {remaining[0]} wins"))
            return remaining[0], rule, trace
        trace.append((rule, f"{label}: tie between sides {remaining}"))
    return None, None, trace


def _lead(indices: list[int], sides: list[Side]) -> int | None:
    """The standing side reported beside the score: highest score, positive first."""
    return min(indices, key=lambda i: (-sides[i].score, sides[i].polarity != "positive", i), default=None)


def aggregate(claims: list[Claim], relations: list[PairRelation], sources: dict[str, Source],
              weights: Weights) -> list[Assertion]:
    claims_by_id = {c.claim_id: c for c in claims}
    admissible_ids = {c.claim_id for c in claims if admissible(c, weights)}
    groups: dict[tuple[str, Predicate], list[Claim]] = defaultdict(list)
    for c in claims:
        groups[(c.subject.key, c.predicate)].append(c)
    relations_by_claim: dict[str, list[PairRelation]] = defaultdict(list)
    for rel in relations:
        relations_by_claim[rel.claim_a].append(rel)
        relations_by_claim[rel.claim_b].append(rel)

    assertions = []
    for (subject_key, predicate), members in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1].value)):
        members = sorted(members, key=lambda c: c.claim_id)
        kept = [c for c in members if c.claim_id in admissible_ids]
        excluded = tuple(c.claim_id for c in members if c.claim_id not in admissible_ids)
        by_side: dict[tuple[str, str], list[Claim]] = defaultdict(list)
        for c in kept:
            by_side[(c.polarity, value_key(c))].append(c)
        sides = [_side(by_side[k], weights) for k in sorted(by_side)]
        side_of = {cid: i for i, s in enumerate(sides) for cid in s.claim_ids}

        member_ids = {c.claim_id for c in members}
        rels = {rel.pair_id: rel for c in members for rel in relations_by_claim[c.claim_id]}
        own = [rel for rel in rels.values() if rel.claim_a in member_ids and rel.claim_b in member_ids]
        conflicts = [rel for rel in own if rel.pair_class in CONFLICT_CLASSES
                     and rel.claim_a in side_of and rel.claim_b in side_of
                     and side_of[rel.claim_a] != side_of[rel.claim_b]]
        # Exclusive GENE_ENCODES conflicts cross subjects (other genes of the same reaction); both
        # claims must be admissible to count.
        cross = [rel for rel in rels.values() if rel not in own and rel.pair_class in CONFLICT_CLASSES
                 and rel.claim_a in admissible_ids and rel.claim_b in admissible_ids
                 and (rel.claim_a in side_of or rel.claim_b in side_of)]
        status = _conflict_status(conflicts + cross, claims_by_id)
        n_scoped = sum(1 for rel in own if rel.pair_class == PairClass.SCOPED_APART)
        trace: list[tuple[str, str]] = [("R1", f"{n_scoped} pair(s) scoped apart")]

        winner, resolved_by, resolved = None, None, True
        if not sides:
            standing: list[int] = []
            trace.append(("R0", "no admissible claims"))
        elif not conflicts:
            standing = list(range(len(sides)))
            trace.append(("R0", "no contradiction between sides; all sides stand"))
        else:
            contenders = sorted({side_of[r.claim_a] for r in conflicts} | {side_of[r.claim_b] for r in conflicts})
            superseded = [r for r in conflicts if r.pair_class == PairClass.SUPERSEDED]
            winner, resolved_by, steps = _resolve(contenders, sides, superseded, claims_by_id, side_of)
            trace += steps
            outside = [i for i in range(len(sides)) if i not in contenders]
            if winner is None:
                resolved, standing = False, list(range(len(sides)))
                trace.append(("UNRESOLVED", "R2–R6 leave a tie"))
            else:
                standing = sorted([winner] + outside)
        if cross:
            trace.append(("X", f"{len(cross)} exclusive conflict(s) with other genes of this reaction"))
        trace.append(("R7", "network consistency not computed in M1 (advisory only)"))

        lead = winner if winner is not None else _lead(standing, sides)   # a resolved contradiction scores from its winner
        score = sides[lead].score if lead is not None else 0
        flag = None
        if not resolved:
            score, flag = min(score, weights.unresolved_cap), "CONFLICTING"
        leaves = [claims_by_id[cid] for i in standing for cid in sides[i].claim_ids]
        homology_only = bool(kept) and not any(
            c.taxon_distance <= 1 and weights.band(c.claim_level) == FUNCTIONAL for c in kept)
        inherited = None
        if kept and all(sources[c.source_id].kind == SourceKind.PRIOR_MODEL for c in kept):
            parents = sorted({sources[c.source_id].properties.get("model_name", c.source_id) for c in kept
                              if sources[c.source_id].properties.get("role") == "parent"})
            inherited = ";".join(parents) or None

        assertions.append(Assertion(
            assertion_id=assertion_id(subject_key, predicate), subject_key=subject_key, predicate=predicate,
            sides=tuple(sides), excluded_claim_ids=excluded,
            pair_ids=tuple(sorted(r.pair_id for r in own + cross)), conflict_status=status,
            resolution_trace=tuple(trace), resolved=resolved, winner=winner, standing=tuple(standing),
            resolved_by=resolved_by, score=score, flag=flag,
            n_independent=sides[lead].n_independent if lead is not None else 0,
            max_taxon_distance_used=max((c.taxon_distance for c in leaves), default=None),
            homology_only=homology_only, inherited_from=inherited, weights_version=weights.version))
    return assertions
