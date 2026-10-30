"""Forward pass L2 → L5 and the backward pass (spec §2.3): trace, load_bearing, chain_quality."""

from collections import defaultdict
from dataclasses import dataclass, field

from .aggregate import aggregate
from .rules import classify_all, derive_existence_claims
from .schemas import (Assertion, Claim, ModelElementState, PairRelation, Passage, Predicate, RejectedClaim, Source,
                      Verdict)
from .verdicts import reaction_verdicts
from .weights import Weights


@dataclass(frozen=True)
class ModelIdentity:
    path: str
    commit: str
    sha256: str


@dataclass
class Bundle:
    sources: dict[str, Source]
    passages: dict[str, Passage]
    claims: dict[str, Claim]                       # includes derived claims
    rejected: list[RejectedClaim]
    relations: dict[str, PairRelation] = field(default_factory=dict)
    assertions: dict[str, Assertion] = field(default_factory=dict)
    verdicts: dict[str, Verdict] = field(default_factory=dict)
    states: dict[str, ModelElementState] = field(default_factory=dict)
    load_bearing: dict[str, tuple[str, ...]] = field(default_factory=dict)


def _subsystem_claims(claims: list[Claim]) -> dict[str, list[tuple[str, str]]]:
    out: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for c in claims:
        if c.predicate == Predicate.SUBSYSTEM_MEMBER and c.polarity == "positive" and c.assertion_strength != "speculative":
            out[c.subject.reaction_id] += [(c.value.ontology, p) for p in sorted(c.value.pathways)]
    return {k: sorted(set(v)) for k, v in out.items()}


def forward(claims: list[Claim], sources: dict[str, Source], states: dict[str, ModelElementState],
            weights: Weights, model: ModelIdentity) -> tuple[list[PairRelation], list[Assertion], list[Verdict]]:
    relations = classify_all(claims, weights)
    assertions = aggregate(claims, relations, sources, weights)
    by_key = {(a.subject_key, a.predicate): a for a in assertions}
    subsystems = _subsystem_claims(claims)
    verdicts = []
    for rid in sorted(states):
        verdicts += reaction_verdicts(states[rid], by_key, weights, model_path=model.path,
                                      model_commit=model.commit, model_sha256=model.sha256,
                                      subsystem_claims=subsystems)
    return relations, assertions, verdicts


def run(sources: dict[str, Source], passages: dict[str, Passage], claims: list[Claim],
        rejected: list[RejectedClaim], states: dict[str, ModelElementState], weights: Weights,
        model: ModelIdentity, *, with_load_bearing: bool = True) -> Bundle:
    unverified = [c for c in claims if not passages.get(c.passage_id) or not passages[c.passage_id].verified]
    if unverified:
        raise ValueError(f"claims without a verified passage: {[c.claim_id for c in unverified][:5]}")
    all_claims = sorted(claims + derive_existence_claims(claims), key=lambda c: c.claim_id)
    for label, ids in (("claim", [c.claim_id for c in all_claims]), ("rejection", [r.reject_id for r in rejected])):
        if len(ids) != len(set(ids)):
            raise ValueError(f"duplicate {label} IDs: {sorted({i for i in ids if ids.count(i) > 1})[:5]}")
    relations, assertions, verdicts = forward(all_claims, sources, states, weights, model)
    bundle = Bundle(sources=sources, passages=passages, claims={c.claim_id: c for c in all_claims},
                    rejected=sorted(rejected, key=lambda r: r.reject_id),
                    relations={r.pair_id: r for r in relations}, assertions={a.assertion_id: a for a in assertions},
                    verdicts={v.verdict_id: v for v in verdicts}, states=states)
    if with_load_bearing:
        for vid in sorted(bundle.verdicts):
            bundle.load_bearing[vid] = load_bearing(bundle, vid, weights, model)
    return bundle


def leaf_claims(bundle: Bundle, verdict_id: str) -> list[str]:
    """Claims under a verdict's assertions (all sides and excluded claims), in order."""
    out: list[str] = []
    for aid in bundle.verdicts[verdict_id].assertion_ids:
        a = bundle.assertions[aid]
        for side in a.sides:
            out += side.claim_ids
        out += a.excluded_claim_ids
    return sorted(set(out))


def trace(bundle: Bundle, verdict_id: str) -> list[tuple[str, str, str]]:
    """Ordered (layer, id, label) nodes from the verdict down to its L0 sources."""
    v = bundle.verdicts[verdict_id]
    nodes = [("L5", v.verdict_id, f"{v.reaction_id} {v.element}: {v.status} {list(v.actions)} — {v.reason}")]
    for aid in v.assertion_ids:
        a = bundle.assertions[aid]
        nodes.append(("L4", aid, f"{a.subject_key} {a.predicate.value}: score {a.score}, n_independent "
                                 f"{a.n_independent}, {a.conflict_status}, resolved_by {a.resolved_by}, "
                                 f"trace {list(a.resolution_trace)}"))
        for pid in a.pair_ids:
            r = bundle.relations[pid]
            nodes.append(("L3", pid, f"{r.claim_a} ~ {r.claim_b}: {r.pair_class.value} {r.reason or ''}".strip()))
    for cid in leaf_claims(bundle, verdict_id):
        c = bundle.claims[cid]
        root = bundle.claims[c.derived_from] if c.derived_from else c
        nodes.append(("L2", cid, f"{c.predicate.value} {c.polarity} {c.claim_level.value} {c.assertion_strength} "
                                 f"taxon {c.taxon_distance}" + (f" derived from {c.derived_from}" if c.derived_from else "")))
        p = bundle.passages[root.passage_id]
        nodes.append(("L1", p.passage_id, f"\"{p.quote}\" [{p.char_start}:{p.char_end}]"))
        s = bundle.sources[root.source_id]
        nodes.append(("L0", s.source_id, f"{s.kind.value}: {s.title} ({s.year})"))
    return nodes


def load_bearing(bundle: Bundle, verdict_id: str, weights: Weights, model: ModelIdentity) -> tuple[str, ...]:
    """Leaves whose removal changes the verdict (status or actions) or an assertion score under it."""
    v = bundle.verdicts[verdict_id]
    state = bundle.states[v.reaction_id]
    bucket_claims = [c for c in bundle.claims.values() if c.subject.bucket == v.reaction_id]
    base_scores = {bundle.assertions[aid].subject_key: bundle.assertions[aid].score for aid in v.assertion_ids}
    changed = []
    for cid in leaf_claims(bundle, verdict_id):
        drop = {cid} | {c.claim_id for c in bucket_claims if c.derived_from == cid}
        kept = [c for c in bucket_claims if c.claim_id not in drop]
        _, assertions, verdicts = forward(kept, bundle.sources, {v.reaction_id: state}, weights, model)
        new = next(x for x in verdicts if x.element == v.element) if any(
            x.element == v.element for x in verdicts) else None
        if new is None:     # the element disappears without the claim (e.g. the only cofactor claim)
            changed.append(cid)
            continue
        scores = {a.subject_key: a.score for a in assertions if a.assertion_id in v.assertion_ids}
        if (new.status, new.actions) != (v.status, v.actions) or any(
                scores.get(k, 0) != s for k, s in base_scores.items()):
            changed.append(cid)
    return tuple(changed)


def chain_quality(assertion: Assertion) -> bool:
    """True when every chain bottoms out below an organism-specific functional leaf (homology_only)."""
    return assertion.homology_only
