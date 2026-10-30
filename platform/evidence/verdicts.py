"""L5 verdict engine (spec §2.2): the model's current state against the L4 assertions.

`reaction_verdicts` is pure: it takes a `ModelElementState` (extracted once from the model by
`element_states`) so verdicts can be recomputed for any model commit and for leave-one-out.
Verdicts read every standing side of an assertion (sides scoped apart all hold). A verdict that
proposes changing the model (CONTRADICTED) needs a standing side with score ≥ 3, i.e.
organism-specific functional evidence; weaker negative evidence is flagged, never acted on.
UNSUPPORTED is never a removal proposal: gap-filled and gene-less reactions may be needed.
"""

import hashlib
import json
from pathlib import Path

from .schemas import Assertion, ClaimLevel, ModelElementState, Predicate, Side, Verdict
from .weights import Weights

GAP_NOTE = "gap-filled or gene-less reactions may be needed; this is never a removal proposal"
STRONG = 3


def verdict_id(reaction_id: str, element: str) -> str:
    return "V-" + hashlib.sha256(f"{reaction_id}\0{element}".encode()).hexdigest()[:12]


def model_direction(state: ModelElementState) -> str:
    lb, ub = state.lower_bound, state.upper_bound
    if lb == 0 and ub == 0:
        return "closed"
    if lb < 0 < ub:
        return "reversible"
    return "forward" if lb >= 0 else "reverse"


def _priority(status: str, state: ModelElementState, weights: Weights) -> tuple[float, str | None]:
    cfg = weights.priority["impact"]
    impact = cfg["base"] + (cfg["phosphatidate_hub"] if state.hub else 0)
    if state.impact_known:
        impact += cfg["essential_in_simulation"] if state.essential_genes else 0
        note = None
    else:
        note = "impact unknown: no knockout table for this model; essentiality term left out"
    return float(impact * weights.priority["uncertainty"][status]), note


def _strong(sides: list[Side]) -> list[Side]:
    return [s for s in sides if s.score >= STRONG]


def _weak_status(a: Assertion, extra: str = "") -> tuple[str, tuple[str, ...], str]:
    """No strong standing side: PARTIAL at score 2, otherwise UNSUPPORTED."""
    positive = max((s.score for s in a.standing_sides if s.polarity == "positive"), default=0)
    negative = max((s.score for s in a.standing_sides if s.polarity == "negative"), default=0)
    note = f"; weak negative evidence (score {negative}) is never acted on" if negative else ""
    if positive == 2:
        return "PARTIAL", ("flag_for_evidence",), f"score 2: below the SUPPORTED threshold of 3{note}{extra}"
    return "UNSUPPORTED", ("flag_for_evidence",), f"score {positive}{note}; {GAP_NOTE}{extra}"


def _unresolved(a: Assertion | None, what: str):
    if a is not None and not a.resolved:
        return "CONFLICTING", ("human_review",), f"{what} claims conflict and R2–R6 leave a tie", (a.assertion_id,)
    return None


def _existence(state, get):
    a = get(state.reaction_id, Predicate.REACTION_EXISTS)
    if a is None or not a.sides:
        return "UNSUPPORTED", ("flag_for_evidence",), f"no admissible claims; {GAP_NOTE}", ()
    if (c := _unresolved(a, "existence")):
        return c
    pos = _strong([s for s in a.standing_sides if s.polarity == "positive"])
    neg = _strong([s for s in a.standing_sides if s.polarity == "negative"])
    if pos:
        scoped = "; a strong negative is scoped apart (strain or level)" if neg else ""
        return "SUPPORTED", (), f"score {max(s.score for s in pos)}{scoped}", (a.assertion_id,)
    if neg:
        how = f"contradiction resolved by {a.resolved_by}" if a.resolved_by else "no contradiction with the negative"
        return "CONTRADICTED", ("propose_remove",), f"negative evidence score {max(s.score for s in neg)}; {how}", \
            (a.assertion_id,)
    return (*_weak_status(a), (a.assertion_id,))


def _gpr(state, gene_assertions):
    model_genes = set(state.gene_ids)
    used = tuple(a.assertion_id for a in gene_assertions)
    if any(not a.resolved for a in gene_assertions):
        return "CONFLICTING", ("human_review",), "a gene-reaction assertion is unresolved", used
    if any(rule == "X" for a in gene_assertions for rule, _ in a.resolution_trace) and len(model_genes) > 1:
        return "CONFLICTING", ("human_review",), "an exclusive gene claim conflicts with another gene's claim", used
    add, remove, supported, flagged = [], [], set(), []
    for a in gene_assertions:
        gene = a.subject_key.split("|")[0]
        standing = a.standing_sides
        isozymes = " or " in f" {state.gpr.lower()} " and model_genes - {gene}
        if isozymes and any(s.polarity == "positive" and s.value_key == "exclusive" and s.score >= STRONG
                            for s in standing):
            return "CONFLICTING", ("human_review",), f"exclusive claim (score ≥ 3) for {gene} vs model isozymes " \
                f"{sorted(model_genes - {gene})}", used
        pos = _strong([s for s in standing if s.polarity == "positive"])
        neg = _strong([s for s in standing if s.polarity == "negative"])
        if pos:
            supported.add(gene)
            if gene not in model_genes:
                add.append(gene)
        elif neg and gene in model_genes:
            # A negative knockout alone is explained by isozymes (spec §1.6): flagged, never a removal.
            if all(s.best_level == ClaimLevel.PHENOTYPE_KO for s in neg):
                flagged.append(gene)
            else:
                remove.append(gene)
    if add or remove:
        actions = tuple([f"add_gene:{g}" for g in sorted(add)] + [f"remove_gene:{g}" for g in sorted(remove)])
        return "PARTIAL", actions, "supported gene missing or model gene contradicted", used
    if model_genes and model_genes <= supported:
        return "SUPPORTED", (), "every model gene has a supported (score ≥ 3) assertion", used
    if not model_genes:
        return "UNSUPPORTED", ("flag_for_evidence",), f"gene-less reaction; {GAP_NOTE}", used
    missing = sorted(model_genes - supported)
    note = f"; negative knockout only (isozymes possible) for {', '.join(flagged)}" if flagged else ""
    return "UNSUPPORTED", ("flag_for_evidence",), f"no score ≥ 3 support for {', '.join(missing)}{note}", used


def _compartment(state, get):
    a = get(state.reaction_id, Predicate.LOCALIZED_IN)
    if a is None or not a.sides:
        return "UNSUPPORTED", ("flag_for_evidence",), "no localization claims", ()
    if (c := _unresolved(a, "localization")):
        return c
    model = set(state.compartments)
    places = lambda side: set(side.value_key.split(",")) if side.value_key else set()  # noqa: E731
    pos = _strong([s for s in a.standing_sides if s.polarity == "positive"])
    neg = _strong([s for s in a.standing_sides if s.polarity == "negative"])
    claimed = set().union(*(places(s) for s in pos)) if pos else set()
    excluded = set().union(*(places(s) for s in neg)) if neg else set()
    if claimed & model:
        return "SUPPORTED", (), f"evidence places it in {sorted(claimed)}", (a.assertion_id,)
    if claimed or excluded & model:
        why = f"evidence places it in {sorted(claimed)}" if claimed else f"evidence excludes {sorted(excluded & model)}"
        return "CONTRADICTED", ("relocate",), f"{why}; model in {sorted(model)}", (a.assertion_id,)
    return (*_weak_status(a), (a.assertion_id,))


def _incompatible_direction(model_dir: str, polarity: str, direction: str, level: str) -> bool:
    """Bounds are a capability. In vivo flux one way does not contradict a reversible reaction, in vitro
    reversibility does not contradict a physiologically irreversible one, and flux not observed in vivo
    does not remove a capability. Evidence that a closed reaction runs does contradict it."""
    if model_dir == "closed":
        return polarity == "positive"
    if polarity == "negative":          # "does not run <direction>": only an enzyme property removes a capability
        if level != "biochemical":
            return False
        return direction == model_dir or model_dir == "reversible"
    if {model_dir, direction} == {"forward", "reverse"}:
        return True
    if direction in ("forward", "reverse") and model_dir == "reversible":
        return level == "biochemical"
    return direction == "reversible" and model_dir in ("forward", "reverse") and level == "physiological"


def _bounds(state, get):
    a = get(state.reaction_id, Predicate.DIRECTIONALITY)
    if a is None or not a.sides:
        return "UNSUPPORTED", ("flag_for_evidence",), "no directionality claims", ()
    if (c := _unresolved(a, "directionality")):
        return c
    model_dir = model_direction(state)
    strong = _strong(a.standing_sides)
    clashing = [s for s in strong if _incompatible_direction(model_dir, s.polarity, *s.value_key.split("@"))]
    if clashing:
        claims = ", ".join(f"{s.polarity} {s.value_key}" for s in clashing)
        return "CONTRADICTED", ("change_bounds",), f"evidence: {claims}; model: {model_dir}", (a.assertion_id,)
    if strong:
        return "SUPPORTED", (), f"directionality evidence consistent with model ({model_dir})", (a.assertion_id,)
    return (*_weak_status(a), (a.assertion_id,))


def _cofactor(state, get):
    a = get(state.reaction_id, Predicate.COFACTOR)
    if a is None or not a.sides:
        return None
    if (c := _unresolved(a, "cofactor")):
        return c
    present = set(state.species_base_ids)
    wrong = []
    for side in _strong(a.standing_sides):
        for role, ids in json.loads(side.value_key or "{}").items():
            hit = bool(set(ids) & present)
            if (side.polarity == "positive" and not hit) or (side.polarity == "negative" and hit):
                wrong.append(f"{side.polarity} {role}")
    if wrong:
        return "CONTRADICTED", ("change_stoichiometry",), f"cofactor evidence disagrees: {wrong}", (a.assertion_id,)
    if _strong(a.standing_sides):
        return "SUPPORTED", (), "cofactor evidence consistent with the model", (a.assertion_id,)
    return (*_weak_status(a), (a.assertion_id,))


def _subsystem(state, get, weights, subsystem_claims):
    """Union of claimed pathways; primary by the committed precedence list; never CONTRADICTED."""
    a = get(state.reaction_id, Predicate.SUBSYSTEM_MEMBER)
    claimed = subsystem_claims.get(state.reaction_id, [])
    if a is None or not a.sides or not claimed:
        return "UNSUPPORTED", ("flag_for_evidence",), "no subsystem claims", ()
    basis = f" (only parent-model inclusion: {a.inherited_from})" if a.inherited_from else ""
    for ontology in weights.subsystem_precedence:
        pathways = sorted(p for o, p in claimed if o == ontology)
        if not pathways:
            continue
        primary = pathways[0]
        if ontology == "model_groups" and primary not in state.groups:
            return "PARTIAL", (f"reassign:{primary}",), \
                f"primary subsystem by precedence is {primary!r}; model groups {list(state.groups)}{basis}", \
                (a.assertion_id,)
        return "SUPPORTED", (), f"primary subsystem {primary!r} ({ontology}) agrees with the model{basis}", \
            (a.assertion_id,)
    return "UNSUPPORTED", ("flag_for_evidence",), "no claim in an ontology on the precedence list", (a.assertion_id,)


def reaction_verdicts(state: ModelElementState, assertions: dict[tuple[str, Predicate], Assertion],
                      weights: Weights, *, model_path: str, model_commit: str, model_sha256: str,
                      subsystem_claims: dict[str, list[tuple[str, str]]] | None = None) -> list[Verdict]:
    """All element verdicts for one reaction. `subsystem_claims` maps reaction → (ontology, pathway)."""

    def get(key: str, predicate: Predicate) -> Assertion | None:
        return assertions.get((key, predicate))

    gene_assertions = sorted((a for (key, p), a in assertions.items()
                              if p == Predicate.GENE_ENCODES and key.endswith(f"|{state.reaction_id}")),
                             key=lambda a: a.subject_key)
    results = {
        "existence": _existence(state, get),
        "gpr": _gpr(state, gene_assertions),
        "compartment": _compartment(state, get),
        "bounds": _bounds(state, get),
        "subsystem": _subsystem(state, get, weights, subsystem_claims or {}),
    }
    cofactor = _cofactor(state, get)
    if cofactor is not None:
        results["cofactor"] = cofactor
    verdicts = []
    for element in sorted(results):
        status, actions, reason, used = results[element]
        priority, note = _priority(status, state, weights)
        verdicts.append(Verdict(
            verdict_id=verdict_id(state.reaction_id, element), reaction_id=state.reaction_id, element=element,
            status=status, actions=actions, reason=reason, assertion_ids=tuple(used), priority=priority,
            priority_note=note, model_path=model_path, model_commit=model_commit, model_sha256=model_sha256,
            weights_version=weights.version))
    return verdicts


def element_states(model, reaction_ids: list[str], canonical_gene, essential_genes: set[str] | None,
                   hub_reactions: set[str]) -> dict[str, ModelElementState]:
    """Extract `ModelElementState` for the given reactions; `essential_genes=None` means unknown."""
    states = {}
    for rid in reaction_ids:
        r = model.reactions.get_by_id(rid)
        genes = tuple(sorted(canonical_gene(g.id) for g in r.genes))
        states[rid] = ModelElementState(
            reaction_id=rid, name=r.name, gene_ids=genes, gpr=r.gene_reaction_rule,
            compartments=tuple(sorted(r.compartments)), lower_bound=r.lower_bound, upper_bound=r.upper_bound,
            groups=tuple(sorted(g.name for g in model.groups if r in g.members)),
            species_base_ids=tuple(sorted({m.id.split("[")[0] for m in r.metabolites})),
            essential_genes=tuple(sorted(g for g in genes if g in (essential_genes or set()))),
            impact_known=essential_genes is not None, hub=rid in hub_reactions)
    return states


def model_identity(model_path: Path, repo_root: Path) -> tuple[str, str, str]:
    """(path relative to the repo, last commit touching it, SHA-256 of the file)."""
    import subprocess

    rel = model_path.resolve().relative_to(repo_root.resolve()).as_posix()
    commit = subprocess.run(["git", "-C", str(repo_root), "log", "-1", "--format=%H", "--", rel],
                            capture_output=True, text=True, check=True).stdout.strip() or "uncommitted"
    sha = hashlib.sha256(model_path.read_bytes()).hexdigest()
    return rel, commit, sha
