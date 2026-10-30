"""Typed records of the evidence chain (spec §1.2–1.4, §1.7, §2, §6).

L0 Source → L1 Passage → L2 Claim → L3 PairRelation → L4 Assertion → L5 Verdict → L6 Decision.
The reviewer contracts (R1–R3) are defined here as types only; nothing calls a reviewer in M1.
"""

from enum import Enum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, field_serializer, field_validator, model_validator

from gem_review.lit.models import short_quote


class Predicate(str, Enum):
    REACTION_EXISTS = "reaction_exists"
    FLUX_ACTIVE_UNDER = "flux_active_under"
    GENE_ENCODES = "gene_encodes"
    LOCALIZED_IN = "localized_in"
    DIRECTIONALITY = "directionality"
    COFACTOR = "cofactor"
    SUBSYSTEM_MEMBER = "subsystem_member"
    ESSENTIAL_UNDER = "essential_under"


class ClaimLevel(str, Enum):
    STRUCTURE = "STRUCTURE"
    ACTIVITY_PURIFIED = "ACTIVITY_PURIFIED"
    ACTIVITY_HETEROLOGOUS = "ACTIVITY_HETEROLOGOUS"
    FLUX_13C = "FLUX_13C"
    ACTIVITY_CRUDE = "ACTIVITY_CRUDE"
    LOCALIZATION_DIRECT = "LOCALIZATION_DIRECT"
    PHENOTYPE_COMPLEMENTATION = "PHENOTYPE_COMPLEMENTATION"
    PHENOTYPE_KO = "PHENOTYPE_KO"
    PHENOTYPE_OE = "PHENOTYPE_OE"
    EXPRESSION_PROTEIN = "EXPRESSION_PROTEIN"
    EXPRESSION_TRANSCRIPT = "EXPRESSION_TRANSCRIPT"
    PHYSIOLOGY = "PHYSIOLOGY"
    ORTHOLOGY = "ORTHOLOGY"
    HOMOLOGY_ANNOTATION = "HOMOLOGY_ANNOTATION"
    GENOMIC_PRESENCE = "GENOMIC_PRESENCE"
    AUTHOR_STATEMENT_TRACEABLE = "AUTHOR_STATEMENT_TRACEABLE"
    AUTHOR_STATEMENT_NONTRACEABLE = "AUTHOR_STATEMENT_NONTRACEABLE"
    MODEL_INCLUSION = "MODEL_INCLUSION"
    MODEL_VALIDATION = "MODEL_VALIDATION"


ACTIVITY_LEVELS = frozenset({ClaimLevel.ACTIVITY_PURIFIED, ClaimLevel.ACTIVITY_HETEROLOGOUS, ClaimLevel.ACTIVITY_CRUDE})


class SourceKind(str, Enum):
    ARTICLE = "article"
    PREPRINT = "preprint"
    THESIS = "thesis"
    PATENT = "patent"
    DATASET = "dataset"
    DATABASE_ENTRY = "database_entry"
    PRIOR_MODEL = "prior_model"
    LEGACY_CASE = "legacy_case"     # D7: the 240 essentiality FN cases, never live evidence


class SourceTier(str, Enum):
    A = "A"
    B = "B"
    C = "C"
    D = "D"
    E = "E"


class PairClass(str, Enum):
    AGREE_INDEPENDENT = "agree_independent"
    AGREE_DERIVATIVE = "agree_derivative"
    SCOPED_APART = "scoped_apart"
    SUPERSEDED = "superseded"
    ANOMALY_CROSS_LEVEL = "anomaly_cross_level"
    CONTRADICTION_STRONG = "contradiction_strong"
    CONTRADICTION_WEAK = "contradiction_weak"


CONFLICT_CLASSES = frozenset({PairClass.CONTRADICTION_STRONG, PairClass.CONTRADICTION_WEAK,
                              PairClass.SUPERSEDED, PairClass.ANOMALY_CROSS_LEVEL})


class Frozen(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")


# ── L0 ───────────────────────────────────────────────────────────────────────────────────────

class Source(Frozen):
    source_id: str
    kind: SourceKind
    title: str
    year: int | None = None
    doi: str | None = None
    pmid: str | None = None
    tier: SourceTier | None = None
    peer_reviewed: bool
    unpublished: bool = False                     # D9
    legacy_scoring: bool = False                  # D7
    retraction_status: Literal["not_retracted", "retracted", "suspected", "unknown"] = "unknown"
    organisms: tuple[str, ...] = ()
    authors: tuple[str, ...] = ()
    corresponding_author: str | None = None
    parent_model: str | None = None               # D8: the next model up the lineage chain
    text_ref: str | None = None                   # where the verified text lives (path or cache key)
    text_sha256: str | None = None
    subjects: tuple[str, ...] = ()                # reactions/genes a source is indexed under
    properties: dict[str, str] = Field(default_factory=dict)  # e.g. dataset cutoffs, legacy verdicts


# ── L1 ───────────────────────────────────────────────────────────────────────────────────────

class Passage(Frozen):
    passage_id: str
    source_id: str
    quote: str                                    # ≤ 25 words, whitespace-normalised (D3)
    section: str | None = None
    char_start: int                               # offsets into the source text named by text_ref
    char_end: int
    quote_sha256: str
    verified: bool

    @field_validator("quote")
    @classmethod
    def _short(cls, value: str) -> str:
        return short_quote(value)


# ── L2 ───────────────────────────────────────────────────────────────────────────────────────

class Condition(Frozen):
    carbon_source: str | None = None
    nitrogen: str | None = None
    oxygen: str | None = None
    growth_phase: str | None = None
    adaptation: str | None = None
    medium: str | None = None                     # D4
    marker_complementation: str | None = None     # D4
    free_text: str | None = None

    def overlaps(self, other: "Condition | None") -> bool:
        """Every field specified in both is equal (free_text is description only)."""
        if other is None:
            return True
        for name in type(self).model_fields:
            if name == "free_text":
                continue
            mine, theirs = getattr(self, name), getattr(other, name)
            if mine is not None and theirs is not None and mine.casefold() != theirs.casefold():
                return False
        return True


class SubjectRef(Frozen):
    """Normalised subject. GENE_ENCODES subjects are (gene, reaction); ESSENTIAL_UNDER may be a gene."""

    kind: Literal["reaction", "gene_reaction", "gene"]
    reaction_id: str | None = None
    gene_id: str | None = None
    status: Literal["resolved", "unresolved", "unmapped"] = "resolved"
    raw: str | None = None                        # the identifier as written in the source
    candidate_subjects: tuple[str, ...] = ()

    @property
    def key(self) -> str:
        """Assertion key (one L4 node per key and predicate)."""
        if self.kind == "gene_reaction":
            return f"{self.gene_id}|{self.reaction_id}"
        return self.reaction_id if self.kind == "reaction" else self.gene_id

    @property
    def bucket(self) -> str:
        """Pairing bucket: GENE_ENCODES claims about different genes of one reaction are compared."""
        return self.reaction_id if self.kind in ("reaction", "gene_reaction") else self.gene_id


FLAG_PREDICATES = (Predicate.REACTION_EXISTS, Predicate.FLUX_ACTIVE_UNDER, Predicate.ESSENTIAL_UNDER)


class FlagValue(Frozen):
    """Boolean predicates: the value is carried by the claim's polarity."""

    predicate: Literal[Predicate.REACTION_EXISTS, Predicate.FLUX_ACTIVE_UNDER, Predicate.ESSENTIAL_UNDER]


class GeneEncodesValue(Frozen):
    predicate: Literal[Predicate.GENE_ENCODES]
    gene_id: str
    reaction_id: str
    exclusive: bool = False


class LocalizedValue(Frozen):
    predicate: Literal[Predicate.LOCALIZED_IN]
    compartments: frozenset[str]                  # model compartment IDs

    @field_serializer("compartments")
    def _sorted(self, value: frozenset[str]) -> list[str]:   # set order depends on the hash seed
        return sorted(value)


class DirectionalityValue(Frozen):
    predicate: Literal[Predicate.DIRECTIONALITY]
    direction: Literal["forward", "reverse", "reversible"]   # relative to the model's equation
    level: Literal["biochemical", "physiological"]


class CofactorValue(Frozen):
    predicate: Literal[Predicate.COFACTOR]
    roles: dict[str, frozenset[str]]              # role -> model species base IDs (e.g. m41)

    @field_serializer("roles")
    def _sorted(self, value: dict[str, frozenset[str]]) -> dict[str, list[str]]:
        return {role: sorted(ids) for role, ids in sorted(value.items())}


class SubsystemValue(Frozen):
    predicate: Literal[Predicate.SUBSYSTEM_MEMBER]
    pathways: frozenset[str]
    ontology: str

    @field_serializer("pathways")
    def _sorted(self, value: frozenset[str]) -> list[str]:
        return sorted(value)


ClaimValue = Annotated[FlagValue | GeneEncodesValue | LocalizedValue | DirectionalityValue | CofactorValue
                       | SubsystemValue, Field(discriminator="predicate")]


class Claim(Frozen):
    claim_id: str
    source_id: str
    passage_id: str
    subject: SubjectRef
    predicate: Predicate
    value: ClaimValue
    polarity: Literal["positive", "negative"]
    assertion_strength: Literal["assertive", "hedged", "speculative"]
    claim_level: ClaimLevel
    level_justification: str | None = None        # stored, never used downstream
    transfer_from_level: ClaimLevel | None = None # ORTHOLOGY only
    organism: str
    strain: str | None
    taxon_distance: int = Field(ge=0, le=5)       # computed in code (taxonomy.py)
    compartment: str | None = None
    condition: Condition | None = None
    strain_specific: bool = False
    method: str
    is_original_finding: bool
    cited_support: tuple[str, ...] = ()
    rebuts: tuple[str, ...] = ()
    eco: str
    extraction_run_id: str
    authors: tuple[str, ...] = ()
    corresponding_author: str | None = None
    year: int
    derived_from: str | None = None               # FLUX_ACTIVE_UNDER → REACTION_EXISTS derivation

    @model_validator(mode="after")
    def _consistent(self) -> "Claim":
        if self.value.predicate != self.predicate:
            raise ValueError(f"value is for {self.value.predicate.value}, claim is {self.predicate.value}")
        if self.transfer_from_level is not None and self.claim_level != ClaimLevel.ORTHOLOGY:
            raise ValueError("transfer_from_level is only for ORTHOLOGY claims")
        return self

    @property
    def predicate_key(self) -> str:
        """Key into the directness/admissibility tables (directionality splits by level)."""
        if self.predicate == Predicate.DIRECTIONALITY:
            return f"directionality.{self.value.level}"
        return self.predicate.value


class RejectedClaim(Frozen):
    reject_id: str
    source_id: str | None
    reason: str
    detail: str
    payload: dict = Field(default_factory=dict)


# ── L3 ───────────────────────────────────────────────────────────────────────────────────────

class PairRelation(Frozen):
    pair_id: str
    claim_a: str
    claim_b: str
    pair_class: PairClass
    reason: str | None = None
    independence_shared: bool = False


# ── L4 ───────────────────────────────────────────────────────────────────────────────────────

class Side(Frozen):
    polarity: Literal["positive", "negative"]
    value_key: str
    claim_ids: tuple[str, ...]
    best_level: ClaimLevel
    best_rank: int
    score: int = Field(ge=0, le=4)
    n_independent: int
    min_taxon_distance: int
    latest_year: int


class Assertion(Frozen):
    assertion_id: str
    subject_key: str
    predicate: Predicate
    sides: tuple[Side, ...]
    excluded_claim_ids: tuple[str, ...] = ()      # inadmissible or speculative, with nothing dropped
    pair_ids: tuple[str, ...] = ()
    conflict_status: Literal["NONE", "CONTRADICTION_STRONG", "CONTRADICTION_WEAK", "ANOMALY"]
    resolution_trace: tuple[tuple[str, str], ...]
    resolved: bool
    winner: int | None                            # side that won a contradiction (None if there was none)
    standing: tuple[int, ...]                     # sides that hold: the winner plus sides outside the contradiction
    resolved_by: str | None = None
    score: int = Field(ge=0, le=4)
    flag: Literal["CONFLICTING"] | None = None
    n_independent: int = 0
    max_taxon_distance_used: int | None = None
    homology_only: bool = False
    inherited_from: str | None = None             # D8
    advisory_flags: tuple[str, ...] = ()          # R7, never decides
    weights_version: str

    @property
    def standing_sides(self) -> list[Side]:
        return [self.sides[i] for i in self.standing]


# ── L5 / L6 ──────────────────────────────────────────────────────────────────────────────────

VerdictStatus = Literal["SUPPORTED", "CONTRADICTED", "UNSUPPORTED", "CONFLICTING", "PARTIAL"]
Element = Literal["existence", "gpr", "compartment", "bounds", "cofactor", "subsystem"]


class ModelElementState(Frozen):
    """What the model says about one reaction; extracted once so verdicts need no cobra."""

    reaction_id: str
    name: str
    gene_ids: tuple[str, ...]                     # canonical IDs
    gpr: str
    compartments: tuple[str, ...]
    lower_bound: float
    upper_bound: float
    groups: tuple[str, ...]
    species_base_ids: tuple[str, ...]
    essential_genes: tuple[str, ...] = ()         # genes essential in simulation (threshold from YAML)
    impact_known: bool = True                     # False when no knockout table exists for this model
    hub: bool = False


class Verdict(Frozen):
    verdict_id: str
    reaction_id: str
    element: Element
    status: VerdictStatus
    actions: tuple[str, ...] = ()
    reason: str
    assertion_ids: tuple[str, ...] = ()
    priority: float
    priority_note: str | None = None
    model_path: str
    model_commit: str
    model_sha256: str
    weights_version: str


class Decision(Frozen):
    decision_id: str
    verdict_id: str
    action: str
    rule_trace: tuple[tuple[str, str], ...]
    reviewer_verdict_id: str | None = None
    curator: str
    date: str
    rationale: str
    commit: str | None = None
    patch_id: str | None = None


# ── Reviewer contracts (spec §6), types only in M1 ───────────────────────────────────────────

class ConflictPacket(Frozen):
    subject: SubjectRef
    predicate: Predicate
    claims: tuple[Claim, ...]
    passages: tuple[Passage, ...]
    pair_relations: tuple[PairRelation, ...]
    resolution_trace: tuple[tuple[str, str], ...]
    proposed_verdict: Verdict | None
    model_state: ModelElementState
    network_flags: tuple[str, ...]


class ReviewVerdict(Frozen):
    agrees_with_rule_engine: bool
    verdict: Literal["side_positive", "side_negative", "unresolved", "scope_error", "extraction_error"]
    rationale: str
    disputed_rule: str | None
    disputed_claims: tuple[str, ...]
    requested_evidence: tuple[str, ...]
    confidence: Literal["low", "medium", "high"]


class ClaimChallenge(Frozen):
    claim_id: str
    challenged_fields: tuple[Literal["polarity", "assertion_strength", "claim_level", "scope"], ...]
    upheld: bool
    rationale: str
    proposed_values: dict[str, str] = Field(default_factory=dict)   # triggers re-extraction, never an edit


class DecisionReview(Frozen):
    verdict_id: str
    action: str
    case_for: str
    case_against: str
    recommendation: Literal["apply", "do_not_apply", "needs_more_evidence"]
    confidence: Literal["low", "medium", "high"]
