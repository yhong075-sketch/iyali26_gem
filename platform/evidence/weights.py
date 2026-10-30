"""The committed rule tables (`evidence_weights.yaml`) as a typed object."""

from functools import cached_property
from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict, model_validator

from scripts.gem_annotate.model_layout import MODEL

from .schemas import ClaimLevel

WEIGHTS_FILE = "evidence_weights.yaml"
CURATION_DIR = MODEL.curation / "evidence_pipeline"


class LevelInfo(BaseModel):
    eco: str
    branch: str
    band: str


class ConfidenceRule(BaseModel):
    score: int
    levels: list[ClaimLevel] | str            # "any" matches every admissible level
    predicates: list[str] | None = None
    max_taxon: int = 1
    transfer_from: list[ClaimLevel] | None = None


class Weights(BaseModel):
    model_config = ConfigDict(extra="forbid")

    version: str
    eco: dict
    levels: dict[ClaimLevel, LevelInfo]
    directness: dict[str, list[list[ClaimLevel]]]
    confidence: list[ConfidenceRule]
    unresolved_cap: int
    taxonomy: dict
    compartments: dict[str, str]
    subsystem_precedence: list[str]
    essentiality_threshold: float
    essentiality_sensitivity: list[float]
    priority: dict
    independence_window_years: int
    independence_min_shared_authors: int

    @model_validator(mode="after")
    def _complete(self) -> "Weights":
        missing = set(ClaimLevel) - set(self.levels)
        if missing:
            raise ValueError(f"levels without ECO/band entry: {sorted(m.value for m in missing)}")
        if self.essentiality_threshold not in self.essentiality_sensitivity:
            raise ValueError("essentiality_threshold must be one of essentiality_sensitivity")
        return self

    @cached_property
    def _ranks(self) -> dict[str, dict[ClaimLevel, int]]:
        return {key: {level: index for index, tier in enumerate(tiers) for level in tier}
                for key, tiers in self.directness.items()}

    def admissible(self, predicate_key: str, level: ClaimLevel) -> bool:
        return level in self._ranks[predicate_key]

    def rank(self, predicate_key: str, level: ClaimLevel) -> int:
        """0 is most direct. Raises KeyError for an inadmissible level."""
        return self._ranks[predicate_key][level]

    def band(self, level: ClaimLevel) -> str:
        return self.levels[level].band

    def eco_for(self, level: ClaimLevel) -> str:
        return self.levels[level].eco

    def claim_score(self, predicate_key: str, level: ClaimLevel, taxon_distance: int,
                    transfer_from: ClaimLevel | None) -> int:
        """Score of one admissible claim: the highest matching confidence rule."""
        for rule in sorted(self.confidence, key=lambda r: -r.score):
            if rule.levels != "any" and level not in rule.levels:
                continue
            if rule.predicates is not None and predicate_key not in rule.predicates:
                continue
            if rule.levels == "any":
                return rule.score
            if taxon_distance > rule.max_taxon:
                continue
            if rule.transfer_from is not None and transfer_from not in rule.transfer_from:
                continue
            return rule.score
        return 0


def load_weights(path: Path | None = None) -> Weights:
    path = path or MODEL.curation_file(WEIGHTS_FILE)
    with Path(path).open() as handle:
        return Weights.model_validate(yaml.safe_load(handle))
