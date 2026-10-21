"""Validated extraction output for one gene–reaction pair.

The agent returns only `ClaimSet`; `PairReview` is assembled in Python so the evidence level is
computed from validated claims and never set by the model.
"""

import re
from enum import IntEnum

from pydantic import BaseModel, Field, field_validator

_WS = re.compile(r"\s+")
MAX_QUOTE_WORDS = 25


def short_quote(value: str) -> str:
    """Whitespace-normalised quote of 1 to MAX_QUOTE_WORDS words."""
    value = _WS.sub(" ", value).strip()
    if not value:
        raise ValueError("quote is empty")
    if len(value.split()) > MAX_QUOTE_WORDS:
        raise ValueError(f"quote has more than {MAX_QUOTE_WORDS} words")
    return value


class EvidenceType(IntEnum):
    """Checklist sheet "Evidence Levels". IntEnum so levels compare and cap numerically."""

    biochemical = 4
    genetic = 3
    physiological = 2
    homology = 1
    none = 0


class Claim(BaseModel):
    pmid: str = Field(description="PMID of the cited article; must be one of the retrieved PMIDs")
    supports_gene_reaction: bool = Field(
        description="True if the quote supports that this gene catalyses this reaction; False if it contradicts it")
    evidence_type: EvidenceType = Field(
        description="4 biochemical (purified enzyme / in vitro activity), 3 genetic (knockout or "
                    "overexpression phenotype in the organism), 2 physiological (growth or secretion "
                    "phenotype), 1 homology (ortholog of a characterised enzyme), 0 none")
    compartment: str | None = Field(default=None, description="Subcellular compartment exactly as stated in the paper, else null")
    organism_is_yarrowia: bool = Field(description="False if the evidence comes from another organism")
    quote: str = Field(description=f"At most {MAX_QUOTE_WORDS} words copied verbatim from the cited chunk")
    chunk_id: str = Field(description="chunk_id of the chunk the quote is copied from")

    @field_validator("quote")
    @classmethod
    def _short_quote(cls, value: str) -> str:
        return short_quote(value)

    @property
    def effective_level(self) -> EvidenceType:
        """Evidence from another organism counts at most as homology."""
        if not self.organism_is_yarrowia:
            return min(self.evidence_type, EvidenceType.homology)
        return self.evidence_type


class ClaimSet(BaseModel):
    """What the agent returns. Empty when the chunks do not establish the link."""

    claims: list[Claim] = Field(default_factory=list)


class PairReview(BaseModel):
    gene: str
    reaction_id: str
    claims: list[Claim]
    best_evidence_level: EvidenceType
    no_organism_specific_evidence: bool

    @classmethod
    def from_claims(cls, gene: str, reaction_id: str, claims: list[Claim]) -> "PairReview":
        supporting = [c for c in claims if c.supports_gene_reaction]
        best = max((c.effective_level for c in supporting), default=EvidenceType.none)
        return cls(gene=gene, reaction_id=reaction_id, claims=claims, best_evidence_level=best,
                   no_organism_specific_evidence=not any(c.organism_is_yarrowia for c in supporting))

    def claims_disagree(self) -> bool:
        """Mixed support/contradiction, or more than one stated compartment."""
        verdicts = {c.supports_gene_reaction for c in self.claims}
        return len(verdicts) > 1 or compartments_disagree(c.compartment for c in self.claims)


def normalise_text(text: str) -> str:
    return _WS.sub(" ", text).strip()


# Spelling variants papers use for one compartment; compound statements ("ER and peroxisome")
# are left as written because they are a different claim.
_COMPARTMENT_SYNONYMS = {
    "peroxisomal": "peroxisome", "peroxisomes": "peroxisome",
    "mitochondrial": "mitochondria", "mitochondrion": "mitochondria",
    "er": "endoplasmic reticulum", "cytosolic": "cytosol", "cytoplasm": "cytosol",
    "cytoplasmic": "cytosol", "lipid bodies": "lipid body", "lipid droplets": "lipid body",
    "lipid droplet": "lipid body",
}


def canonical_compartment(value: str | None) -> str | None:
    if not value:
        return None
    key = normalise_text(value).lower().rstrip(".")
    return _COMPARTMENT_SYNONYMS.get(key, key)


def compartments_disagree(values) -> bool:
    return len({canonical_compartment(v) for v in values if v}) > 1
