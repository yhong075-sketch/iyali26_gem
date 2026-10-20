"""Search terms and Europe PMC queries for one gene–reaction pair."""

import re
from dataclasses import dataclass

ORGANISM = '"Yarrowia lipolytica"'
# Case-insensitive: the locus crosswalk returns YALI0 tags in lower case.
_LOCUS = re.compile(r"^(YALI[01])_?([A-F]\d{5})([a-z])$", re.IGNORECASE)


@dataclass(frozen=True, slots=True)
class PairTerms:
    gene_terms: tuple[str, ...]      # locus spellings, YALI0 counterpart, gene symbols
    reaction_terms: tuple[str, ...]  # EC codes and the cleaned reaction name

    def gene_query(self) -> str:
        return f"({_or(self.gene_terms)}) AND {ORGANISM}"

    def pair_query(self) -> str | None:
        """None when the reaction has no usable term, so the pair query would equal the gene query."""
        if not self.reaction_terms:
            return None
        return f"{self.gene_query()} AND ({_or(self.reaction_terms)})"


def _or(terms) -> str:
    return " OR ".join(f'"{term}"' for term in terms)


def locus_spellings(locus: str) -> list[str]:
    match = _LOCUS.match(locus)
    if not match:
        return [locus]
    assembly, body, suffix = match.groups()
    canonical = f"{assembly.upper()}{body.upper()}{suffix.lower()}"
    return [canonical, f"{assembly.upper()}_{body.upper()}{suffix.lower()}"]


def clean_reaction_name(name: str) -> str | None:
    """Drop parenthetical qualifiers such as '(octadecanoyl-CoA)'; skip names too short to search."""
    cleaned = re.sub(r"\s*\([^)]*\)", "", name or "").strip(" ,;")
    return cleaned if len(cleaned) >= 4 else None


def build_terms(gene: str, yali0: str | None, symbols, reaction_name: str, ec_codes) -> PairTerms:
    gene_terms = locus_spellings(gene)
    if yali0:
        gene_terms += locus_spellings(yali0)
    for symbol in symbols:
        if symbol not in gene_terms and not _LOCUS.match(symbol):
            gene_terms.append(symbol)
    reaction_terms = list(ec_codes)
    name = clean_reaction_name(reaction_name)
    if name:
        reaction_terms.append(name)
    return PairTerms(tuple(dict.fromkeys(gene_terms)), tuple(dict.fromkeys(reaction_terms)))
