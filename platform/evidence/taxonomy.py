"""taxon_distance relative to Y. lipolytica W29, computed in code (spec §1.1, §8; D5)."""

from .weights import Weights


def _norm(value: str) -> str:
    return "".join(value.split()).casefold()


def taxon_distance(organism: str, strain: str | None, weights: Weights) -> int:
    """0 canonical strain set, 1 other/unknown Y. lipolytica strain, 2 other Yarrowia,
    3 Saccharomycotina, 4 other fungi, 5 everything else (an unrecognised organism included)."""
    tax = weights.taxonomy
    name = " ".join(organism.split()).casefold()
    if name in {s.casefold() for s in tax["species_names"]}:
        canonical = {_norm(s) for s in tax["canonical_strains"]}
        return 0 if strain is not None and _norm(strain) in canonical else 1
    genus = name.split(" ")[0] if name else ""
    if genus == tax["genus"].casefold():
        return 2
    if genus in {g.casefold() for g in tax["saccharomycotina_genera"]}:
        return 3
    if genus in {g.casefold() for g in tax["fungi_genera"]}:
        return 4
    return 5
