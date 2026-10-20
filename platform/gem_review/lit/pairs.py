"""Select the reactions and gene–reaction pairs a literature review runs on."""

import csv
from dataclasses import dataclass

from scripts.gem_annotate.model_layout import MODEL

# Tier → curated subsystem list (model/curation/gem_review/). Only tier 1 is defined so far.
TIER_FILES = {1: "tier1_lipid_subsystems.csv"}


@dataclass(frozen=True, slots=True)
class Pair:
    gene: str
    reaction_id: str
    reaction_name: str
    subsystem: str
    ec_codes: tuple[str, ...]
    tier: int | None


def load_tier_subsystems(tier: int) -> set[str]:
    """Read the subsystems marked include=yes for this tier."""
    if tier not in TIER_FILES:
        raise ValueError(f"No subsystem list is defined for tier {tier}; known: {sorted(TIER_FILES)}")
    path = MODEL.curation_file(TIER_FILES[tier])
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    subsystems = {row["subsystem"] for row in rows if row["include"].strip().lower() == "yes"}
    if not subsystems:
        raise ValueError(f"{path} marks no subsystem as include=yes")
    return subsystems


def reaction_subsystems(reaction) -> set[str]:
    return {part.strip() for part in (reaction.subsystem or "").split(";") if part.strip()}


def reaction_ec_codes(reaction) -> tuple[str, ...]:
    codes = reaction.annotation.get("ec-code", [])
    if isinstance(codes, str):
        codes = [codes]
    return tuple(sorted({code for code in codes if code}))


def select_reactions(model, *, reaction_ids=None, tier: int | None = None) -> list:
    """Explicit reaction IDs and/or a tier. Refuses to run on the whole model."""
    if not reaction_ids and tier is None:
        raise ValueError("Refusing to run on the whole model: pass reaction IDs or a tier.")
    chosen = []
    if reaction_ids:
        missing = [rid for rid in reaction_ids if rid not in model.reactions]
        if missing:
            raise KeyError(f"Reactions not in model: {missing}")
        chosen.extend(model.reactions.get_by_id(rid) for rid in reaction_ids)
    if tier is not None:
        wanted = load_tier_subsystems(tier)
        chosen.extend(r for r in model.reactions if reaction_subsystems(r) & wanted)
    seen: set[str] = set()
    unique = []
    for reaction in chosen:
        if reaction.id not in seen:
            seen.add(reaction.id)
            unique.append(reaction)
    return unique


def list_pairs(reactions, tier: int | None = None) -> list[Pair]:
    """One Pair per (gene, reaction); reactions without a GPR yield nothing."""
    pairs = []
    for reaction in reactions:
        for gene in sorted(reaction.genes, key=lambda g: g.id):
            pairs.append(Pair(gene.id, reaction.id, reaction.name or "", reaction.subsystem or "",
                              reaction_ec_codes(reaction), tier))
    return pairs


def geneless(reactions) -> list[str]:
    return [reaction.id for reaction in reactions if not reaction.genes]
