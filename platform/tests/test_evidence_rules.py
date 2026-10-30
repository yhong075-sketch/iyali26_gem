"""Contradiction engine: the twelve required cases of spec §5 plus the rules they rest on."""

from evidence.converters import Builder
from evidence.id_map import IdMap
from evidence.rules import classify_all, classify_pair
from evidence.schemas import ClaimLevel as L, Condition, PairClass, Predicate as P
from evidence.taxonomy import taxon_distance

from evidence_helpers import W, make_claim


def cls(a, b):
    rel = classify_pair(a, b, W)
    return None if rel is None else (rel.pair_class, rel.reason)


def test_01_ferment_vs_crude_activity_not_comparable():
    a = make_claim("a", P.FLUX_ACTIVE_UNDER, "negative", L.PHYSIOLOGY, reaction="ADH",
                   condition=Condition(oxygen="aerobic"))
    b = make_claim("b", P.REACTION_EXISTS, "positive", L.ACTIVITY_CRUDE, reaction="ADH")
    assert cls(a, b) is None
    assert classify_all([a, b], W) == []


def test_02_xylose_adaptation_scoped_apart_by_condition():
    a = make_claim("a", P.FLUX_ACTIVE_UNDER, "negative", L.PHYSIOLOGY, reaction="XR",
                   condition=Condition(carbon_source="xylose", adaptation="none"))
    b = make_claim("b", P.FLUX_ACTIVE_UNDER, "positive", L.PHYSIOLOGY, reaction="XR",
                   condition=Condition(carbon_source="xylose", adaptation="after serial passage on xylose"))
    assert cls(a, b) == (PairClass.SCOPED_APART, "condition")


def test_03_yali0_and_yali1_ids_agree_only_after_normalisation():
    id_map = IdMap.from_pairs([("YALI1A00001g", "YALI0A00002g")])
    builder = Builder(W, id_map, "test")
    old = builder.gene_subject("YALI0A00002g", "R1")
    new = builder.gene_subject("YALI1_A00001g", "R1")
    assert old.gene_id == new.gene_id == "YALI1A00001g"
    a = make_claim("a", P.GENE_ENCODES, subject=old, exclusive=True)
    b = make_claim("b", P.GENE_ENCODES, subject=new, exclusive=True)
    assert cls(a, b) == (PairClass.AGREE_INDEPENDENT, None)
    # Without normalisation the same assertion looks like two exclusive genes: a contradiction.
    raw_a = make_claim("a", P.GENE_ENCODES, gene="YALI0A00002g", exclusive=True)
    raw_b = make_claim("b", P.GENE_ENCODES, gene="YALI1A00001g", exclusive=True)
    assert cls(raw_a, raw_b) == (PairClass.CONTRADICTION_STRONG, "exclusive")


def test_04_ko_no_effect_vs_heterologous_activity_redundancy():
    a = make_claim("a", P.GENE_ENCODES, "negative", L.PHENOTYPE_KO, gene="YALI1A00001g")
    b = make_claim("b", P.GENE_ENCODES, "positive", L.ACTIVITY_HETEROLOGOUS, gene="YALI1A00001g")
    assert cls(a, b) == (PairClass.SCOPED_APART, "redundancy_possible")


def test_05_two_genes_neither_exclusive_both_support():
    a = make_claim("a", P.GENE_ENCODES, gene="YALI1A00001g")
    b = make_claim("b", P.GENE_ENCODES, gene="YALI1B00002g")
    assert cls(a, b) == (PairClass.AGREE_INDEPENDENT, None)


def test_06_peroxisome_vs_mitochondria_direct_strong():
    a = make_claim("a", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_pe"})
    b = make_claim("b", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_mi"})
    assert cls(a, b) == (PairClass.CONTRADICTION_STRONG, None)


def test_07_direct_vs_predicted_localisation_level_mismatch():
    a = make_claim("a", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_pe"})
    b = make_claim("b", P.LOCALIZED_IN, level=L.HOMOLOGY_ANNOTATION, compartments={"C_cy"},
                   method="targeting-signal prediction")
    assert cls(a, b) == (PairClass.SCOPED_APART, "level_mismatch")


def test_08_cerevisiae_ortholog_not_comparable():
    a = make_claim("a", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_cy"},
                   organism="Saccharomyces cerevisiae", strain=None)
    b = make_claim("b", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_pe"})
    assert a.taxon_distance == 3
    assert cls(a, b) is None


def test_09_review_restating_paper_is_derivative():
    paper = make_claim("x", P.REACTION_EXISTS, source="S-paper", year=2001)
    review = make_claim("r", P.REACTION_EXISTS, level=L.AUTHOR_STATEMENT_TRACEABLE, source="S-review", year=2010,
                        original=False, cited=("S-paper",))
    assert cls(paper, review) == (PairClass.AGREE_DERIVATIVE, None)


def test_10_explicit_rebuttal_supersedes():
    early = make_claim("e", P.REACTION_EXISTS, "positive", L.ACTIVITY_CRUDE, source="S-early", year=1990)
    late = make_claim("l", P.REACTION_EXISTS, "negative", L.ACTIVITY_CRUDE, source="S-late", year=2015,
                      rebuts=("S-early",))
    pair_class, reason = cls(early, late)
    assert pair_class == PairClass.SUPERSEDED
    assert reason == "l rebuts S-early"


def test_11_no_gene_annotated_vs_activity_is_anomaly():
    a = make_claim("a", P.REACTION_EXISTS, "negative", L.HOMOLOGY_ANNOTATION)
    b = make_claim("b", P.REACTION_EXISTS, "positive", L.ACTIVITY_CRUDE)
    assert cls(a, b) == (PairClass.ANOMALY_CROSS_LEVEL, None)


def test_12_speculative_claim_never_paired():
    a = make_claim("a", P.REACTION_EXISTS, "negative", strength="speculative")
    b = make_claim("b", P.REACTION_EXISTS, "positive")
    assert cls(a, b) is None
    assert classify_all([a, b], W) == []


# ── supporting rules ─────────────────────────────────────────────────────────────────────────

def test_hedged_contradiction_is_weak():
    a = make_claim("a", P.REACTION_EXISTS, "negative", strength="hedged")
    b = make_claim("b", P.REACTION_EXISTS, "positive")
    assert cls(a, b) == (PairClass.CONTRADICTION_WEAK, None)


def test_strain_specific_result_scoped_apart():
    a = make_claim("a", P.REACTION_EXISTS, "negative", strain="H222", strain_specific=True)
    b = make_claim("b", P.REACTION_EXISTS, "positive", strain="W29")
    assert cls(a, b) == (PairClass.SCOPED_APART, "strain")


def test_essentiality_keys_on_medium_and_marker():
    a = make_claim("a", P.ESSENTIAL_UNDER, "positive", L.PHENOTYPE_KO, gene="YALI1A00001g",
                   condition=Condition(medium="SD-Leu", marker_complementation="LEU2"))
    b = make_claim("b", P.ESSENTIAL_UNDER, "negative", L.PHENOTYPE_KO, gene="YALI1A00001g",
                   condition=Condition(medium="YPD"))
    c = make_claim("c", P.ESSENTIAL_UNDER, "negative", L.PHENOTYPE_KO, gene="YALI1A00001g",
                   condition=Condition(medium="SD-Leu", marker_complementation="none"))
    d = make_claim("d", P.ESSENTIAL_UNDER, "negative", L.PHENOTYPE_KO, gene="YALI1A00001g",
                   condition=Condition(medium="sd-leu"))
    assert cls(a, b) == (PairClass.SCOPED_APART, "condition")
    assert cls(a, c) == (PairClass.SCOPED_APART, "condition")
    assert cls(a, d) == (PairClass.CONTRADICTION_STRONG, None)


def test_physiological_vs_biochemical_direction_scoped():
    a = make_claim("a", P.DIRECTIONALITY, level=L.FLUX_13C, direction="forward", dlevel="physiological")
    b = make_claim("b", P.DIRECTIONALITY, level=L.ACTIVITY_PURIFIED, direction="reverse", dlevel="biochemical")
    assert cls(a, b) == (PairClass.SCOPED_APART, "direction_level")


def test_reversible_vs_forward_needs_both_assertive():
    a = make_claim("a", P.DIRECTIONALITY, direction="reversible")
    b = make_claim("b", P.DIRECTIONALITY, direction="forward")
    c = make_claim("c", P.DIRECTIONALITY, direction="forward", strength="hedged")
    assert cls(a, b) == (PairClass.CONTRADICTION_STRONG, None)
    assert cls(a, c) == (PairClass.AGREE_INDEPENDENT, None)


def test_cofactor_subset_agrees_disjoint_contradicts():
    prefers = make_claim("a", P.COFACTOR, roles={"donor": {"m41"}})
    both = make_claim("b", P.COFACTOR, roles={"donor": {"m41", "m122"}})
    nad = make_claim("c", P.COFACTOR, roles={"donor": {"m122"}})
    predicted = make_claim("d", P.COFACTOR, level=L.HOMOLOGY_ANNOTATION, roles={"donor": {"m122"}})
    assert cls(prefers, both) == (PairClass.AGREE_INDEPENDENT, None)
    assert cls(prefers, nad) == (PairClass.CONTRADICTION_STRONG, None)
    assert cls(prefers, predicted) == (PairClass.SCOPED_APART, "level_mismatch")


def test_shared_authors_mark_independence_class():
    a = make_claim("a", P.REACTION_EXISTS, authors=("Smith J", "Lee K", "Ng P"), year=2010)
    b = make_claim("b", P.REACTION_EXISTS, authors=("Lee K", "Ng P"), year=2014)
    c = make_claim("c", P.REACTION_EXISTS, authors=("Lee K", "Ng P"), year=2020)
    assert classify_pair(a, b, W).independence_shared
    assert not classify_pair(a, c, W).independence_shared


def test_unresolved_subject_never_paired():
    id_map = IdMap.from_pairs([("YALI1A00001g", "YALI0A00002g"), ("YALI1A00003g", "YALI0A00002g")])
    subject = Builder(W, id_map, "test").gene_subject("YALI0A00002g", "R1")
    assert subject.status == "unresolved"
    a = make_claim("a", P.GENE_ENCODES, subject=subject)
    b = make_claim("b", P.GENE_ENCODES, gene="YALI1A00001g")
    assert cls(a, b) is None


def test_taxon_distance_computed_in_code():
    assert taxon_distance("Yarrowia lipolytica", "W29", W) == 0
    assert taxon_distance("Yarrowia lipolytica", "PO1f::Cas9", W) == 0
    assert taxon_distance("Yarrowia lipolytica", "po1F", W) == 0
    assert taxon_distance("Yarrowia lipolytica", None, W) == 1          # null strain is never W29
    assert taxon_distance("Yarrowia lipolytica", "CLIB122", W) == 1
    assert taxon_distance("Candida lipolytica", None, W) == 1           # historical name
    assert taxon_distance("Yarrowia deformans", None, W) == 2
    assert taxon_distance("Saccharomyces cerevisiae", "S288C", W) == 3
    assert taxon_distance("Aspergillus niger", None, W) == 4
    assert taxon_distance("Escherichia coli", None, W) == 5


def test_id_map_cardinality():
    id_map = IdMap.from_pairs([("YALI1A00001g", "YALI0A00001g"),
                               ("YALI1B00001g", "YALI0B00001g"), ("YALI1B00001g", "YALI0B00002g"),
                               ("YALI1C00001g", "YALI0C00001g"), ("YALI1C00002g", "YALI0C00001g")])
    assert id_map.normalize("YALI0A00001g").canonical == "YALI1A00001g"
    assert id_map.normalize("YALI0B00001g").status == "unresolved"
    assert id_map.normalize("YALI1B00001g").cardinality == "1:n"
    assert id_map.normalize("YALI0C00001g").cardinality == "1:n"
    assert id_map.normalize("YALI1C00001g").cardinality == "n:1"
    assert id_map.normalize("YALI0F99999g").status == "unmapped"
    assert id_map.normalize("YALIUNK1").status == "unmapped"
