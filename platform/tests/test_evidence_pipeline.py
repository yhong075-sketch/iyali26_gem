"""L4 aggregation, L5 verdicts and the backward pass on small in-memory claim sets (no LLM, no model file)."""

import pytest

from evidence import store
from evidence.converters import Builder, find_quote
from evidence.id_map import IdMap
from evidence.pipeline import ModelIdentity, leaf_claims, run, trace
from evidence.schemas import ClaimLevel as L, Condition, ModelElementState, Predicate as P, SourceKind

from evidence_helpers import W, make_claim, support_for

MODEL = ModelIdentity(path="model/candidates/test.xml", commit="abc123", sha256="0" * 64)


def state(rid="R1", genes=("YALI1A00001g",), lb=0.0, ub=1000.0, compartments=("C_cy",), groups=("Glycolysis",)):
    return ModelElementState(reaction_id=rid, name=rid, gene_ids=tuple(genes), gpr=" or ".join(genes),
                             compartments=compartments, lower_bound=lb, upper_bound=ub, groups=groups,
                             species_base_ids=("m41", "m10"))


def build(claims, states=None, kind=SourceKind.ARTICLE):
    sources, passages = support_for(claims, kind)
    return run(sources, passages, claims, [], states or {"R1": state()}, W, MODEL)


def verdict(bundle, element, rid="R1"):
    return next(v for v in bundle.verdicts.values() if v.reaction_id == rid and v.element == element)


def assertion(bundle, subject, predicate):
    return next(a for a in bundle.assertions.values() if a.subject_key == subject and a.predicate == predicate)


def test_direct_activity_supports_existence_with_score_4():
    b = build([make_claim("a", P.REACTION_EXISTS, level=L.ACTIVITY_PURIFIED)])
    a = assertion(b, "R1", P.REACTION_EXISTS)
    assert (a.score, a.n_independent, a.conflict_status) == (4, 1, "NONE")
    assert verdict(b, "existence").status == "SUPPORTED"


def test_contradiction_resolved_by_directness():
    claims = [make_claim("pos", P.REACTION_EXISTS, "positive", L.ACTIVITY_PURIFIED),
              make_claim("neg", P.REACTION_EXISTS, "negative", L.PHENOTYPE_KO)]
    b = build(claims)
    a = assertion(b, "R1", P.REACTION_EXISTS)
    assert a.conflict_status == "CONTRADICTION_STRONG"
    assert (a.resolved, a.resolved_by, a.score) == (True, "R2", 4)
    assert a.sides[a.winner].polarity == "positive" and a.standing == (a.winner,)


def test_unresolved_contradiction_caps_score_and_conflicts():
    claims = [make_claim("pos", P.REACTION_EXISTS, "positive", L.ACTIVITY_PURIFIED, year=2010),
              make_claim("neg", P.REACTION_EXISTS, "negative", L.ACTIVITY_PURIFIED, year=2010)]
    b = build(claims)
    a = assertion(b, "R1", P.REACTION_EXISTS)
    assert (a.resolved, a.flag, a.score) == (False, "CONFLICTING", 2)
    assert a.resolution_trace[-2][0] == "UNRESOLVED"
    v = verdict(b, "existence")
    assert v.status == "CONFLICTING"
    assert v.priority == W.priority["uncertainty"]["CONFLICTING"]


def test_recency_breaks_only_a_full_tie():
    claims = [make_claim("old", P.REACTION_EXISTS, "positive", L.ACTIVITY_PURIFIED, year=2001),
              make_claim("new", P.REACTION_EXISTS, "negative", L.ACTIVITY_PURIFIED, year=2019)]
    a = assertion(build(claims), "R1", P.REACTION_EXISTS)
    assert (a.resolved_by, a.sides[a.winner].polarity) == ("R6", "negative")


def test_assertive_rebuttal_wins_regardless_of_directness():
    # R2 alone would pick the purified-enzyme claim; the assertive rebuttal wins instead.
    claims = [make_claim("e", P.REACTION_EXISTS, "positive", L.ACTIVITY_PURIFIED, source="S-early", year=1990),
              make_claim("l", P.REACTION_EXISTS, "negative", L.PHENOTYPE_KO, source="S-late", year=2015,
                         rebuts=("S-early",))]
    b = build(claims)
    a = assertion(b, "R1", P.REACTION_EXISTS)
    assert (a.resolved_by, a.sides[a.winner].polarity) == ("R5", "negative")
    v = verdict(b, "existence")
    assert (v.status, v.actions) == ("CONTRADICTED", ("propose_remove",))


def test_direct_negative_beats_model_inclusion_without_contradiction():
    claims = [make_claim("incl", P.REACTION_EXISTS, "positive", L.MODEL_INCLUSION),
              make_claim("neg", P.REACTION_EXISTS, "negative", L.ACTIVITY_PURIFIED)]
    b = build(claims)
    a = assertion(b, "R1", P.REACTION_EXISTS)
    assert a.conflict_status == "NONE"          # different bands: scoped apart, both sides stand
    assert {s.polarity for s in a.standing_sides} == {"positive", "negative"}
    assert verdict(b, "existence").status == "CONTRADICTED"


def test_model_inclusion_only_is_unsupported_inherited_and_homology_only():
    claims = [make_claim("incl", P.REACTION_EXISTS, level=L.MODEL_INCLUSION)]
    sources, passages = support_for(claims, SourceKind.PRIOR_MODEL)
    sources = {k: v.model_copy(update={"properties": {"role": "parent", "model_name": "iYli21"}})
               for k, v in sources.items()}
    b = run(sources, passages, claims, [], {"R1": state()}, W, MODEL)
    a = assertion(b, "R1", P.REACTION_EXISTS)
    assert (a.score, a.inherited_from, a.homology_only) == (1, "iYli21", True)
    v = verdict(b, "existence")
    assert v.status == "UNSUPPORTED" and "never a removal" in v.reason


def test_derivative_claims_do_not_add_independence():
    claims = [make_claim("x", P.REACTION_EXISTS, level=L.PHENOTYPE_KO, source="S-x"),
              make_claim("r", P.REACTION_EXISTS, level=L.AUTHOR_STATEMENT_TRACEABLE, source="S-r",
                         original=False, cited=("S-x",)),
              make_claim("y", P.REACTION_EXISTS, level=L.PHENOTYPE_KO, source="S-y", authors=("A B",))]
    a = assertion(build(claims), "R1", P.REACTION_EXISTS)
    assert a.n_independent == 2


def test_transfer_from_distant_species_scores_two_only_from_activity():
    act = make_claim("o", P.REACTION_EXISTS, level=L.ACTIVITY_PURIFIED, organism="Saccharomyces cerevisiae",
                     strain=None)
    assert assertion(build([act]), "R1", P.REACTION_EXISTS).score == 2
    far = make_claim("f", P.REACTION_EXISTS, level=L.ACTIVITY_PURIFIED, organism="Aspergillus niger", strain=None)
    assert assertion(build([far]), "R1", P.REACTION_EXISTS).score == 1


def test_flux_claim_derives_existence_but_negative_flux_does_not():
    pos = make_claim("f", P.FLUX_ACTIVE_UNDER, level=L.FLUX_13C, condition=Condition(carbon_source="glucose"))
    neg = make_claim("g", P.FLUX_ACTIVE_UNDER, "negative", L.FLUX_13C, condition=Condition(carbon_source="oleic acid"))
    b = build([pos, neg])
    derived = b.claims["f~exists"]
    assert (derived.predicate, derived.derived_from, derived.condition) == (P.REACTION_EXISTS, "f", None)
    assert "g~exists" not in b.claims
    assert verdict(b, "existence").status == "SUPPORTED"


def test_gpr_add_gene_and_model_as_input():
    claims = [make_claim("g1", P.GENE_ENCODES, level=L.ACTIVITY_HETEROLOGOUS, gene="YALI1A00001g"),
              make_claim("g2", P.GENE_ENCODES, level=L.PHENOTYPE_COMPLEMENTATION, gene="YALI1B00002g")]
    with_one = build(claims, {"R1": state(genes=("YALI1A00001g",))})
    v = verdict(with_one, "gpr")
    assert (v.status, v.actions) == ("PARTIAL", ("add_gene:YALI1B00002g",))
    with_both = build(claims, {"R1": state(genes=("YALI1A00001g", "YALI1B00002g"))})
    assert verdict(with_both, "gpr").status == "SUPPORTED"
    assert v.model_commit == MODEL.commit and v.model_sha256 == MODEL.sha256


def test_negative_knockout_alone_never_removes_a_gene():
    claims = [make_claim("ko", P.GENE_ENCODES, "negative", L.PHENOTYPE_KO, gene="YALI1A00001g")]
    v = verdict(build(claims), "gpr")
    assert v.status == "UNSUPPORTED" and not any(a.startswith("remove_gene") for a in v.actions)


def test_exclusive_claim_against_isozymes_conflicts():
    claims = [make_claim("x", P.GENE_ENCODES, level=L.ACTIVITY_PURIFIED, gene="YALI1A00001g", exclusive=True),
              make_claim("y", P.GENE_ENCODES, level=L.ACTIVITY_PURIFIED, gene="YALI1B00002g")]
    v = verdict(build(claims, {"R1": state(genes=("YALI1A00001g", "YALI1B00002g"))}), "gpr")
    assert v.status == "CONFLICTING"


def test_localisation_and_direction_contradict_model_state():
    claims = [make_claim("loc", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_pe"}),
              make_claim("dir", P.DIRECTIONALITY, level=L.ACTIVITY_PURIFIED, direction="forward")]
    b = build(claims, {"R1": state(lb=-1000.0, compartments=("C_cy",))})
    assert (verdict(b, "compartment").status, verdict(b, "compartment").actions) == ("CONTRADICTED", ("relocate",))
    assert (verdict(b, "bounds").status, verdict(b, "bounds").actions) == ("CONTRADICTED", ("change_bounds",))


def test_subsystem_reassign_never_contradicted():
    v = verdict(build([make_subsystem("s", "TCA cycle")]), "subsystem")
    assert (v.status, v.actions) == ("PARTIAL", ("reassign:TCA cycle",))


def make_subsystem(cid, pathway):
    from evidence.schemas import SubsystemValue
    base = make_claim(cid, P.REACTION_EXISTS, level=L.MODEL_INCLUSION)
    return base.model_copy(update={"predicate": P.SUBSYSTEM_MEMBER, "value": SubsystemValue(
        predicate=P.SUBSYSTEM_MEMBER, pathways=frozenset({pathway}), ontology="model_groups")})


def test_trace_reaches_sources_for_every_verdict():
    claims = [make_claim("pos", P.REACTION_EXISTS, level=L.ACTIVITY_PURIFIED),
              make_claim("neg", P.REACTION_EXISTS, "negative", L.PHENOTYPE_KO),
              make_claim("f", P.FLUX_ACTIVE_UNDER, level=L.FLUX_13C)]
    b = build(claims)
    for vid, v in b.verdicts.items():
        nodes = trace(b, vid)
        assert nodes[0][0] == "L5"
        if leaf_claims(b, vid):
            layers = {layer for layer, _, _ in nodes}
            assert {"L4", "L2", "L1", "L0"} <= layers


def test_load_bearing_finds_the_claim_that_carries_the_verdict():
    claims = [make_claim("strong", P.REACTION_EXISTS, level=L.ACTIVITY_PURIFIED),
              make_claim("weak", P.REACTION_EXISTS, level=L.MODEL_INCLUSION, source="S-model")]
    b = build(claims)
    v = verdict(b, "existence")
    assert b.load_bearing[v.verdict_id] == ("strong",)


def test_claim_without_verified_passage_does_not_exist():
    claims = [make_claim("a", P.REACTION_EXISTS)]
    sources, passages = support_for(claims)
    passages = {k: p.model_copy(update={"verified": False}) for k, p in passages.items()}
    with pytest.raises(ValueError, match="verified passage"):
        run(sources, passages, claims, [], {"R1": state()}, W, MODEL)


def test_quote_verification_is_whitespace_normalised_exact():
    text = "PS species were not\n   detected in the lipid extracts"
    start, end = find_quote(text, "species were not detected in")
    assert " ".join(text[start:end].split()) == "species were not detected in" and start == 3
    assert find_quote(text, "species were not found") is None
    builder = Builder(W, IdMap.from_pairs([]), "test")
    from evidence.schemas import Source
    src = Source(source_id="s", kind=SourceKind.ARTICLE, title="t", peer_reviewed=True)
    assert builder.passage(src, text, "lipid extracts of", "sec") is None
    assert builder.out.rejected[0].reason == "quote_not_found"


def test_store_is_deterministic_and_traceable(tmp_path):
    claims = [make_claim("pos", P.REACTION_EXISTS, level=L.ACTIVITY_PURIFIED),
              make_claim("neg", P.REACTION_EXISTS, "negative", L.ACTIVITY_PURIFIED, year=2020)]
    first = store.write(tmp_path / "a.sqlite", build(claims), [], W.version, [])
    second = store.write(tmp_path / "b.sqlite", build(list(reversed(claims))), [], W.version, [])
    assert first == second
    loaded = store.load(tmp_path / "a.sqlite")
    vid = next(v.verdict_id for v in loaded.verdicts.values() if v.element == "existence")
    assert any(layer == "L0" for layer, _, _ in trace(loaded, vid))


# ── review findings (M1 code review) ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("negative", [
    make_claim("hom", P.REACTION_EXISTS, "negative", L.HOMOLOGY_ANNOTATION, method="no gene annotated"),
    make_claim("gen", P.REACTION_EXISTS, "negative", L.GENOMIC_PRESENCE),
    make_claim("sce", P.REACTION_EXISTS, "negative", L.PHENOTYPE_KO, organism="Saccharomyces cerevisiae", strain=None),
])
def test_weak_or_transferred_negative_never_proposes_removal(negative):
    incl = make_claim("incl", P.REACTION_EXISTS, level=L.MODEL_INCLUSION, source="S-model")
    v = verdict(build([incl, negative]), "existence")
    assert v.status == "UNSUPPORTED" and "propose_remove" not in v.actions


def test_negative_localisation_and_direction_are_not_reported_as_support():
    claims = [make_claim("loc", P.LOCALIZED_IN, "negative", L.LOCALIZATION_DIRECT, compartments={"C_mi"}),
              make_claim("dir", P.DIRECTIONALITY, "negative", L.ACTIVITY_PURIFIED, direction="reversible")]
    b = build(claims, {"R1": state(lb=-1000.0, compartments=("C_mi",))})
    assert verdict(b, "compartment").status == "CONTRADICTED"
    assert verdict(b, "bounds").status == "CONTRADICTED"


@pytest.mark.parametrize("years", [(2024, 2010), (2010, 2024)])
def test_scoped_apart_sides_both_stand_whatever_the_year(years):
    claims = [make_claim("flux", P.DIRECTIONALITY, level=L.FLUX_13C, direction="forward", dlevel="physiological",
                         year=years[0]),
              make_claim("enz", P.DIRECTIONALITY, level=L.ACTIVITY_PURIFIED, direction="reversible",
                         dlevel="biochemical", year=years[1])]
    b = build(claims, {"R1": state(lb=-1000.0)})
    a = assertion(b, "R1", P.DIRECTIONALITY)
    assert len(a.standing) == 2
    assert verdict(b, "bounds").status == "SUPPORTED"


def test_set_values_serialise_in_a_fixed_order_under_any_hash_seed():
    import subprocess
    import sys
    code = ("from evidence.converters import Builder; from evidence.id_map import IdMap; "
            "from evidence.schemas import *; from evidence.weights import load_weights; "
            "b=Builder(load_weights(), IdMap.from_pairs([]), 't'); "
            "s=Source(source_id='s', kind='article', title='t', peer_reviewed=True); "
            "p=b.passage(s, 'in four places', 'four places', 'x'); "
            "c=b.claim(s, p, SubjectRef(kind='reaction', reaction_id='R1'), Predicate.LOCALIZED_IN, "
            "LocalizedValue(predicate=Predicate.LOCALIZED_IN, compartments=frozenset({'C_cy','C_mi','C_pe','C_er'})), "
            "'positive', ClaimLevel.LOCALIZATION_DIRECT, method='GFP'); print(c.claim_id)")
    ids = {subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                          env={"PYTHONHASHSEED": str(seed), "PYTHONPATH": "."}).stdout for seed in (1, 2, 3, 4)}
    assert len(ids) == 1


def test_exclusive_claim_alone_conflicts_with_model_isozymes():
    claims = [make_claim("x", P.GENE_ENCODES, level=L.ACTIVITY_PURIFIED, gene="YALI1A00001g", exclusive=True)]
    v = verdict(build(claims, {"R1": state(genes=("YALI1A00001g", "YALI1B00002g"))}), "gpr")
    assert v.status == "CONFLICTING"


def test_inadmissible_exclusive_claim_does_not_create_a_conflict():
    claims = [make_claim("x", P.GENE_ENCODES, level=L.ACTIVITY_CRUDE, gene="YALI1A00001g", exclusive=True),
              make_claim("y", P.GENE_ENCODES, level=L.ACTIVITY_PURIFIED, gene="YALI1B00002g")]
    b = build(claims, {"R1": state(genes=("YALI1A00001g", "YALI1B00002g"))})
    assert assertion(b, "YALI1B00002g|R1", P.GENE_ENCODES).conflict_status == "NONE"
    assert verdict(b, "gpr").status != "CONFLICTING"


def test_flux_derived_existence_counts_one_independent_class():
    claims = [make_claim("f", P.FLUX_ACTIVE_UNDER, level=L.FLUX_13C, condition=Condition(carbon_source="glucose"))]
    a = assertion(build(claims), "R1", P.REACTION_EXISTS)
    assert (a.score, a.n_independent) == (4, 1)


def test_claims_differing_only_in_organism_get_different_ids():
    from evidence.schemas import FlagValue, Source, SubjectRef
    builder = Builder(W, IdMap.from_pairs([]), "test")
    src = Source(source_id="s", kind=SourceKind.ARTICLE, title="t", peer_reviewed=True, year=2020)
    p = builder.passage(src, "activity was found", "activity was found", "sec")
    args = (src, p, SubjectRef(kind="reaction", reaction_id="R1"), P.REACTION_EXISTS,
            FlagValue(predicate=P.REACTION_EXISTS), "positive", L.ACTIVITY_PURIFIED)
    a = builder.claim(*args, method="assay")
    b = builder.claim(*args, method="assay", organism="Saccharomyces cerevisiae")
    assert a.claim_id != b.claim_id


def test_over_long_quote_is_rejected_not_truncated():
    from evidence.schemas import Source
    builder = Builder(W, IdMap.from_pairs([]), "test")
    src = Source(source_id="s", kind=SourceKind.ARTICLE, title="t", peer_reviewed=True)
    text = " ".join(f"w{i}" for i in range(24)) + " the activity was not detected"
    assert builder.passage(src, text, text, "sec") is None
    assert builder.out.rejected[0].reason == "quote_too_long"


def test_hedged_rebuttal_is_a_weak_contradiction():
    claims = [make_claim("e", P.REACTION_EXISTS, "positive", L.ACTIVITY_CRUDE, source="S-early", year=1990),
              make_claim("l", P.REACTION_EXISTS, "negative", L.ACTIVITY_CRUDE, source="S-late", year=2015,
                         rebuts=("S-early",), strength="hedged")]
    a = assertion(build(claims), "R1", P.REACTION_EXISTS)
    assert a.conflict_status == "CONTRADICTION_WEAK"
    assert a.resolved_by != "R5"            # a hedged rebuttal does not override R2–R4


def test_unknown_impact_is_stated_on_the_verdict():
    s = state().model_copy(update={"impact_known": False})
    v = verdict(build([make_claim("a", P.REACTION_EXISTS)], {"R1": s}), "existence")
    assert v.priority_note and "impact unknown" in v.priority_note


def test_uniprot_and_ncbi_ids_normalise_through_model_annotations():
    id_map = IdMap.from_pairs([])
    id_map.add_model_genes({"YALI1A00001g": {"uniprot": "Q6C000", "ncbigene": "2906000"}})
    assert id_map.normalize("Q6C000").canonical == "YALI1A00001g"
    assert id_map.normalize("2906000").canonical == "YALI1A00001g"
    assert id_map.normalize("P99999").status == "unmapped"


def test_parent_reaction_with_different_species_is_not_claimed(tmp_path):
    from evidence.converters import parent_model_claims
    xml = tmp_path / "parent.xml"
    xml.write_text('<model><reaction id="R1" name="enzyme"><listOfReactants>'
                   '<speciesReference species="m1__91__C_cy__93__"/></listOfReactants></reaction></model>')
    lineage = {"iYli21": {"title": "p", "year": 2022, "parent": None, "publication": "", "origin_tags": ""}}
    same = parent_model_claims(xml, ["R1"], W, IdMap.from_pairs([]), lineage, {"R1": {"m1[C_cy]"}})
    other = parent_model_claims(xml, ["R1"], W, IdMap.from_pairs([]), lineage, {"R1": {"m1[C_mi]"}})
    assert len(same.claims) == 1 and not other.claims
    assert other.rejected[0].reason == "parent_reaction_differs"
    extra_proton = parent_model_claims(xml, ["R1"], W, IdMap.from_pairs([]), lineage, {"R1": {"m1[C_cy]", "m10[C_cy]"}},
                                       frozenset({"m10[C_cy]"}))
    assert len(extra_proton.claims) == 1


# ── re-review findings ───────────────────────────────────────────────────────────────────────

def test_evidence_that_a_closed_reaction_runs_contradicts_its_bounds():
    claims = [make_claim("dir", P.DIRECTIONALITY, level=L.ACTIVITY_PURIFIED, direction="forward")]
    v = verdict(build(claims, {"R1": state(lb=0.0, ub=0.0)}), "bounds")
    assert (v.status, v.actions) == ("CONTRADICTED", ("change_bounds",))


def test_negative_physiological_direction_does_not_remove_a_capability():
    phys = [make_claim("n", P.DIRECTIONALITY, "negative", L.FLUX_13C, direction="reverse", dlevel="physiological")]
    bio = [make_claim("n", P.DIRECTIONALITY, "negative", L.ACTIVITY_PURIFIED, direction="reverse")]
    assert verdict(build(phys, {"R1": state(lb=-1000.0)}), "bounds").status == "SUPPORTED"
    assert verdict(build(bio, {"R1": state(lb=-1000.0)}), "bounds").status == "CONTRADICTED"


def test_accession_of_an_unmapped_model_gene_does_not_resolve():
    id_map = IdMap.from_pairs([])
    id_map.add_model_genes({"YALI0Z99999g": {"uniprot": "Q6XXXX"}})
    assert id_map.normalize("Q6XXXX").status == "unmapped"


def test_weak_exclusive_claim_or_and_complex_does_not_conflict():
    weak = [make_claim("x", P.GENE_ENCODES, level=L.HOMOLOGY_ANNOTATION, gene="YALI1A00001g", exclusive=True,
                       strength="hedged")]
    isozymes = state(genes=("YALI1A00001g", "YALI1B00002g"))
    assert verdict(build(weak, {"R1": isozymes}), "gpr").status != "CONFLICTING"
    strong = [make_claim("x", P.GENE_ENCODES, level=L.ACTIVITY_PURIFIED, gene="YALI1A00001g", exclusive=True)]
    complex_ = isozymes.model_copy(update={"gpr": "YALI1A00001g and YALI1B00002g"})
    assert verdict(build(strong, {"R1": complex_}), "gpr").status != "CONFLICTING"


def test_identical_rejections_keep_one_row_each():
    from evidence.schemas import Source
    builder = Builder(W, IdMap.from_pairs([]), "test")
    src = Source(source_id="s", kind=SourceKind.ARTICLE, title="t", peer_reviewed=True)
    builder.passage(src, "some text", "missing words", "sec")
    builder.passage(src, "some text", "missing words", "sec")
    assert len({r.reject_id for r in builder.out.rejected}) == 2


def test_resolved_contradiction_scores_from_its_winner():
    claims = [make_claim("pe", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_pe"}, year=2020),
              make_claim("mi", P.LOCALIZED_IN, level=L.LOCALIZATION_DIRECT, compartments={"C_mi"}, year=2010),
              make_claim("cy", P.LOCALIZED_IN, level=L.HOMOLOGY_ANNOTATION, compartments={"C_pe", "C_cy"},
                         method="prediction")]
    a = assertion(build(claims), "R1", P.LOCALIZED_IN)
    assert a.winner is not None and a.score == a.sides[a.winner].score
