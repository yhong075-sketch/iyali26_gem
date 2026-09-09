"""Existing Q9 build chain, migrated from the saved 2026-09-07 pipeline.

Source identity and integration scope: data/quinone_pipeline_provenance.json.
This chain runs in every CoQ9 curation mode, including off.
"""

import copy
import logging
import math
from cobra.core.gene import GPR

from .coq9 import boolean_key

logger = logging.getLogger(__name__)

_LEGACY_EXTERNAL_NDH2_GPR = 'YALI1A21711g and YALI1B00908g and YALI1B19507g and YALI1B26679g and YALI1C04281g and YALI1D00766g and YALI1D06302g and YALI1D07089g and YALI1D09203g and YALI1D18037g and YALI1D24109g and YALI1D32550g and YALI1E27218g and YALI1E37603g and YALI1F01456g and YALI1F03371g and YALI1F09003g and YALI1F22993g and YALI1F24343g and YALI1M00056g and YALI1M00064g and YALI1M00296g and YALI1M00335g and YALI1M00338r and YALI1M00458g and YALIfMp29 and YALI1A21711g and YALI1B00908g and YALI1B19507g and YALI1B26679g and YALI1C04281g and YALI1D00766g and YALI1D06302g and YALI1D07089g and YALI1D09203g and YALI1D18037g and YALI1D24109g and YALI1D32550g and YALI1E06573g and YALI1E27218g and YALI1E37603g and YALI1F01456g and YALI1F03371g and YALI1F09003g and YALI1F22993g and YALI1F24343g and YALI1F32476g'


_R570_NDH2_GENE = "YALI1F32476g"


def _upsert_gene_annotation(gene, *, name: str, annotation: dict) -> None:
    """Set curated identity fields without discarding prior identifier mapping."""
    gene.name = name
    current = dict(gene.annotation) if isinstance(gene.annotation, dict) else {}
    current.update(annotation)
    gene.annotation = current


def _stoichiometry_by_metabolite_id(reaction) -> dict[str, float]:
    return {met.id: float(coefficient) for met, coefficient in reaction.metabolites.items()}


def correct_external_ndh2_gpr_and_remove_duplicate(model) -> int:
    """Correct R570 to the verified external NDH2 and remove duplicate R2063.

    The function verifies reaction identity and bounds before deleting R2063.
    R1889 (complex I) is intentionally untouched because a structurally present
    subunit is not automatically a required Boolean GPR component.
    """
    try:
        r570 = model.reactions.get_by_id("R570")
    except KeyError as exc:
        raise ValueError("R570 is required for the external-NDH2 correction") from exc

    target_rule = _R570_NDH2_GENE
    old_rule = r570.gene_reaction_rule.strip()
    if boolean_key(r570.gpr.body) not in {
        boolean_key(GPR.from_string(_LEGACY_EXTERNAL_NDH2_GPR).body),
        boolean_key(GPR.from_string(target_rule).body),
    }:
        raise ValueError(f"R570 has an unexpected GPR; refusing replacement: {old_rule!r}")
    # Validate the duplicate before any mutation, including gene annotations.
    if "R2063" in model.reactions:
        duplicate = model.reactions.get_by_id("R2063")
        if (_stoichiometry_by_metabolite_id(duplicate) != _stoichiometry_by_metabolite_id(r570)
                or duplicate.bounds != r570.bounds
                or boolean_key(duplicate.gpr.body) not in {
                    boolean_key(GPR.from_string(_LEGACY_EXTERNAL_NDH2_GPR).body),
                    boolean_key(GPR.from_string(target_rule).body)}):
            raise ValueError("R2063 differs from the reviewed external NDH2 duplicate")
    if old_rule != target_rule:
        r570.gene_reaction_rule = target_rule
        gpr_changed = 1
        logger.info("  Reviewed GPR correction: R570 = %s (NDH2)", target_rule)
    else:
        gpr_changed = 0

    gene = model.genes.get_by_id(_R570_NDH2_GENE)
    _upsert_gene_annotation(
        gene,
        name="NDH2",
        annotation={
            "sbo": "SBO:0000243",
            "uniprot": "F2Z699",
            "ec-code": "1.6.5.9",
        },
    )
    notes = dict(r570.notes) if isinstance(r570.notes, dict) else {}
    notes["curated_gpr_correction"] = (
        "YALI1F32476g (NDH2), external alternative NADH:ubiquinone "
        "oxidoreductase; experimentally verified in Yarrowia lipolytica: "
        "https://pubmed.ncbi.nlm.nih.gov/11719558/"
    )
    notes["complex_i_scope"] = (
        "R1889 has no GPR by design in this patch; see "
        "docs/curation/complex_i_gpr_evidence.csv before assigning a complex-I AND rule."
    )
    r570.notes = notes

    try:
        r2063 = model.reactions.get_by_id("R2063")
    except KeyError:
        return gpr_changed
    if _stoichiometry_by_metabolite_id(r570) != _stoichiometry_by_metabolite_id(r2063):
        raise ValueError("R2063 is not stoichiometrically identical to R570; refusing removal")
    if r570.bounds != r2063.bounds:
        raise ValueError("R2063 bounds differ from R570; refusing duplicate removal")
    model.remove_reactions([r2063], remove_orphans=False)
    logger.info("  Reviewed duplicate removal: R2063 duplicates R570 external NDH2")
    return gpr_changed + 1


_SPURIOUS_QUINONE_REACTION_GPRS = {
    "R189": "YALI1D17983g and YALI1B21088g",
    "R2242": "YALI1E01159g",
    "R2247": "YALI1E07601g",
    "R2248": "",
    "R2249": "",
    "R2250": (
        "(YALI1E11415g and YALI1B21088g) or "
        "(YALI1E16694g and YALI1E33302g) or YALI1D21543g or "
        "(YALI1D17983g and YALI1B21088g)"
    ),
}


_REVIEWED_QUINONE_REACTION_VARIANTS = {
    "R189": (
        {
            "stoichiometry": {
                "m203[C_cy]": 1.0,
                "m366[C_cy]": -1.0,
                "m367[C_cy]": -1.0,
                "m368[C_cy]": 1.0,
            },
            "bounds": (-1000.0, 1000.0),
            "compartments": frozenset({"C_cy"}),
            "reversible": True,
        },
        {
            "stoichiometry": {
                "m10[C_cy]": 2.0,
                "m203[C_cy]": 1.0,
                "m366[C_cy]": -1.0,
                "m367[C_cy]": -1.0,
                "m368[C_cy]": 1.0,
            },
            "bounds": (-1000.0, 1000.0),
            "compartments": frozenset({"C_cy"}),
            "reversible": True,
        },
    ),
    "R2242": (
        {
            "stoichiometry": {
                "m1923[C_nu]": -1.0,
                "m1924[C_nu]": 1.0,
                "m1925[C_nu]": -1.0,
                "m1926[C_nu]": 1.0,
                "m627[C_nu]": 1.0,
            },
            "bounds": (0.0, 1000.0),
            "compartments": frozenset({"C_nu"}),
            "reversible": False,
        },
    ),
    "R2247": (
        {
            "stoichiometry": {
                "m10[C_cy]": -1.0,
                "m1928[C_cy]": 1.0,
                "m1929[C_cy]": -1.0,
                "m82[C_cy]": 1.0,
            },
            "bounds": (0.0, 1000.0),
            "compartments": frozenset({"C_cy"}),
            "reversible": False,
        },
    ),
    "R2248": (
        {
            "stoichiometry": {
                "m109[C_cy]": -0.5,
                "m1927[C_cy]": 1.0,
                "m1928[C_cy]": -1.0,
            },
            "bounds": (0.0, 1000.0),
            "compartments": frozenset({"C_cy"}),
            "reversible": False,
        },
    ),
    "R2249": (
        {
            "stoichiometry": {
                "m1923[C_nu]": 1.0,
                "m1927[C_cy]": -1.0,
            },
            "bounds": (-1000.0, 1000.0),
            "compartments": frozenset({"C_cy", "C_nu"}),
            "reversible": True,
        },
    ),
    "R2250": (
        {
            "stoichiometry": {
                "m1929[C_cy]": 1.0,
                "m1930[C_cy]": -1.0,
                "m203[C_cy]": 1.0,
                "m366[C_cy]": -1.0,
            },
            "bounds": (0.0, 1000.0),
            "compartments": frozenset({"C_cy"}),
            "reversible": False,
        },
    ),
}


_RETAINED_COQ2_REACTION = "R407"


_RETAINED_COQ2_GENE = "YALI1F08349g"


_RETAINED_COQ3_GENE = "YALI1B20835g"


_REJECTED_QUINONE_ORPHAN_METABOLITES = {
    "m367[C_cy]",
    "m368[C_cy]",
    "m1923[C_nu]",
    "m1924[C_nu]",
    "m1927[C_cy]",
    "m1928[C_cy]",
    "m1929[C_cy]",
    "m1930[C_cy]",
}


def remove_spurious_quinone_branches(model) -> int:
    """Remove six inert, mis-annotated quinone reactions as one atomic patch.

    The patch fails closed on a partially removed branch, an unexpected GPR,
    changed reaction signature, or loss/change of the retained mitochondrial
    COQ2 reaction. Branch-only metabolites are removed, but orphan genes are
    deliberately retained: five are in the positive-only essentiality
    reference and must remain visible as unresolved FNs rather than silently
    leaving the evaluation denominator. The operation is idempotent and
    returns the number of reactions removed.
    """

    present = {
        reaction_id
        for reaction_id in _SPURIOUS_QUINONE_REACTION_GPRS
        if reaction_id in model.reactions
    }
    expected = set(_SPURIOUS_QUINONE_REACTION_GPRS)
    if present and present != expected:
        missing = sorted(expected - present)
        raise ValueError(
            "Spurious quinone branch is only partially present; refusing an "
            f"atomic cleanup (missing {missing})"
        )

    try:
        retained_coq2 = model.reactions.get_by_id(_RETAINED_COQ2_REACTION)
    except KeyError as exc:
        raise ValueError(
            f"{_RETAINED_COQ2_REACTION} is required before removing quinone duplicates"
        ) from exc
    if retained_coq2.gene_reaction_rule.strip() != _RETAINED_COQ2_GENE:
        raise ValueError(
            f"{_RETAINED_COQ2_REACTION} has unexpected GPR "
            f"{retained_coq2.gene_reaction_rule!r}"
        )
    retained_metabolites = {met.id for met in retained_coq2.metabolites}
    retained_markers = {
        "m138[C_mi]",
        "m640[C_mi]",
        "m641[C_mi]",
        "m204[C_mi]",
    }
    if not retained_markers <= retained_metabolites:
        raise ValueError(
            f"{_RETAINED_COQ2_REACTION} no longer matches the mitochondrial COQ2 reaction"
        )

    # Preserve the two real mitochondrial candidates and make their evidence
    # status explicit. These are curated annotations, not direct Yarrowia
    # knockout/biochemical validation of the current model reactions.
    coq2 = model.genes.get_by_id(_RETAINED_COQ2_GENE)
    _upsert_gene_annotation(
        coq2,
        name="COQ2",
        annotation={
            "sbo": "SBO:0000243",
            "uniprot": ["A0A1H6PM88", "Q6C2S2"],
            "ncbigene": "2907969",
            "kegg.genes": "yli:2907969",
            "refseq": "XP_505040.1",
            "ec-code": "2.5.1.39",
        },
    )
    coq3 = model.genes.get_by_id(_RETAINED_COQ3_GENE)
    _upsert_gene_annotation(
        coq3,
        name="COQ3",
        annotation={
            "sbo": "SBO:0000243",
            "uniprot": ["A0A1D8N802", "Q6CEG2"],
            "ncbigene": "2907025",
            "kegg.genes": "yli:2907025",
            "refseq": "XP_500950.3",
            "ec-code": ["2.1.1.64", "2.1.1.114"],
        },
    )
    retained_notes = (
        dict(retained_coq2.notes)
        if isinstance(retained_coq2.notes, dict)
        else {}
    )
    retained_notes["curated_quinone_branch_cleanup"] = (
        "Retained mitochondrial COQ2 reaction. Removed inert/mis-annotated "
        "R189 and R2242/R2247/R2248/R2249/R2250 branches; see "
        "docs/curation/quinone_branch_cleanup.md."
    )
    retained_notes["gpr_evidence_status"] = "curated_annotation"
    if model.metabolites.get_by_id("m640[C_mi]").formula == "C45H76O7P2":
        retained_notes.pop("remaining_chain_length_gate", None)
        retained_notes["coq9_chain_status"] = (
            "The formal pipeline has replaced the legacy CoQ6 identities with "
            "the curated, mass-balanced CoQ9 main-chain representation."
        )
    else:
        retained_notes["remaining_chain_length_gate"] = (
            "The legacy main route uses hexaprenyl/CoQ6 intermediates; native "
            "Yarrowia CoQ9 chemistry is not claimed repaired by this patch."
        )
    retained_coq2.notes = retained_notes

    if not present:
        return 0

    reactions = []
    for reaction_id, expected_gpr in _SPURIOUS_QUINONE_REACTION_GPRS.items():
        reaction = model.reactions.get_by_id(reaction_id)
        if reaction.gene_reaction_rule.strip() != expected_gpr:
            raise ValueError(
                f"{reaction_id} has unexpected GPR; refusing quinone cleanup: "
                f"{reaction.gene_reaction_rule!r}"
            )
        actual_signature = {
            "stoichiometry": _stoichiometry_by_metabolite_id(reaction),
            "bounds": tuple(float(value) for value in reaction.bounds),
            "compartments": frozenset(
                metabolite.compartment for metabolite in reaction.metabolites
            ),
            "reversible": bool(reaction.reversibility),
        }
        reviewed_variants = _REVIEWED_QUINONE_REACTION_VARIANTS[reaction_id]
        if actual_signature not in reviewed_variants:
            raise ValueError(
                f"{reaction_id} no longer matches a reviewed quinone "
                "stoichiometry/bounds/compartment variant"
            )
        reactions.append(reaction)

    model.remove_reactions(reactions, remove_orphans=False)

    branch_metabolites = []
    for metabolite_id in sorted(_REJECTED_QUINONE_ORPHAN_METABOLITES):
        try:
            metabolite = model.metabolites.get_by_id(metabolite_id)
        except KeyError as exc:
            raise ValueError(
                f"Reviewed quinone branch metabolite {metabolite_id} disappeared "
                "before explicit orphan cleanup"
            ) from exc
        if metabolite.reactions:
            connected = sorted(reaction.id for reaction in metabolite.reactions)
            raise ValueError(
                f"Reviewed quinone branch metabolite {metabolite_id} remains "
                f"connected after reaction cleanup: {connected}"
            )
        branch_metabolites.append(metabolite)
    model.remove_metabolites(branch_metabolites, destructive=False)
    logger.info(
        "  Spurious quinone cleanup: removed %s; retained orphan genes for "
        "honest essentiality accounting",
        ", ".join(sorted(expected)),
    )
    return len(reactions)


_COQ9_ROUTE_IDS = (
    "R763",
    "R407",
    "R969",
    "R39",
    "R808",
    "R715",
    "R40",
    "R19",
    "R18",
    "R695",
    "R385",
)


_COQ9_CONNECTED_REACTIONS = {
    *_COQ9_ROUTE_IDS,
    "R1889",
    "R1977",
    "R2062",
    "R262",
    "R305",
    "R570",
    "R573",
    "R740",
}


_COQ9_OPTIONAL_CONNECTED_REACTION = "R2063"


_COQ9_METABOLITES = {
    "m640[C_mi]": (
        "nonaprenyl diphosphate_C45H76O7P2",
        "C30H52O7P2",
        "C45H76O7P2",
        0,
        {"chebi": "CHEBI:53044", "metanetx.chemical": "MNXM1372137"},
    ),
    "m641[C_mi]": (
        "nonaprenyl 4-hydroxybenzoate_C52H78O3",
        "C37H54O3",
        "C52H78O3",
        0,
        {"chebi": "CHEBI:18162", "metanetx.chemical": "MNXM733461"},
    ),
    "m108[C_cy]": (
        "nonaprenyl 4-hydroxybenzoate_C52H78O3",
        "C37H54O3",
        "C52H78O3",
        0,
        {"chebi": "CHEBI:18162", "metanetx.chemical": "MNXM733461"},
    ),
    "m110[C_cy]": (
        "3-nonaprenyl-4,5-dihydroxybenzoate_C52H77O4",
        "C37H53O4",
        "C52H77O4",
        -1,
        {
            "chebi": "CHEBI:62789",
            "metanetx.chemical": "MNXM10069",
            "metacyc.compound": "CPD-9896",
            "seed.compound": "cpd25895",
        },
    ),
    "m939[C_mi]": (
        "3-nonaprenyl-4,5-dihydroxybenzoate_C52H77O4",
        "C37H53O4",
        "C52H77O4",
        -1,
        {
            "chebi": "CHEBI:62789",
            "metanetx.chemical": "MNXM10069",
            "metacyc.compound": "CPD-9896",
            "seed.compound": "cpd25895",
        },
    ),
    "m111[C_mi]": (
        "3-nonaprenyl-4-hydroxy-5-methoxybenzoate_C53H79O4",
        "C38H55O4",
        "C53H79O4",
        -1,
        {
            "chebi": "CHEBI:62791",
            "metanetx.chemical": "MNXM10070",
            "metacyc.compound": "CPD-9898",
            "seed.compound": "cpd25897",
        },
    ),
    "m63[C_mi]": (
        "2-methoxy-6-(all-trans-nonaprenyl)phenol_C52H80O2",
        "C37H56O2",
        "C52H80O2",
        0,
        {
            "chebi": "CHEBI:84522",
            "metanetx.chemical": "MNXM8068",
            "metacyc.compound": "CPD-9866",
            "seed.compound": "cpd25882",
        },
    ),
    "m59[C_mi]": (
        "2-nonaprenyl-6-methoxy-1,4-benzoquinone_C52H78O3",
        "C37H54O3",
        "C52H78O3",
        0,
        {
            "chebi": "CHEBI:203861",
            "metanetx.chemical": "MNXM9872",
            "metacyc.compound": "CPD-11661",
            "seed.compound": "cpd16766",
        },
    ),
    "m61[C_mi]": (
        "2-nonaprenyl-3-methyl-6-methoxy-1,4-benzoquinone_C53H80O3",
        "C38H56O3",
        "C53H80O3",
        0,
        {
            "chebi": "CHEBI:183116",
            "metanetx.chemical": "MNXM9870",
            "metacyc.compound": "CPD-11662",
            "seed.compound": "cpd16764",
        },
    ),
    "m611[C_mi]": (
        "3-demethylubiquinone-9_C53H80O4",
        "C38H56O4",
        "C53H80O4",
        0,
        {
            "chebi": "CHEBI:18238",
            "kegg.compound": "C03226",
            "metanetx.chemical": "MNXM1370748",
        },
    ),
    "m468[C_mi]": (
        "ubiquinone-9_C54H82O4",
        "C39H58O4",
        "C54H82O4",
        0,
        {
            "bigg.metabolite": "q9",
            "chebi": "CHEBI:18160",
            "kegg.compound": "C01967",
            "lipidmaps": "LMPR02010004",
            "metacyc.compound": "UBIQUINONE-9",
            "metanetx.chemical": "MNXM1363635",
            "seed.compound": "cpd01351",
        },
    ),
    "m471[C_mi]": (
        "ubiquinol-9_C54H84O4",
        "C39H60O4",
        "C54H84O4",
        0,
        {
            "bigg.metabolite": "q9h2",
            "chebi": "CHEBI:84424",
            "metacyc.compound": "CPD-9957",
            "metanetx.chemical": "MNXM1094084",
            "seed.compound": "cpd25914",
        },
    ),
}


_COQ9_REACTION_NAMES = {
    "R763": "all-trans-nonaprenyl-diphosphate synthase (four-IPP lump)",
    "R407": "4-hydroxybenzoate nonaprenyltransferase",
    "R969": "nonaprenyl 4-hydroxybenzoate transport",
    "R39": "nonaprenyl 4-hydroxybenzoate hydroxylase",
    "R808": "3-nonaprenyl-4,5-dihydroxybenzoate transport",
    "R715": "SAM:3-nonaprenyl-4,5-dihydroxybenzoate O-methyltransferase",
    "R40": "3-nonaprenyl-4-hydroxy-5-methoxybenzoate decarboxylase",
    "R19": "2-methoxy-6-(all-trans-nonaprenyl)phenol monooxygenase",
    "R18": "2-nonaprenyl-6-methoxy-1,4-benzoquinone methyltransferase",
    "R695": "2-nonaprenyl-3-methyl-6-methoxy-1,4-benzoquinone hydroxylase",
    "R385": "3-demethylubiquinone-9 3-O-methyltransferase",
}


_COQ9_LEGACY_R763 = {
    "m984[C_mi]": -1.0,
    "m985[C_mi]": -1.0,
    "m204[C_mi]": 1.0,
    "m640[C_mi]": 1.0,
}


_COQ9_TARGET_R763 = {
    "m984[C_mi]": -4.0,
    "m985[C_mi]": -1.0,
    "m204[C_mi]": 4.0,
    "m640[C_mi]": 1.0,
}


_COQ9_LEGACY_R385 = {
    "m28[C_mi]": -1.0,
    "m60[C_mi]": -1.0,
    "m611[C_mi]": -1.0,
    "m471[C_mi]": 1.0,
    "m62[C_mi]": 1.0,
}


_COQ9_TARGET_R385 = {
    "m60[C_mi]": -1.0,
    "m611[C_mi]": -1.0,
    "m468[C_mi]": 1.0,
    "m62[C_mi]": 1.0,
}


def _replace_reaction_stoichiometry(model, reaction, target: dict[str, float]) -> None:
    reaction.add_metabolites(
        {metabolite: -coefficient for metabolite, coefficient in reaction.metabolites.items()}
    )
    reaction.add_metabolites(
        {
            model.metabolites.get_by_id(metabolite_id): coefficient
            for metabolite_id, coefficient in target.items()
        }
    )


def _coq_balance(reaction) -> dict[str, float]:
    return {
        key: float(value)
        for key, value in reaction.check_mass_balance().items()
        if not math.isclose(float(value), 0.0, abs_tol=1e-9)
    }


def _coq9_reaction_annotation(reaction, **overrides) -> dict:
    current = reaction.annotation if isinstance(reaction.annotation, dict) else {}
    annotation = {
        "sbo": current.get(
            "sbo",
            "SBO:0000185"
            if len({met.compartment for met in reaction.metabolites}) > 1
            else "SBO:0000176",
        )
    }
    if "ec-code" in current:
        annotation["ec-code"] = copy.deepcopy(current["ec-code"])
    annotation.update(overrides)
    return annotation


def replace_coq6_route_with_coq9(model) -> int:
    """Atomically replace the legacy CoQ6 identities with a balanced CoQ9 route.

    This patch changes chain chemistry only.  It preserves every GPR, bound,
    compartment, biomass coefficient and boundary reaction; no CoQ demand or
    sink is introduced.  The Q9 intermediate identities are homolog-series
    curation, while the oxidized R385 endpoint is an explicit model convention.
    """

    try:
        metabolites = {
            metabolite_id: model.metabolites.get_by_id(metabolite_id)
            for metabolite_id in _COQ9_METABOLITES
        }
        route = {
            reaction_id: model.reactions.get_by_id(reaction_id)
            for reaction_id in _COQ9_ROUTE_IDS
        }
    except KeyError as exc:
        raise ValueError(f"CoQ9 route requires {exc.args[0]}") from exc

    connected = {
        reaction.id
        for metabolite in metabolites.values()
        for reaction in metabolite.reactions
    }
    allowed_connected = (
        _COQ9_CONNECTED_REACTIONS,
        _COQ9_CONNECTED_REACTIONS | {_COQ9_OPTIONAL_CONNECTED_REACTION},
    )
    if connected not in allowed_connected:
        raise ValueError(
            "CoQ9 connected component changed; refusing identity replacement: "
            f"{sorted(connected)}"
        )

    states = set()
    for metabolite_id, (_, legacy_formula, target_formula, charge, _) in _COQ9_METABOLITES.items():
        metabolite = metabolites[metabolite_id]
        pair = (metabolite.formula, metabolite.charge)
        if pair == (legacy_formula, charge):
            states.add("legacy")
        elif pair == (target_formula, charge):
            states.add("coq9")
        else:
            raise ValueError(
                f"{metabolite_id} has unexpected formula/charge {pair!r}; "
                "refusing a partial CoQ9 identity replacement"
            )
    if len(states) != 1:
        raise ValueError("CoQ6/CoQ9 metabolite identities are only partially migrated")
    state = states.pop()

    expected_r763 = _COQ9_LEGACY_R763 if state == "legacy" else _COQ9_TARGET_R763
    expected_r385 = _COQ9_LEGACY_R385 if state == "legacy" else _COQ9_TARGET_R385
    if _stoichiometry_by_metabolite_id(route["R763"]) != expected_r763:
        raise ValueError("R763 no longer matches the reviewed CoQ chain-length signature")
    if _stoichiometry_by_metabolite_id(route["R385"]) != expected_r385:
        raise ValueError("R385 no longer matches the reviewed terminal signature")
    if any(
        not metabolite.formula or metabolite.charge is None
        for reaction in route.values()
        for metabolite in reaction.metabolites
    ):
        raise ValueError("CoQ9 route must be fully formula/charge annotated before replacement")

    impacted_reactions = {
        reaction_id: model.reactions.get_by_id(reaction_id)
        for reaction_id in connected
    }
    counts_before = (len(model.reactions), len(model.metabolites), len(model.genes))
    demands_before = {reaction.id for reaction in model.demands}
    sinks_before = {reaction.id for reaction in model.sinks}
    balance_before = {
        reaction_id: _coq_balance(reaction)
        for reaction_id, reaction in impacted_reactions.items()
    }
    gpr_bounds_before = {
        reaction_id: (reaction.gene_reaction_rule, tuple(reaction.bounds))
        for reaction_id, reaction in impacted_reactions.items()
    }
    metabolite_before = {
        metabolite_id: (
            metabolite.name,
            metabolite.formula,
            metabolite.charge,
            copy.deepcopy(metabolite.annotation),
            copy.deepcopy(metabolite.notes),
        )
        for metabolite_id, metabolite in metabolites.items()
    }
    reaction_before = {
        reaction_id: (
            reaction.name,
            _stoichiometry_by_metabolite_id(reaction),
            copy.deepcopy(reaction.annotation),
            copy.deepcopy(reaction.notes),
        )
        for reaction_id, reaction in impacted_reactions.items()
    }

    try:
        for metabolite_id, (name, _, formula, charge, annotation) in _COQ9_METABOLITES.items():
            metabolite = metabolites[metabolite_id]
            metabolite.name = name
            metabolite.formula = formula
            metabolite.charge = charge
            metabolite.annotation = {"sbo": "SBO:0000247", **annotation}
            notes = dict(metabolite.notes) if isinstance(metabolite.notes, dict) else {}
            notes["curated_coq9_identity"] = (
                "Q9 chain identity supported in Yarrowia; exact intermediate "
                "formula assigned by the homologous +C15H24 series."
            )
            metabolite.notes = notes

        _replace_reaction_stoichiometry(model, route["R763"], _COQ9_TARGET_R763)
        _replace_reaction_stoichiometry(model, route["R385"], _COQ9_TARGET_R385)

        for reaction_id, reaction in route.items():
            reaction.name = _COQ9_REACTION_NAMES[reaction_id]
            overrides = {}
            if reaction_id == "R763":
                overrides = {"ec-code": "2.5.1.85"}
            elif reaction_id == "R407":
                overrides = {"ec-code": "2.5.1.39", "kegg.reaction": "R07273"}
            elif reaction_id == "R385":
                overrides = {
                    "ec-code": "2.1.1.64",
                    "kegg.reaction": "R08781",
                }
            reaction.annotation = _coq9_reaction_annotation(reaction, **overrides)
            notes = dict(reaction.notes) if isinstance(reaction.notes, dict) else {}
            notes.pop("remaining_chain_length_gate", None)
            notes["curated_coq9_chemistry"] = (
                "Legacy C30/CoQ6 identities replaced with the balanced C45/CoQ9 "
                "main chain; GPRs, bounds, compartments and demand are unchanged."
            )
            if reaction_id == "R763":
                notes["coq9_stoichiometry_scope"] = (
                    "The four-IPP lump is balanced model bookkeeping; direct "
                    "Yarrowia evidence supports Q9 chain length, not this exact lump."
                )
            elif reaction_id == "R385":
                notes["terminal_redox_convention"] = (
                    "Balanced oxidized convention: SAM + 3-demethylubiquinone-9 "
                    "-> SAH + ubiquinone-9. Native Yarrowia terminal redox form "
                    "is unresolved; KEGG R08781 is the closest neutral convention."
                )
            reaction.notes = notes

        # These three inherited cross-references explicitly encode Q6.  The
        # reactions remain; only their now-false chain-specific xrefs are removed.
        for reaction_id in {"R262", "R570", "R740", "R2063"} & connected:
            reaction = impacted_reactions[reaction_id]
            reaction.annotation = _coq9_reaction_annotation(reaction)
            notes = dict(reaction.notes) if isinstance(reaction.notes, dict) else {}
            notes["coq9_identity_update"] = (
                "Uses the curated mitochondrial ubiquinone-9/ubiquinol-9 pair; "
                "reaction chemistry, bounds and GPR are otherwise unchanged."
            )
            reaction.notes = notes
        impacted_reactions["R740"].name = "succinate dehydrogenase (ubiquinone-9)"

        for reaction_id in _COQ9_ROUTE_IDS:
            imbalance = _coq_balance(route[reaction_id])
            if imbalance:
                raise ValueError(f"CoQ9 route reaction {reaction_id} is imbalanced: {imbalance}")
        for reaction_id, reaction in impacted_reactions.items():
            if reaction_id in {"R763", "R385"}:
                continue
            after = _coq_balance(reaction)
            if after != balance_before[reaction_id]:
                raise ValueError(
                    f"CoQ9 identity replacement changed {reaction_id} balance: "
                    f"{balance_before[reaction_id]} -> {after}"
                )
        if any(
            (reaction.gene_reaction_rule, tuple(reaction.bounds))
            != gpr_bounds_before[reaction_id]
            for reaction_id, reaction in impacted_reactions.items()
        ):
            raise ValueError("CoQ9 identity replacement changed a GPR or bound")
        if (len(model.reactions), len(model.metabolites), len(model.genes)) != counts_before:
            raise ValueError("CoQ9 identity replacement changed model object counts")
        if {reaction.id for reaction in model.demands} != demands_before:
            raise ValueError("CoQ9 identity replacement changed model demands")
        if {reaction.id for reaction in model.sinks} != sinks_before:
            raise ValueError("CoQ9 identity replacement changed model sinks")

        stale_tokens = ("q6", "u6", "hexaprenyl", "octaprenyl", "ubiquinone-6", "ubiquinol-6")
        for metabolite in metabolites.values():
            if any(token in str(metabolite.annotation).lower() for token in stale_tokens):
                raise ValueError(f"{metabolite.id} retains a CoQ6-specific annotation")
        for reaction in impacted_reactions.values():
            text = f"{reaction.name} {reaction.annotation}".lower()
            if any(token in text for token in stale_tokens):
                raise ValueError(f"{reaction.id} retains a CoQ6-specific name/annotation")
    except Exception:
        for metabolite_id, snapshot in metabolite_before.items():
            metabolite = metabolites[metabolite_id]
            (
                metabolite.name,
                metabolite.formula,
                metabolite.charge,
                metabolite.annotation,
                metabolite.notes,
            ) = snapshot
        for reaction_id, snapshot in reaction_before.items():
            reaction = impacted_reactions[reaction_id]
            reaction.name = snapshot[0]
            _replace_reaction_stoichiometry(model, reaction, snapshot[1])
            reaction.annotation = snapshot[2]
            reaction.notes = snapshot[3]
        raise

    changed_metabolites = sum(
        (
            metabolite.name,
            metabolite.formula,
            metabolite.charge,
            metabolite.annotation,
            metabolite.notes,
        )
        != metabolite_before[metabolite_id]
        for metabolite_id, metabolite in metabolites.items()
    )
    changed_reactions = sum(
        (
            reaction.name,
            _stoichiometry_by_metabolite_id(reaction),
            reaction.annotation,
            reaction.notes,
        )
        != reaction_before[reaction_id]
        for reaction_id, reaction in impacted_reactions.items()
    )
    return changed_metabolites + changed_reactions


_LEGACY_QUINONE_SYNTHOME_GPR = (
    "YALI1F34625g and YALI1B20527g and YALI1A08781g and "
    "YALI1F34675g and YALI1C25352g and YALI1B20835g and YALI1E18269g"
)


_LEGACY_QUINONE_STEP_GPRS = {
    "R715": _LEGACY_QUINONE_SYNTHOME_GPR,
    "R385": _LEGACY_QUINONE_SYNTHOME_GPR,
    "R18": _LEGACY_QUINONE_SYNTHOME_GPR,
    "R695": _LEGACY_QUINONE_SYNTHOME_GPR,
    "R40": "",
    "R19": _LEGACY_QUINONE_SYNTHOME_GPR,
}


_REVIEWED_QUINONE_STEP_GPRS = {
    "R715": "YALI1B20835g",
    "R385": "YALI1B20835g",
    "R18": "YALI1C25352g",
    "R695": "YALI1E18269g",
    "R40": "YALI1F34625g",
    "R19": "",
}


_REVIEWED_QUINONE_STEP_EVIDENCE = {
    "R715": (
        "YALI1B20835g (COQ3 candidate; active UniProt Q6CEG2), CoQ "
        "O-methyltransferase. Cross-species step evidence: E. coli UbiG "
        "structure 4KDC and reconstructed ancestral tetrapod COQ3."
    ),
    "R385": (
        "YALI1B20835g (COQ3 candidate; active UniProt Q6CEG2), CoQ "
        "O-methyltransferase. Cross-species evidence supports the second "
        "CoQ O-methylation; the native Yarrowia substrate redox state remains unresolved."
    ),
    "R18": (
        "YALI1C25352g (COQ5 candidate; active UniProt Q6CBJ6), CoQ-ring "
        "C-methyltransferase. Cross-species step evidence: S. cerevisiae "
        "SAM-bound Coq5 structure 4OBW and reconstructed ancestral COQ5 activity."
    ),
    "R695": (
        "YALI1E18269g (COQ7 candidate; active UniProt Q6C5T9), "
        "demethoxyubiquinone hydroxylase. Cross-species step evidence: human "
        "COQ7:COQ9 structure 7SSS and reconstructed COQ7 activity."
    ),
    "R40": (
        "YALI1F34625g (COQ4; UniProt Q6C074), CoQ-ring C1 decarboxylase/"
        "synthome-organising protein. Cross-species experiments support C1 "
        "decarboxylation, while oxidative versus sequential decarboxylation/"
        "hydroxylation remains unresolved."
    ),
}


def apply_reviewed_quinone_step_gprs(model) -> int:
    """Replace the inherited seven-gene AND with reviewed step-specific GPRs.

    Cross-species biochemical/structural evidence plus compatible AlphaFold
    cores support the COQ3/4/5/7 assignments.  R19 is deliberately left
    GPR-less because the exact COQ6 regioselectivity, product redox state and
    electron-transfer partners do not yet match the model reaction.  An empty
    R19 rule means unknown, not spontaneous.  No demand, bound or chemistry is
    changed.
    """

    try:
        reactions = {
            reaction_id: model.reactions.get_by_id(reaction_id)
            for reaction_id in _REVIEWED_QUINONE_STEP_GPRS
        }
        model.metabolites.get_by_id("m468[C_mi]")
        for gene_id in set(_REVIEWED_QUINONE_STEP_GPRS.values()) - {""}:
            model.genes.get_by_id(gene_id)
    except KeyError as exc:
        raise ValueError(f"Reviewed quinone GPR patch requires {exc.args[0]}") from exc

    if model.metabolites.get_by_id("m468[C_mi]").formula != "C54H82O4":
        raise ValueError("Reviewed quinone GPR patch requires the formal CoQ9 chemistry")
    imbalanced = [
        reaction_id
        for reaction_id, reaction in reactions.items()
        if _coq_balance(reaction)
    ]
    if imbalanced:
        raise ValueError(
            "Reviewed quinone GPR targets are not mass/charge balanced: "
            f"{sorted(imbalanced)}"
        )

    actual = {
        reaction_id: reaction.gene_reaction_rule.strip()
        for reaction_id, reaction in reactions.items()
    }
    if actual == _REVIEWED_QUINONE_STEP_GPRS:
        changed = 0
    elif actual == _LEGACY_QUINONE_STEP_GPRS:
        changed = len(_REVIEWED_QUINONE_STEP_GPRS)
    else:
        raise ValueError(
            "Quinone step GPRs are partially migrated or unexpected; refusing "
            f"an atomic replacement: {actual}"
        )

    for reaction_id, target_rule in _REVIEWED_QUINONE_STEP_GPRS.items():
        reaction = reactions[reaction_id]
        reaction.gene_reaction_rule = target_rule
        notes = dict(reaction.notes) if isinstance(reaction.notes, dict) else {}
        if reaction_id == "R19":
            notes["curated_gpr_correction"] = (
                "Removed the unsupported inherited seven-gene synthome AND. "
                "YALI1A08781g (COQ6 candidate; active UniProt F2Z6J4) remains "
                "deferred for this exact reaction; an empty GPR denotes unknown "
                "catalyst identity, not a spontaneous reaction."
            )
            notes["gpr_evidence_status"] = "deferred_reaction_identity_unresolved"
            notes["gpr_evidence_limit"] = (
                "COQ6 family/AlphaFold compatibility is supported, but R19 "
                "regioselectivity, product redox state and ferredoxin/reductase "
                "electron transfer remain unresolved."
            )
        else:
            notes["curated_gpr_correction"] = _REVIEWED_QUINONE_STEP_EVIDENCE[
                reaction_id
            ]
            notes["gpr_evidence_status"] = (
                "cross_species_experimental_plus_alphafold_compatible; "
                "native_yarrowia_biochemistry_unverified"
            )
        reaction.notes = notes

    return changed



def run_quinone_step(model, operation):
    """Preflight legacy mutations on a copy so a local conflict preserves input.

    The saved branch cleanup writes annotations before all its preconditions;
    preflight avoids partial writes without a second implementation of its rules.
    """
    try:
        operation(model.copy())
    except (ValueError, KeyError) as error:
        logger.warning("Existing Q9 step %s: %s", operation.__name__, error)
        return {"item": operation.__name__, "status": "conflict", "reason": str(error)}
    changed = operation(model)
    return {"item": operation.__name__, "status": "applied" if changed else "already_correct",
            "changed_objects": changed}
