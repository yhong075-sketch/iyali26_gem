"""
main.py — orchestration entry point for the iYali26 annotation pipeline.
"""

import csv
import hashlib
import json
import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

from cobra.io import read_sbml_model

from .coq9 import CURATION_PATH, GENE_EVIDENCE_PATH, apply_coq9_curation, gene_evidence
from .biomass import fix_biomass_reaction
from .config import (
    CACHE_DIR,
    REPO_ROOT,
    PROJECT_PATHS,
    CURATION_DATA_DIR,
    MNX_DIR,
    OUTPUT_MODEL_PATH,
    STARTING_MODEL_PATH,
)
from .exchange import configure_medium, set_exchange_bounds
from .gap_fill_direction import DEFAULT_GAP_FILL_DIRECTION_TABLE
from .gaps import DUPLICATE_PAIRS, add_gap_fill_reactions, find_gaps, merge_duplicate_metabolites, report_gaps
from .annotate_reactions_extended import annotate_remaining_reactions
from .ec_annotation import enrich_genes_with_ec
from .essentiality_evidence import sha256_file
from .genes import annotate_genes, apply_curated_gene_annotation_overrides
from .idmapping import _enrich_via_idmapping
from .io import load_chem_prop, load_chem_xref, load_mnxm_depr, load_reac_prop, load_reac_xref
from .metabolites import annotate_metabolites, fix_proton_water_balance, normalize_all_annotations
from .microspecies import apply_curated_microspecies, normalize_hydroxide_reactions
from .patches import add_direct_enzyme_like_gprs, add_isozyme_gprs, add_r612_ura3_gpr, annotate_isozyme_genes, apply_all_patches, apply_curated_essentiality_patches, apply_reviewed_quinone_step_gprs, clean_ec_overload, correct_external_ndh2_gpr_and_remove_duplicate, extend_acyl_pool_c161, fill_neutral_formulas, fix_activex_names, fix_charge_stage1, fix_charge_stage2, fix_ec_code_format, move_tcdb_out_of_ec, remove_misannotated_gprs, remove_spurious_quinone_branches, remove_spurious_transport_reactions, remove_stale_adp_atp_transporter_ec_codes, replace_coq6_route_with_coq9, split_trna_charging_from_biomass
from .provisional_capacity import apply_provisional_isozyme_capacities
from .reactions import annotate_reactions, backfill_reaction_xrefs
from .sbml import write_deterministic_sbml_model
from .execution import guarded_execution
from .r608 import apply_r608_candidate
from .reaction_selection import SELECTION_PATH, apply_metadata_reaction_selection
from .r1159_direction import apply_r1159_direction

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def _deterministic_model_sha256(model) -> str:
    """Hash the exact in-memory canonical-stage model before an overlay."""
    with TemporaryDirectory(prefix="iyali26-provisional-reference-") as directory:
        reference_path = Path(directory) / "model.xml"
        write_deterministic_sbml_model(model, reference_path)
        return sha256_file(reference_path)


@guarded_execution
def build_reference_chain(
    provisional_capacity_path: Path | None = None,
    trna_biomass_mode: str | None = None,
    output_model_path: Path = OUTPUT_MODEL_PATH,
    allow_network: bool = True,
    canonical_copy: bool = False,
    no_solve: bool = False,
    r608_curation_path: Path | None = None,
    starting_model_path: Path = STARTING_MODEL_PATH,
    coq9_mode: str = "metadata",
):
    starting_model_path = Path(starting_model_path)
    if coq9_mode not in {"off", "metadata", "qcycle"}:
        raise ValueError(f"Unknown CoQ9 mode: {coq9_mode}")
    if starting_model_path.resolve() == Path(output_model_path).resolve():
        raise ValueError("Input and output must be separate files")
    if r608_curation_path is not None:
        no_solve, allow_network = True, False
    output_model_path = Path(output_model_path)
    if canonical_copy and (provisional_capacity_path is not None or trna_biomass_mode is not None):
        raise ValueError("a canonical copy cannot include experimental overlays")
    if canonical_copy and (
        output_model_path.exists()
        or output_model_path.resolve() in {OUTPUT_MODEL_PATH.resolve(), starting_model_path.resolve()}
    ):
        raise FileExistsError("a canonical copy requires a new, separate output path")
    if provisional_capacity_path is not None:
        provisional_capacity_path = Path(provisional_capacity_path)
        if output_model_path.resolve() == OUTPUT_MODEL_PATH.resolve():
            raise ValueError(
                "a provisional capacity profile cannot overwrite canonical "
                "model.xml; provide a separate output_model_path"
            )
    if trna_biomass_mode not in {None, "split"}:
        raise ValueError(f"unsupported tRNA biomass mode: {trna_biomass_mode}")
    if (
        trna_biomass_mode is not None
        and output_model_path.resolve() == OUTPUT_MODEL_PATH.resolve()
    ):
        raise ValueError(
            "an experimental tRNA biomass overlay cannot overwrite canonical "
            "model.xml; provide a separate output_model_path"
        )

    if no_solve or r608_curation_path is not None:
        if output_model_path.exists() or output_model_path.is_symlink() or output_model_path.resolve() in {OUTPUT_MODEL_PATH.resolve(), starting_model_path.resolve()}:
            raise FileExistsError("no-solve and R608 builds require a new, separate output path")
    if r608_curation_path is not None and (provisional_capacity_path is not None or trna_biomass_mode is not None):
        raise ValueError("R608 cannot be combined with experimental overlays")

    if not starting_model_path.exists():
        raise FileNotFoundError(f"Could not find starting model at {starting_model_path}")

    logger.info(f"Loading raw model: {starting_model_path.name}")
    model = read_sbml_model(str(starting_model_path))

    mnx_ok = MNX_DIR.exists() and (MNX_DIR / "chem_xref.tsv").exists()

    if mnx_ok:
        # Load MetaNetX tables once, reuse across functions
        chem_xref = load_chem_xref(MNX_DIR / "chem_xref.tsv")
        chem_prop_data = load_chem_prop(MNX_DIR / "chem_prop.tsv")
        reac_xref = load_reac_xref(MNX_DIR / "reac_xref.tsv")

        reac_prop_path = MNX_DIR / "reac_prop.tsv"
        reac_prop = load_reac_prop(reac_prop_path) if reac_prop_path.exists() else None
        if reac_prop is None:
            logger.warning("reac_prop.tsv not found — Strategy C (fingerprint) disabled")

        mnxm_depr = load_mnxm_depr(MNX_DIR / "chem_depr.tsv")
        if not mnxm_depr:
            logger.info(
                "  chem_depr.tsv not found in data/metanetx/ — "
                "download from https://www.metanetx.org/mnxdoc/mnxref.html "
                "to improve fingerprint match rates for deprecated MNXM IDs"
            )

        # Pre-annotation cleanup: strip Excel "ActiveX VT_ERROR" corruption from
        # metabolite names so the annotation step can match them to MetaNetX.
        logger.info("=== Pre-annotation: clean ActiveX-corrupted names ===")
        n_activex = fix_activex_names(model)
        logger.info(f"  ActiveX name cleanup: {n_activex} metabolite name(s) cleaned")

        # Priority 1 + 2a
        logger.info("=== Priority 1+2a: metabolite annotation + formulas ===")
        annotate_metabolites(model, chem_xref, chem_prop_data)

        # Patches: fix known data bugs (NADP+ formula, ceramide formulas)
        # Must run before H+/H2O balance so balanced reactions see correct formulas.
        logger.info("=== Patches: known data-bug fixes ===")
        apply_all_patches(model)

    else:
        logger.warning(
            f"MetaNetX files not found in {MNX_DIR}. "
            "Download chem_xref.tsv, chem_prop.tsv, reac_xref.tsv from "
            "https://www.metanetx.org/mnxdoc/mnxref.html"
        )

    # Curated chemical identities override raw/name-derived and MetaNetX values.
    # Formula and charge are applied atomically from the pH 7.3 Rhea/ChEBI table.
    # Families that still need a connected-component migration remain recorded as
    # component_review and are never changed by this pipeline stage.
    logger.info("=== Chemical convention: Rhea/ChEBI pH 7.3 microspecies ===")
    microspecies_report = apply_curated_microspecies(model)
    logger.info(
        "  Microspecies: %d metabolite(s) changed; %d already canonical; "
        "%d family/families deferred",
        microspecies_report["changed_metabolites"],
        microspecies_report["already_canonical"],
        microspecies_report["deferred_families"],
    )

    # MetaNetX canonicalisation: OH- is represented as H2O - H+.  Only exact,
    # single-compartment hydroxide reactions are normalised; transport and
    # boundary cases require explicit curation of ion coupling.
    hydroxide_report = normalize_hydroxide_reactions(model)
    logger.info(
        "  Hydroxide normalisation: %d reaction(s) changed; %d rejected",
        hydroxide_report["changed_reactions"],
        hydroxide_report["rejected_reactions"],
    )

    if mnx_ok:
        # Priority 2b.  The balancer is charge-aware: it may add H+/H2O only
        # when mass and charge can be brought to zero simultaneously.
        logger.info("=== Priority 2b: charge-aware H+/H2O balance ===")
        fix_proton_water_balance(model)

        # Priority 4a  (Strategy C needs metabolite annotations from 1+2a above)
        logger.info("=== Priority 4a: reaction annotation ===")
        name_exclusions = ({"R305": "ferrocytochrome-c:oxygen oxidoreductase"}
                           if coq9_mode != "off" else {})
        annotate_reactions(model, reac_xref, reac_prop, name_exclusions=name_exclusions)

    # Stoichiometric consistency: merge known duplicate metabolite pairs
    # These pairs were identified by MIS analysis (metabolites appearing in
    # stoichiometrically inconsistent reactions with identical formulas/charges
    # but separate model IDs).  Merge before any FBA/FVA to make the network solvable.
    logger.info("=== Stoichiometric consistency: merging duplicate metabolites ===")
    merge_duplicate_metabolites(model, DUPLICATE_PAIRS)

    # Priority 2c: exchange-bound calibration (must run before any FBA/FVA)
    logger.info("=== Priority 2c: exchange bounds (minimal medium) ===")
    set_exchange_bounds(model)   # uses MINIMAL_MEDIUM_BIGG (Tier 1) + MINIMAL_MEDIUM_NAMES (Tier 2)

    logger.info("=== Priority 2c+: mineral salts + vitamins medium extension ===")
    configure_medium(model)      # adds Mg, K, Na, biotin, thiamine, pyridoxine

    # Priority 3 (always run — independent of MetaNetX)
    logger.info("=== Priority 3: biomass reaction R1372 ===")
    fix_biomass_reaction(model, no_solve=no_solve)

    # Priority 4b — network required, skip if offline
    logger.info("=== Priority 4b: gene annotation via UniProt ===")
    annotate_genes(model, allow_network=allow_network)

    # Curated assembly-identity corrections run after broad automatic matching
    # but before ID-mapping and EC enrichment, so stale accessions cannot seed
    # downstream annotations. These edits affect metadata only, never GPRs.
    logger.info("=== Priority 4b+: curated gene identity overrides ===")
    n_gene_identity_overrides = apply_curated_gene_annotation_overrides(model)
    logger.info(
        "  Curated gene identity overrides: %d gene(s) corrected",
        n_gene_identity_overrides,
    )

    # Priority 4c — ncbigene → UniProt ID-mapping for genes still missing uniprot
    logger.info("=== Priority 4c: UniProt ID-mapping (ncbigene → UniProtKB) ===")
    _enrich_via_idmapping(model, allow_network=allow_network)

    # Priority 4d — enrich genes with EC numbers via UniProt stream API
    logger.info("=== Priority 4d: gene EC number enrichment via UniProt ===")
    enrich_genes_with_ec(model, allow_network=allow_network)

    # Priority 4e — extended reaction annotation (exchange / transport / EC→MNXR)
    # Runs after 4d so gene EC numbers are already populated.
    if mnx_ok:
        logger.info("=== Priority 4e: extended reaction annotation ===")
        annotate_remaining_reactions(model, reac_xref, reac_prop, mnxm_depr=mnxm_depr,
                                     name_exclusions=name_exclusions)

    # === EC backfill: copy gene EC numbers to reaction annotations ===
    logger.info("=== EC backfill: gene ec-code → reaction annotation ===")
    ec_backfill_count = 0
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        if "ec-code" in ann:
            continue
        ec_set = set()
        for gene in rxn.genes:
            g_ann = gene.annotation if isinstance(gene.annotation, dict) else {}
            ec_raw = g_ann.get("ec-code", [])
            if isinstance(ec_raw, str):
                ec_raw = [ec_raw]
            for ec in ec_raw:
                ec = ec.strip()
                if ec:
                    ec_set.add(ec)
        if ec_set:
            if not isinstance(rxn.annotation, dict):
                rxn.annotation = {}
            rxn.annotation["ec-code"] = sorted(ec_set)
            ec_backfill_count += 1
    logger.info(f"  EC backfill: ec-code added to {ec_backfill_count} reactions")

    # Safety net: fill missing cross-refs for all reactions that already have MNXR
    if mnx_ok:
        logger.info("=== Reaction xref backfill: fill missing bigg/kegg/rhea/ec-code from MNXR ===")
        backfill_reaction_xrefs(model, reac_xref, reac_prop)

    # Priority 2b (second pass): re-run H+/H2O balance now that more metabolites
    # may have formulas (via reaction annotation pulling in new MNXM → chem_prop).
    # Catches reactions skipped in the first pass for missing formula.
    logger.info("=== Priority 2b (second pass): H+/H2O balance after annotation ===")
    fix_proton_water_balance(model)

    # Replace the inherited C30/CoQ6 identities before FVA, so every downstream
    # diagnostic uses the balanced C45/CoQ9 main chain.  This identity patch adds
    # no sink, demand, biomass coefficient, bound or GPR.
    logger.info("=== Curated quinone chemistry: balanced CoQ9 main chain ===")
    n_coq9_changes = replace_coq6_route_with_coq9(model)
    logger.info("  CoQ9 identity/chemistry changes: %d object(s)", n_coq9_changes)

    if not no_solve:
        # Priority 5: gap analysis — FVA before gap-fill
        logger.info("=== Priority 5: gap analysis (FVA, post-medium) ===")
        gaps = find_gaps(model)
        report_gaps(gaps)
        blocked_before_medium = len(gaps["blocked_reactions"])
        logger.info(f"  Blocked reactions after medium extension: {blocked_before_medium}")

    else:
        logger.info("Gap FVA diagnostics not run (no-solve mode)")

    # Priority 6: gap-fill — insert P0 reactions from gap_fill_prioritized.csv
    gap_fill_csv = CURATION_DATA_DIR / "gap_fill_prioritized.csv"
    if gap_fill_csv.exists():
        logger.info("=== Priority 6: gap-fill reaction insertion (P0) ===")
        add_gap_fill_reactions(
            model,
            csv_path=gap_fill_csv,
            mnx_dir=MNX_DIR if mnx_ok else None,
            cache_dir=CACHE_DIR,
            direction_curation_path=DEFAULT_GAP_FILL_DIRECTION_TABLE,
        )
        if not no_solve:
            logger.info("=== Priority 6b: post-gap-fill FVA ===")
            gaps_after = find_gaps(model)
            before = len(gaps["blocked_reactions"])
            after  = len(gaps_after["blocked_reactions"])
            logger.info(
                f"  Blocked reactions: {before} → {after}  "
                f"(Δ {before - after:+d} unblocked)"
            )
            report_gaps(gaps_after)
    else:
        logger.warning(f"gap_fill_prioritized.csv not found at {gap_fill_csv} — skipping")

    # Migrate c_va (lowercase typo) → C_va before compartment naming
    cva_migrated = 0
    for met in model.metabolites:
        if met.compartment == "c_va":
            met.compartment = "C_va"
            cva_migrated += 1
    if cva_migrated:
        logger.info(f"Compartment fix: migrated {cva_migrated} metabolite(s) from c_va → C_va")

    # Fix compartment names so COBRApy can identify the external compartment
    _COMPARTMENT_NAMES = {
        "C_cy": "cytoplasm",
        "C_ex": "extracellular",
        "C_mi": "mitochondria",
        "C_nu": "nucleus",
        "C_er": "endoplasmic reticulum",
        "C_go": "Golgi apparatus",
        "C_va": "vacuole",
        "C_lp": "lipid particle",
        "C_pe": "peroxisome",
        "C_em": "endosomal membrane",
        "C_en": "endosome",
        "C_gm": "Golgi membrane",
        "C_mm": "mitochondrial membrane",
        "C_vm": "vacuolar membrane",
    }
    existing_comps = set(met.compartment for met in model.metabolites)
    model.compartments = {c: _COMPARTMENT_NAMES.get(c, c) for c in existing_comps}
    logger.info(f"Compartments set: {model.compartments}")

    # Verify COBRApy can resolve the external compartment via the "extracellular" name.
    # COBRApy checks compartment names for "extra" substring to identify the boundary.
    _ex_name = model.compartments.get("C_ex", "")
    if _ex_name != "extracellular":
        logger.warning(
            f"C_ex compartment name is {_ex_name!r}, expected 'extracellular' — "
            "COBRApy may not auto-detect the external compartment"
        )
    else:
        logger.info("  C_ex compartment verified as 'extracellular'")

    # === SBO annotation (full coverage) ===
    import re as _re
    logger.info("=== SBO annotation: all model objects ===")

    _BIOMASS_RE     = _re.compile(r"BIOMASS|biomass|newBiom|R1372")
    _MAINTENANCE_RE = _re.compile(r"MAINTENANCE|ATPM")

    def _set_sbo(obj, term: str) -> bool:
        ann = obj.annotation if isinstance(obj.annotation, dict) else {}
        if "sbo" in ann:
            return False
        if not isinstance(obj.annotation, dict):
            obj.annotation = {}
        obj.annotation["sbo"] = term
        return True

    met_set = sum(_set_sbo(m, "SBO:0000247") for m in model.metabolites)
    logger.info(f"  Metabolites : set={met_set}  already_had_sbo={len(model.metabolites) - met_set}")

    gene_set = sum(_set_sbo(g, "SBO:0000243") for g in model.genes)
    logger.info(f"  Genes       : set={gene_set}  already_had_sbo={len(model.genes) - gene_set}")

    _exchanges = set(model.exchanges)
    _demands   = set(model.demands)
    _sinks     = set(model.sinks)
    sbo_counts = {k: 0 for k in ("exchange", "demand", "sink", "biomass", "maintenance", "pool", "transport", "biochemical", "skip")}

    for rxn in model.reactions:
        rid = rxn.id
        if rxn in _exchanges:
            term, kind = "SBO:0000627", "exchange"
        elif rxn in _demands:
            term, kind = "SBO:0000628", "demand"
        elif rxn in _sinks:
            term, kind = "SBO:0000632", "sink"
        elif _BIOMASS_RE.search(rid) or rid.startswith(("xBIOMASS", "newBiom", "biomass_C")):
            term, kind = "SBO:0000629", "biomass"
        elif _MAINTENANCE_RE.search(rid) or rid.startswith("xMAINTENANCE"):
            term, kind = "SBO:0000630", "maintenance"
        elif rid.startswith(("xLIPID", "xAMINOACID", "xPOOL_")):
            term, kind = "SBO:0000395", "pool"
        elif len({m.compartment for m in rxn.metabolites}) >= 2:
            term, kind = "SBO:0000185", "transport"
        else:
            term, kind = "SBO:0000176", "biochemical"

        if _set_sbo(rxn, term):
            sbo_counts[kind] += 1
        else:
            sbo_counts["skip"] += 1

    logger.info(
        f"  Reactions   : "
        f"exchange={sbo_counts['exchange']}  demand={sbo_counts['demand']}  "
        f"sink={sbo_counts['sink']}  biomass={sbo_counts['biomass']}  "
        f"maintenance={sbo_counts['maintenance']}  pool={sbo_counts['pool']}  "
        f"transport={sbo_counts['transport']}  biochemical={sbo_counts['biochemical']}  "
        f"already_had_sbo={sbo_counts['skip']}"
    )

    logger.info("=== Final pass: normalise all annotation keys/values ===")
    normalize_all_annotations(model)

    # EC-code format compliance: pad three-segment EC codes (e.g. "3.1.3" →
    # "3.1.3.-") for identifiers.org validity.  Must run after all EC codes
    # are populated (gene enrichment, EC backfill, xref backfill) and after
    # normalisation, so it is invoked here rather than in apply_all_patches.
    logger.info("=== EC-code format compliance ===")
    n_tcdb = move_tcdb_out_of_ec(model)
    logger.info(f"  TCDB cleanup: moved {n_tcdb} TCDB number(s) from ec-code → tcdb")
    n_ec_padded = fix_ec_code_format(model)
    logger.info(f"  EC format: padded {n_ec_padded} three-segment EC code(s) with '.-'")

    # EC-overload cleanup must run LAST: after all EC back-fill and format
    # steps, so the back-fill cannot re-introduce the dropped EC numbers.
    n_ec_cleaned = clean_ec_overload(model)
    logger.info(f"  EC overload cleanup: removed polluting EC from {n_ec_cleaned} reaction(s)")

    # R815/R816 are transporters; their inherited EC 2.7.4.6 is a kinase code.
    # This must follow all generic EC backfills.
    n_adp_atp_ec_cleaned = remove_stale_adp_atp_transporter_ec_codes(model)
    logger.info(
        "  ADP/ATP transporter EC cleanup: corrected %d reaction(s)",
        n_adp_atp_ec_cleaned,
    )

    # Isozyme GPR additions: add curated CLIB89 isozyme genes to existing
    # reactions' GPR via 'or' (safe subset only; see patches.add_isozyme_gprs).
    n_gpr_added = add_isozyme_gprs(model)
    logger.info(f"  Isozyme GPR additions: {n_gpr_added} (reaction, gene) pair(s) added")

    # Remove mis-annotated genes from GPRs (E07744g trehalase wrongly in R765/R766
    # transketolase; E11370g GatB wrongly in R671 prephenate DH). Safe 'or' cases,
    # true partner retained -> 0 growth impact. See patches.remove_misannotated_gprs.
    n_gpr_removed = remove_misannotated_gprs(model)
    logger.info(f"  Mis-annotation GPR removals: {n_gpr_removed} (reaction, gene) pair(s) removed")

    # Directly reviewed GPR corrections are deliberately separate from the
    # automated isozyme expansion: R612 receives the verified URA3 gene, while
    # R570 is corrected to external NDH2 and duplicate R2063 is removed. R1889
    # remains GPR-less pending a complex-I subunit requiredness review.
    n_r612_gpr = add_r612_ura3_gpr(model)
    logger.info("  Reviewed R612 URA3 GPR correction: %d reaction(s) changed", n_r612_gpr)
    n_ndh2_changes = correct_external_ndh2_gpr_and_remove_duplicate(model)
    logger.info("  Reviewed external-NDH2 correction: %d change(s)", n_ndh2_changes)

    # Three direct enzyme-like assignments: LIP2 on extracellular triolein
    # hydrolysis and the existing EC 2.7.1.35 gene on the two remaining B6
    # vitamer kinase reactions.  See patches.add_direct_enzyme_like_gprs.
    n_direct_enzyme_like_gprs = add_direct_enzyme_like_gprs(model)
    logger.info(
        "  Direct enzyme-like GPR corrections: %d reaction(s) changed",
        n_direct_enzyme_like_gprs,
    )

    # Remove two dead, mis-annotated quinone alternatives while preserving the
    # mitochondrial R407/COQ2 route. This is a topology/identity cleanup, not
    # an essentiality-recall patch.
    n_quinone_removed = remove_spurious_quinone_branches(model)
    logger.info(
        "  Spurious quinone branch cleanup: %d reaction(s) removed",
        n_quinone_removed,
    )
    n_quinone_gprs = apply_reviewed_quinone_step_gprs(model)
    logger.info(
        "  Reviewed step-specific quinone GPR corrections: %d reaction(s) changed",
        n_quinone_gprs,
    )

    # R1172 is retained by the user-requested policy; the historical removal
    # rationale is recorded in the versioned retention data and output notes.
    n_transport_removed = remove_spurious_transport_reactions(model)
    logger.info(f"  Spurious transport removals: {n_transport_removed} reaction(s) removed")

    # Evidence-gated essentiality corrections. Schema-v2 rows require a direct
    # Y. lipolytica evidence dossier, skeptic pass, explicit human approval and
    # a live target fingerprint. EG-GPR-001 is the only legacy-v1 exception.
    essentiality_patches = apply_curated_essentiality_patches(model)
    logger.info(f"  Curated essentiality patches: {len(essentiality_patches)} applied")

    # Annotate those newly added genes (they entered after the main gene
    # annotation + SBO steps, so they need sbo / ncbigene / kegg / uniprot here).
    n_gene_annot = annotate_isozyme_genes(model, network=allow_network)
    logger.info(f"  Isozyme gene annotation: {n_gene_annot} gene(s) annotated")

    # Fill formulas for definite-neutral metabolites (charge=0, unambiguous).
    n_form = fill_neutral_formulas(model)
    logger.info(f"  Neutral formula fill: {n_form} metabolite copy(ies) filled")

    # Lipid chain-menu extension: add C16:1 palmitoleoyl-CoA to the acyl-CoA pools
    # (Y. lipolytica makes ~8% but the pool omitted it). Idempotent.
    n_c161 = extend_acyl_pool_c161(model)
    logger.info(f"  C16:1 acyl-CoA pool extension: {n_c161} pool(s) extended")

    # Charge fix Stage 1: 4 free anionic metabolites left at charge 0 (formula
    # already deprotonated) -> -1. Fixes 5 reactions, breaks 0 (verified). Idempotent.
    n_chg1 = fix_charge_stage1(model)
    logger.info(f"  Charge Stage 1: {n_chg1} metabolite(s) set to -1")

    # Charge fix Stage 2: 7 anion-stored metabolites (InChI dH < 0) -> their dH charge,
    # formula unchanged. Fixes 7 reactions, breaks 0 (verified). Idempotent.
    n_chg2 = fix_charge_stage2(model)
    logger.info(f"  Charge Stage 2: {n_chg2} metabolite(s) set to InChI-dH charge")

    # Final chemistry gate.  Formula fills and curated charge patches above run
    # after the historical second balance pass, so validate the exact common-ion
    # target set again and give newly checkable reactions one strict balance pass.
    # This currently repairs R2065, R2097 and R2274; the simultaneous mass/charge
    # condition prevents an H+/H2O adjustment from hiding an ion-form mismatch.
    logger.info("=== Final chemistry gate: microspecies + charge-aware balance ===")
    final_microspecies_report = apply_curated_microspecies(model)
    logger.info(
        "  Final microspecies validation: %d changed; target set %s",
        final_microspecies_report["changed_metabolites"],
        final_microspecies_report["target_set_fingerprint"],
    )
    final_hydroxide_report = normalize_hydroxide_reactions(model)
    logger.info(
        "  Final hydroxide normalisation: %d changed; %d rejected",
        final_hydroxide_report["changed_reactions"],
        final_hydroxide_report["rejected_reactions"],
    )
    final_balance_report = fix_proton_water_balance(model)
    logger.info(
        "  Final H+/H2O balance: %d changed; %d rejected",
        final_balance_report["changed_reactions"],
        final_balance_report["rejected_reactions"],
    )
    if final_balance_report["changes"]:
        logger.info(
            "  Final balanced reactions: %s",
            ", ".join(
                change["reaction_id"] for change in final_balance_report["changes"]
            ),
        )

    # Optional experimental capacity profile. Apply it after every canonical
    # structural and chemistry operation so later pipeline steps cannot make a
    # copied backup branch diverge from its parent. This remains separate from
    # the evidence-gated patch table and can never target canonical model.xml.
    if provisional_capacity_path is not None:
        provisional_reference_sha256 = _deterministic_model_sha256(model)
        capacity_audit = apply_provisional_isozyme_capacities(
            model,
            provisional_capacity_path,
            reference_model_sha256=provisional_reference_sha256,
        )
        logger.warning(
            "  Applied %d provisional isozyme-capacity split(s); experimental "
            "profile only",
            len(capacity_audit),
        )

    # The experimental biomass requirement remains coupled through all 20 tRNAs.
    # Final metadata selection checks and preserves this reference representation.
    canonical_build = canonical_copy or output_model_path.resolve() == OUTPUT_MODEL_PATH.resolve()
    if canonical_build or trna_biomass_mode == "split":
        trna_audit = split_trna_charging_from_biomass(model)
        if canonical_build:
            model.reactions.get_by_id("biomass_C").notes[
                "canonical_trna_biomass_mode"
            ] = "split_v1"
            logger.info(
                "  Applied %d fully split tRNA-biomass coupling reactions; "
                "canonical B-group representation",
                len(trna_audit),
            )
        else:
            logger.warning(
                "  Applied %d fully split tRNA-biomass coupling reactions; "
                "experimental B-group overlay only",
                len(trna_audit),
            )

    # R608 is opt-in and follows all annotation, chemistry and tRNA steps.
    if r608_curation_path is not None:
        import json
        model, _ = apply_r608_candidate(
            model, json.loads(Path(r608_curation_path).read_text()),
            source_model_sha256=_deterministic_model_sha256(model),
        )

    selection = apply_metadata_reaction_selection(model)
    if not selection["complete"]:
        logger.warning("Metadata reaction selection has local conflicts: inspect the build record")
    retention = annotate_retained_reactions(model)
    coq9 = apply_coq9_curation(model, coq9_mode)
    if not coq9["requested_mode_complete"]:
        logger.warning("CoQ9 local conflicts: inspect the build record")
    selection["post_selection_r1159_direction"] = apply_r1159_direction(model)
    logger.info("R1159 direction: %s", selection["post_selection_r1159_direction"])
    output_model_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Saving updated model to: {output_model_path.name}")
    # COBRApy stores Group.members as sets, so its stock writer emits pathway
    # members in hash-random order.  Canonicalise that serialization detail so
    # the byte-level SHA used by the evidence ledger is reproducible.
    write_deterministic_sbml_model(model, output_model_path)
    logger.info("Model build complete.")
    return model, coq9, retention, selection


RETENTION_PATH = REPO_ROOT / "data" / "reference_build" / "retained_reactions.json"


def annotate_retained_reactions(model):
    """Record explicit retention without inventing a carrier or changing GPRs."""
    policy = json.loads(RETENTION_PATH.read_text())
    records = []
    for rid, row in policy["reactions"].items():
        reaction = model.reactions.get_by_id(rid)
        target = {"retention_previous_policy": row["previous_policy"],
                  "retention_current_status": row["current_status"],
                  "retention_authorization": policy["authorization"]}
        changed = any(reaction.notes.get(k) != v for k, v in target.items())
        reaction.notes.update(target)
        records.append({"item": rid, "status": "applied" if changed else "already_correct",
                        **row, "stoichiometry": {m.id: c for m, c in reaction.metabolites.items()},
                        "bounds": list(reaction.bounds), "gpr": reaction.gene_reaction_rule})
    return records


def build_model(args):
    """Shared CLI/API route for the reference chain and explicit reaction selection."""
    def sha(path):
        with Path(path).open("rb") as handle:
            return hashlib.file_digest(handle, "sha256").hexdigest()
    def git(*arguments):
        result = subprocess.run(["git", "-C", str(REPO_ROOT), *arguments], capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None
    started = datetime.now(timezone.utc).isoformat()
    input_sha = sha(args.starting_model)
    source_sha = {str(p.relative_to(REPO_ROOT)): sha(p) for p in sorted(Path(__file__).parent.glob("*.py"))}
    data_paths = [args.starting_model, CURATION_PATH, GENE_EVIDENCE_PATH, RETENTION_PATH, SELECTION_PATH]
    data_paths += [p for p in (REPO_ROOT / "data" / "reference_build").rglob("*") if p.is_file()]
    for folder in (MNX_DIR, CACHE_DIR, PROJECT_PATHS.locus_map,
                   PROJECT_PATHS.research_root / "reference" / "ncbi",
                   PROJECT_PATHS.research_root / "reference" / "kegg"):
        if folder.is_dir():
            data_paths += [p for p in folder.iterdir() if p.is_file() and not p.name.startswith(".")]
    data_sha = {str(p.resolve()): sha(p) for p in sorted(set(data_paths))}
    ordinary = args.provisional_capacity_profile is None and args.trna_biomass_mode is None
    model, coq9, retention, selection = build_reference_chain(
        provisional_capacity_path=args.provisional_capacity_profile,
        trna_biomass_mode=args.trna_biomass_mode,
        output_model_path=args.output_model,
        allow_network=not args.offline,
        canonical_copy=args.canonical_copy or ordinary,
        no_solve=args.no_solve,
        r608_curation_path=args.r608_curation,
        starting_model_path=args.starting_model,
        coq9_mode=args.coq9_curation,
    )
    output = args.output_model
    evidence = gene_evidence(model)
    evidence_path = output.with_suffix(".coq9_genes.tsv")
    with evidence_path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(evidence)
    record = {
        "started_at_utc": started, "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "argv": sys.argv, "options": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "head": git("rev-parse", "HEAD"), "dirty_status": git("status", "--porcelain=v1"),
        "input": {"path": str(args.starting_model.resolve()), "sha256": input_sha},
        "input_unchanged": sha(args.starting_model) == input_sha,
        "output": {"path": str(output.resolve()), "sha256": sha(output)},
        "reference_sha256": "bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee",
        "source_sha256": source_sha, "data_sha256": data_sha,
        "source_unchanged_during_build": all(sha(REPO_ROOT / p) == v for p, v in source_sha.items()),
        "model_counts": {"reactions": len(model.reactions), "species": len(model.metabolites), "genes": len(model.genes)},
        "runtime_strain_overlay": None, "medium": dict(model.medium),
        "solver": model.solver.interface.__name__, "python": sys.version,
        "diagnostic_solves": "not run; execution guard enabled" if args.no_solve else "existing diagnostics enabled",
        "network_gene_enrichment": "cache only; network guard enabled" if args.offline else "enabled",
        "coq9": coq9, "retained_reactions": retention, "reaction_selection": selection,
        "requested_build_complete": coq9["requested_mode_complete"] and selection["complete"],
        "gene_evidence": {"path": str(evidence_path.resolve()), "sha256": sha(evidence_path)},
    }
    output.with_suffix(".build.json").write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n")
    return record


def main(argv=None):
    from .cli import main as run_cli
    return run_cli(argv)


if __name__ == "__main__":
    main()
