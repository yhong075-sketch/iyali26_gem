"""
main.py — orchestration entry point for the iYli21 annotation pipeline.
"""

import argparse
import csv
import hashlib
import json
import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from cobra.io import read_sbml_model, write_sbml_model

from .biomass import fix_biomass_reaction
from .config import CACHE_DIR, MNX_DIR, OUTPUT_MODEL_PATH, REPO_ROOT, STARTING_MODEL_PATH
from .coq9 import CURATION_PATH, GENE_EVIDENCE_PATH, apply_coq9_curation, gene_evidence
from .exchange import configure_medium, set_exchange_bounds
from .gaps import DUPLICATE_PAIRS, add_gap_fill_reactions, find_gaps, merge_duplicate_metabolites, report_gaps
from .annotate_reactions_extended import annotate_remaining_reactions
from .ec_annotation import enrich_genes_with_ec
from .genes import annotate_genes
from .idmapping import _enrich_via_idmapping
from .io import load_chem_prop, load_chem_xref, load_mnxm_depr, load_reac_prop, load_reac_xref
from .metabolites import annotate_metabolites, fix_proton_water_balance, normalize_all_annotations
from .patches import add_isozyme_gprs, annotate_isozyme_genes, apply_all_patches, clean_ec_overload, extend_acyl_pool_c161, fill_neutral_formulas, fix_activex_names, fix_ec_code_format, move_tcdb_out_of_ec
from .reactions import annotate_reactions, backfill_reaction_xrefs
from .quinone import (replace_coq6_route_with_coq9, correct_external_ndh2_gpr_and_remove_duplicate,
                      remove_spurious_quinone_branches, apply_reviewed_quinone_step_gprs,
                      run_quinone_step)
from .vlcfa_stereochemistry import correct_er_vlcfa_3r_stereochemistry, verify_er_vlcfa_3r_stereochemistry_target

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starting-model", type=Path, default=STARTING_MODEL_PATH)
    parser.add_argument("--output-model", type=Path, default=OUTPUT_MODEL_PATH)
    parser.add_argument("--coq9-curation", choices=("off", "metadata", "qcycle"), default="metadata",
                        help="metadata (default): guarded annotations/Boolean deduplication; qcycle: explicit R305 proton candidate; off: bypass only these new rules")
    parser.add_argument("--offline", action="store_true", help="Use local data; skip network gene enrichment")
    parser.add_argument("--no-solve", action="store_true", help="Skip FVA and biomass precursor diagnostics; retain construction steps")
    parser.add_argument("--mnx-dir", type=Path, default=MNX_DIR)
    parser.add_argument("--cache-dir", type=Path, default=CACHE_DIR)
    args = parser.parse_args(argv)
    return build_model(args)


def build_model(args):
    starting_model_path, output_model_path = args.starting_model, args.output_model
    mnx_dir, cache_dir = args.mnx_dir, args.cache_dir
    started_at = datetime.now(timezone.utc).isoformat()
    if starting_model_path.resolve() == output_model_path.resolve():
        raise ValueError("Input and output must be separate files")
    if not starting_model_path.exists():
        logger.error(f"Could not find starting model at {starting_model_path}")
        raise FileNotFoundError(starting_model_path)

    logger.info(f"Loading raw model: {starting_model_path.name}")
    input_sha256 = hashlib.sha256(starting_model_path.read_bytes()).hexdigest()
    source_sha256 = {str(path.relative_to(REPO_ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in sorted(Path(__file__).parent.glob("*.py"))}
    model = read_sbml_model(str(starting_model_path))

    mnx_ok = mnx_dir.exists() and (mnx_dir / "chem_xref.tsv").exists()

    if mnx_ok:
        # Load MetaNetX tables once, reuse across functions
        chem_xref = load_chem_xref(mnx_dir / "chem_xref.tsv")
        chem_prop_data = load_chem_prop(mnx_dir / "chem_prop.tsv")
        reac_xref = load_reac_xref(mnx_dir / "reac_xref.tsv")

        reac_prop_path = mnx_dir / "reac_prop.tsv"
        reac_prop = load_reac_prop(reac_prop_path) if reac_prop_path.exists() else None
        if reac_prop is None:
            logger.warning("reac_prop.tsv not found — Strategy C (fingerprint) disabled")

        mnxm_depr = load_mnxm_depr(mnx_dir / "chem_depr.tsv")
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

        # Priority 2b
        logger.info("=== Priority 2b: H+/H2O balance ===")
        fix_proton_water_balance(model)

        # Priority 4a  (Strategy C needs metabolite annotations from 1+2a above)
        logger.info("=== Priority 4a: reaction annotation ===")
        # The supplied R305 name is known to identify the wrong enzyme. Do
        # not let that label create new MNXR/EC assertions before correction.
        name_exclusions = {}
        if args.coq9_curation != "off":
            name_exclusions["R305"] = json.loads(CURATION_PATH.read_text())["rules"]["R305"]["fields"]["name"]["before"]
        annotate_reactions(model, reac_xref, reac_prop, name_exclusions=name_exclusions)
    else:
        logger.warning(
            f"MetaNetX files not found in {mnx_dir}. "
            "Download chem_xref.tsv, chem_prop.tsv, reac_xref.tsv from "
            "https://www.metanetx.org/mnxdoc/mnxref.html"
        )

    # Independent of MetaNetX: the raw source's erroneous (S) identity must
    # not survive an offline rebuild.  Its curated contract fills the exact
    # neutral source tuple before applying the identity correction.
    logger.info("=== Curation: ER VLCFA (3R)-3-hydroxyhexacosanoyl-CoA ===")
    # Resolve case aliases from MetaNetX before authoritative curation removes
    # old identity keys; otherwise final normalization can resurrect them.
    normalize_all_annotations(model)
    n_vlcfa_stereo = correct_er_vlcfa_3r_stereochemistry(model)
    logger.info(
        "  ER VLCFA stereochemistry: metabolites=%d reactions=%d annotations=%d",
        n_vlcfa_stereo["metabolites"],
        n_vlcfa_stereo["reactions"],
        n_vlcfa_stereo["annotations"],
    )

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
    fix_biomass_reaction(model, diagnose=not args.no_solve)

    # Priority 4b — network required, skip if offline
    logger.info("=== Priority 4b: gene annotation via UniProt ===")
    if not args.offline:
        annotate_genes(model)

    # Priority 4c — ncbigene → UniProt ID-mapping for genes still missing uniprot
    logger.info("=== Priority 4c: UniProt ID-mapping (ncbigene → UniProtKB) ===")
    if not args.offline:
        _enrich_via_idmapping(model)

    # Priority 4d — enrich genes with EC numbers via UniProt stream API
    logger.info("=== Priority 4d: gene EC number enrichment via UniProt ===")
    if not args.offline:
        enrich_genes_with_ec(model)

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

    # Existing Q9 chemistry belongs before FVA, independently of the new mode.
    quinone_steps = [run_quinone_step(model, replace_coq6_route_with_coq9)]

    # Priority 5: gap analysis — FVA before gap-fill
    if not args.no_solve:
        logger.info("=== Priority 5: gap analysis (FVA, post-medium) ===")
        gaps = find_gaps(model)
        report_gaps(gaps)
        blocked_before_medium = len(gaps["blocked_reactions"])
        logger.info(f"  Blocked reactions after medium extension: {blocked_before_medium}")

    # Priority 6: gap-fill — insert P0 reactions from gap_fill_prioritized.csv
    gap_fill_csv = REPO_ROOT / "data" / "gap_fill_prioritized.csv"
    if gap_fill_csv.exists():
        logger.info("=== Priority 6: gap-fill reaction insertion (P0) ===")
        add_gap_fill_reactions(
            model,
            csv_path=gap_fill_csv,
            mnx_dir=mnx_dir if mnx_ok else None,
            cache_dir=cache_dir,
        )
        if not args.no_solve:
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

    ec_check = sum(1 for r in model.reactions if isinstance(r.annotation, dict) and 'ec-code' in r.annotation)
    logger.info(f"  DEBUG: reactions with ec-code BEFORE normalize: {ec_check}")

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

    # Isozyme GPR additions: add curated CLIB89 isozyme genes to existing
    # reactions' GPR via 'or' (safe subset only; see patches.add_isozyme_gprs).
    n_gpr_added = add_isozyme_gprs(model)
    logger.info(f"  Isozyme GPR additions: {n_gpr_added} (reaction, gene) pair(s) added")

    # Annotate those newly added genes (they entered after the main gene
    # annotation + SBO steps, so they need sbo / ncbigene / kegg / uniprot here).
    n_gene_annot = annotate_isozyme_genes(model, network=not args.offline)
    logger.info(f"  Isozyme gene annotation: {n_gene_annot} gene(s) annotated")

    # Existing reviewed quinone GPR assembly follows generic annotation so
    # automated enrichment cannot restore the inherited synthome rule.
    for operation in (correct_external_ndh2_gpr_and_remove_duplicate,
                      remove_spurious_quinone_branches, apply_reviewed_quinone_step_gprs):
        quinone_steps.append(run_quinone_step(model, operation))

    # Fill formulas for definite-neutral metabolites (charge=0, unambiguous).
    n_form = fill_neutral_formulas(model)
    logger.info(f"  Neutral formula fill: {n_form} metabolite copy(ies) filled")

    # Lipid chain-menu extension: add C16:1 palmitoleoyl-CoA to the acyl-CoA pools
    # (Y. lipolytica makes ~8% but the pool omitted it). Idempotent.
    n_c161 = extend_acyl_pool_c161(model)
    logger.info(f"  C16:1 acyl-CoA pool extension: {n_c161} pool(s) extended")

    ec_check2 = sum(1 for r in model.reactions if isinstance(r.annotation, dict) and 'ec-code' in r.annotation)
    logger.info(f"  DEBUG: reactions with ec-code AFTER normalize: {ec_check2}")

    logger.info("=== Final gate: ER VLCFA (3R) stereochemistry ===")
    verify_er_vlcfa_3r_stereochemistry_target(model)

    coq9 = apply_coq9_curation(model, args.coq9_curation)
    if not coq9["requested_mode_complete"]:
        logger.warning("CoQ9 mode %s has local conflicts; see build record", args.coq9_curation)
    output_model_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Saving updated model to: {output_model_path.name}")
    write_sbml_model(model, str(output_model_path))
    evidence_path = output_model_path.with_suffix(".coq9_genes.tsv")
    evidence = gene_evidence(model)
    with evidence_path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    def git(*arguments):
        result = subprocess.run(["git", "-C", str(REPO_ROOT), *arguments], capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None
    inputs = [starting_model_path, CURATION_PATH, GENE_EVIDENCE_PATH, REPO_ROOT / "data" / "quinone_pipeline_provenance.json"]
    inputs.extend(REPO_ROOT / "data" / name for name in (
        "er_vlcfa_3r_stereochemistry.json", "ec_overload_audit.csv", "gpr_isozyme_additions.csv",
        "missing_formula_fill.csv", "gap_fill_prioritized.csv", "ncbi/clib89_feature_table.txt", "kegg/yli_genes.tsv"))
    if mnx_ok:
        inputs.extend(mnx_dir / name for name in ("chem_xref.tsv", "chem_prop.tsv", "reac_xref.tsv", "reac_prop.tsv", "chem_depr.tsv"))
    inputs.extend(cache_dir / name for name in ("mnxm_to_bigg_metabolite.json", "gene_locus_tag_map.json"))
    record = {
        "started_at_utc": started_at, "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "argv": sys.argv, "options": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "input": {"path": str(starting_model_path.resolve()), "sha256": input_sha256},
        "input_unchanged": sha(starting_model_path) == input_sha256,
        "output": {"path": str(output_model_path.resolve()), "sha256": sha(output_model_path)},
        "head": git("rev-parse", "HEAD"), "dirty_status": git("status", "--porcelain=v1"),
        "source_sha256": source_sha256,
        "source_unchanged_during_build": all(sha(REPO_ROOT / path) == value for path, value in source_sha256.items()),
        "data_sha256": {str(path.resolve()): sha(path) for path in inputs if path.is_file()},
        "runtime_strain_overlay": None, "medium": "existing set_exchange_bounds/configure_medium",
        "solver": model.solver.interface.__name__, "python": sys.version,
        "diagnostic_solves": "skipped by --no-solve" if args.no_solve else "existing diagnostics enabled",
        "network_gene_enrichment": "skipped by --offline" if args.offline else "enabled",
        "metanetx_available": mnx_ok, "coq9": coq9,
        "existing_quinone_chain": quinone_steps,
        "requested_build_complete": coq9["requested_mode_complete"] and all(s["status"] != "conflict" for s in quinone_steps),
        "gene_evidence": {"path": str(evidence_path.resolve()), "sha256": sha(evidence_path)},
    }
    output_model_path.with_suffix(".build.json").write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n")
    logger.info("Model build complete; CoQ9 requested_mode_complete=%s", coq9["requested_mode_complete"])
    return record


if __name__ == "__main__":
    main()
