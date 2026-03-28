import logging
import re
from pathlib import Path

import pandas as pd
from cobra.io import read_sbml_model, write_sbml_model

# Setup clean logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Dynamically resolve paths so this script works anywhere
REPO_ROOT = Path(__file__).resolve().parent.parent
STARTING_MODEL_PATH = REPO_ROOT / "data" / "iyli21.xml"
OUTPUT_MODEL_PATH = REPO_ROOT / "model.xml"
DATA_DIR = REPO_ROOT / "data"


# ──────────────────────────────────────────
# Step 0: Normalize compartment IDs
# ──────────────────────────────────────────

COMPARTMENT_MAP = {
    "C_cy": "c",   # cytoplasm
    "C_mi": "m",   # mitochondria
    "C_ex": "e",   # extracellular
    "C_er": "r",   # endoplasmic reticulum
    "C_nu": "n",   # nucleus
    "C_go": "g",   # Golgi apparatus
    "C_pe": "x",   # peroxisome
    "C_va": "v",   # vacuole
    "C_lp": "l",   # lipid droplet
    "C_en": "ce",  # cell envelope / cell wall
    "C_em": "em",  # ER membrane
    "C_gm": "gm",  # Golgi membrane
    "C_mm": "mm",  # mitochondrial membrane
    "C_vm": "vm",  # vacuolar membrane
    "c_va": "v",   # typo for C_va (vacuole)
}

COMPARTMENT_NAMES = {
    "c":  "cytoplasm",
    "m":  "mitochondria",
    "e":  "extracellular",
    "r":  "endoplasmic reticulum",
    "n":  "nucleus",
    "g":  "Golgi apparatus",
    "x":  "peroxisome",
    "v":  "vacuole",
    "l":  "lipid droplet",
    "ce": "cell envelope",
    "em": "ER membrane",
    "gm": "Golgi membrane",
    "mm": "mitochondrial membrane",
    "vm": "vacuolar membrane",
}


def normalize_compartments(model):
    """Rename non-standard compartment IDs to BiGG-style short IDs using COBRApy."""
    logger.info("Step 0: Normalizing compartment IDs")

    # Rename metabolite IDs: COBRApy loads m1__91__C_cy__93__ as m1[C_cy],
    # so we strip the bracketed suffix and append the new BiGG compartment ID.
    met_count = 0
    for met in model.metabolites:
        old_comp = met.compartment
        new_comp = COMPARTMENT_MAP.get(old_comp, old_comp)
        met.id = re.sub(r'\[.*\]$', '', met.id) + '_' + new_comp
        met.compartment = new_comp
        met_count += 1

    # Update compartment definitions
    new_compartments = {}
    for cid, name in model.compartments.items():
        new_cid = COMPARTMENT_MAP.get(cid, cid)
        new_compartments[new_cid] = COMPARTMENT_NAMES.get(new_cid, name or "")
    model.compartments = new_compartments

    logger.info(f"  Metabolites renamed: {met_count}, compartments: {list(model.compartments.keys())}")
    return model


# ──────────────────────────────────────────
# Step 1: Write metabolite annotations back
# ──────────────────────────────────────────

def apply_metabolite_annotations(model):
    """Apply BiGG/KEGG/LLM annotations from final_annotation_table.csv."""
    logger.info("Step 1: Applying metabolite annotations")

    annot_path = DATA_DIR / "final_annotation_table.csv"
    raw_path   = DATA_DIR / "metabolites_raw.csv"

    if not annot_path.exists() or not raw_path.exists():
        logger.warning("  Annotation CSV files not found, skipping step 1")
        return model

    annot_df = pd.read_csv(annot_path).set_index("base_id")
    raw_df   = pd.read_csv(raw_path)
    met_to_base = dict(zip(raw_df["met_id"], raw_df["base_id"]))

    updated = 0
    for met in model.metabolites:
        base_id = met_to_base.get(met.id)
        if base_id is None or base_id not in annot_df.index:
            continue
        row = annot_df.loc[base_id]
        if pd.notna(row.get("name")):
            met.name = str(row["name"])
        if pd.notna(row.get("formula")):
            met.formula = str(row["formula"])
        if pd.notna(row.get("charge")):
            met.charge = int(row["charge"])
        if pd.notna(row.get("bigg_id")):
            met.annotation["bigg.metabolite"] = str(row["bigg_id"])
        if pd.notna(row.get("kegg_id")):
            met.annotation["kegg.compound"] = str(row["kegg_id"])
        updated += 1

    logger.info(f"  Updated {updated} / {len(model.metabolites)} metabolites")
    return model


# ──────────────────────────────────────────
# Step 2: Quick fixes for memote score gains
# ──────────────────────────────────────────

ID_PATTERNS = {
    "kegg.compound":      re.compile(r'^C\d{5}$'),
    "kegg.reaction":      re.compile(r'^R\d{5}$'),
    "kegg.genes":         re.compile(r'^[a-z]{3,4}:\S+$'),
    "chebi":              re.compile(r'^CHEBI:\d+$'),
    "hmdb":               re.compile(r'^HMDB\d+$'),
    "metanetx.chemical":  re.compile(r'^MNXM\d+$'),
    "metanetx.reaction":  re.compile(r'^MNXR\d+$'),
    "seed.compound":      re.compile(r'^cpd\d+$'),
    "seed.reaction":      re.compile(r'^rxn\d+$'),
    "bigg.metabolite":    re.compile(r'^[a-zA-Z0-9_]+$'),
    "bigg.reaction":      re.compile(r'^[a-zA-Z0-9_]+$'),
    "uniprot":            re.compile(r'^[A-Z0-9]{6,10}$'),
    "ncbigene":           re.compile(r'^\d+$'),
    "ec-code":            re.compile(r'^\d+\.\d+\.\d+\.\d+$'),
    "inchikey":           re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'),
    "reactome":           re.compile(r'^\d+$|^R-[A-Z]{3}-\d+$'),
    "biocyc":             re.compile(r'^.+:.+$'),
    "rhea":               re.compile(r'^\d+$'),
}


def apply_quick_fixes(model):
    """Add SBO terms to biomass/exchange reactions and clean invalid annotation IDs."""
    logger.info("Step 2: Applying quick fixes")

    # Add SBO to biomass reactions
    for rxn in model.reactions:
        if "biomass" in rxn.id.lower() or "biomass" in (rxn.name or "").lower():
            rxn.annotation["sbo"] = "SBO:0000629"

    # Add SBO to exchange/demand/sink reactions
    sbo_added = 0
    for rxn in model.reactions:
        if "sbo" in rxn.annotation:
            continue
        rid = rxn.id.lower()
        if rid.startswith("ex_") or rid.startswith("ex "):
            rxn.annotation["sbo"] = "SBO:0000627"
            sbo_added += 1
        elif rid.startswith("dm_"):
            rxn.annotation["sbo"] = "SBO:0000628"
            sbo_added += 1
        elif rid.startswith("sk_") or rid.startswith("sink_"):
            rxn.annotation["sbo"] = "SBO:0000632"
            sbo_added += 1

    # Clean invalid annotation IDs
    removed = 0
    for obj in list(model.metabolites) + list(model.reactions) + list(model.genes):
        keys_to_remove = [
            key for key, val in obj.annotation.items()
            if key in ID_PATTERNS and not ID_PATTERNS[key].match(str(val))
        ]
        for key in keys_to_remove:
            del obj.annotation[key]
            removed += 1

    logger.info(f"  SBO added to {sbo_added} exchange/demand/sink reactions")
    logger.info(f"  Removed {removed} invalid annotation IDs")
    return model


# ──────────────────────────────────────────
# Main
# ──────────────────────────────────────────

def main():
    if not STARTING_MODEL_PATH.exists():
        logger.error(f"Could not find starting model at {STARTING_MODEL_PATH}")
        return

    # Load model with COBRApy
    logger.info(f"Loading model: {STARTING_MODEL_PATH.name}")
    model = read_sbml_model(str(STARTING_MODEL_PATH))

    # Step 0: normalize compartment IDs
    model = normalize_compartments(model)

    # Step 1: write metabolite annotations
    model = apply_metabolite_annotations(model)

    # Step 2: quick fixes
    model = apply_quick_fixes(model)

    logger.info(f"Saving updated model to: {OUTPUT_MODEL_PATH.name}")
    write_sbml_model(model, str(OUTPUT_MODEL_PATH))
    logger.info("Model build complete.")


if __name__ == "__main__":
    main()
