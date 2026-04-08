"""
update_model.py — full annotation pipeline for iYli21 GEM.

Priority order:
  1. Metabolite cross-database annotation via MetaNetX chem_xref + chem_prop
  2a. Fill missing chemical formulas from metabolite names / chem_prop
  2b. Balance H+/H2O in reactions
  3. Diagnose & fix biomass reaction R1372
  4a. Reaction annotation via MetaNetX reac_xref
  4b. Gene annotation via UniProt REST API
"""

import logging
import re
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
from cobra.io import read_sbml_model, write_sbml_model

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent
STARTING_MODEL_PATH = REPO_ROOT / "data" / "iyli21.xml"
OUTPUT_MODEL_PATH = REPO_ROOT / "model.xml"
MNX_DIR = REPO_ROOT / "data" / "metanetx"

# ─────────────────────────────────────────────
# MetaNetX loaders  (polars-fast if available, else pandas)
# ─────────────────────────────────────────────

def _read_tsv(path: Path, names: list[str]) -> pd.DataFrame:
    """Read a MetaNetX TSV (skip # comment lines) into a pandas DataFrame."""
    return pd.read_csv(
        path, sep="\t", comment="#", header=None,
        names=names, dtype=str, low_memory=False,
    ).fillna("")


def load_chem_xref(path: Path) -> dict[str, list[tuple[str, str]]]:
    """
    Returns two indexes from chem_xref.tsv:
      by_source[source_key]  → list of mnx_ids
      by_mnxid[mnx_id]      → list of (db_prefix, db_id)
    Packed into a single dict with those two keys.
    """
    logger.info("Loading chem_xref.tsv …")
    df = _read_tsv(path, ["source", "mnx_id", "description"])
    # only keep rows that map to real MNXM IDs (not BIOMASS / EMPTY)
    df = df[df["mnx_id"].str.startswith("MNXM")]

    by_source: dict[str, str] = {}          # source_key → mnx_id (first hit)
    by_mnxid: dict[str, list] = defaultdict(list)  # mnx_id → [(prefix, id), …]

    for source, mnx_id, _ in df.itertuples(index=False):
        if source not in by_source:
            by_source[source] = mnx_id
        if ":" in source:
            prefix, sid = source.split(":", 1)
            by_mnxid[mnx_id].append((prefix, sid))

    logger.info(f"  {len(by_source):,} source entries, {len(by_mnxid):,} MNXM IDs")
    return {"by_source": by_source, "by_mnxid": dict(by_mnxid)}


def load_chem_prop(path: Path) -> dict[str, dict]:
    """
    Returns prop[mnx_id] = {name, formula, charge, inchi, inchikey, smiles}
    Also builds name_index: lower(name) → mnx_id  (first hit per name)
    """
    logger.info("Loading chem_prop.tsv …")
    df = _read_tsv(path, ["mnx_id", "name", "reference", "formula", "charge",
                           "mass", "inchi", "inchikey", "smiles"])
    df = df[df["mnx_id"].str.startswith("MNXM")]

    prop: dict[str, dict] = {}
    name_index: dict[str, str] = {}   # lower(name) → mnx_id

    for row in df.itertuples(index=False):
        prop[row.mnx_id] = {
            "name": row.name, "formula": row.formula, "charge": row.charge,
            "inchi": row.inchi, "inchikey": row.inchikey, "smiles": row.smiles,
        }
        key = row.name.lower().strip()
        if key and key not in name_index:
            name_index[key] = row.mnx_id

    logger.info(f"  {len(prop):,} compounds, {len(name_index):,} unique names")
    return {"prop": prop, "name_index": name_index}


def load_reac_xref(path: Path) -> dict[str, list[tuple[str, str]]]:
    """
    Returns:
      by_mnxr[mnx_id]         → list of (db_prefix, db_id)
      bigg_to_mnxr[bigg_id]   → mnx_id
      desc_index[lower(name)] → mnx_id   (name = part before '||' in description)
    """
    logger.info("Loading reac_xref.tsv …")
    df = _read_tsv(path, ["source", "mnx_id", "description"])
    df = df[df["mnx_id"].str.startswith("MNXR")]

    by_mnxr: dict[str, list] = defaultdict(list)
    bigg_to_mnxr: dict[str, str] = {}
    desc_index: dict[str, str] = {}    # lower(short name) → mnx_id

    for source, mnx_id, desc in df.itertuples(index=False):
        if ":" in source:
            prefix, sid = source.split(":", 1)
            by_mnxr[mnx_id].append((prefix, sid))
            if prefix == "bigg.reaction" and sid not in bigg_to_mnxr:
                bigg_to_mnxr[sid] = mnx_id

        # description format: "Short name||equation" — index the short name
        short = desc.split("||")[0].strip().lower() if "||" in desc else desc.strip().lower()
        if short and short not in desc_index:
            desc_index[short] = mnx_id

    logger.info(f"  {len(by_mnxr):,} MNXR IDs, {len(bigg_to_mnxr):,} BiGG reactions, {len(desc_index):,} names indexed")
    return {"by_mnxr": dict(by_mnxr), "bigg_to_mnxr": bigg_to_mnxr, "desc_index": desc_index}


# ─────────────────────────────────────────────
# Helper: parse name_formula from iYli21 name strings
# ─────────────────────────────────────────────

_FORMULA_RE = re.compile(r"^[A-Z][A-Za-z0-9]*$")   # no *, +, -, spaces

def _parse_name_formula(raw_name: str) -> tuple[str, str]:
    """
    iYli21 name format: 'chemical name_FORMULA'
    Returns (chem_name, formula).  formula may be '' if not parseable.
    """
    if "_" not in raw_name:
        return raw_name.strip(), ""
    chem_name, maybe_formula = raw_name.rsplit("_", 1)
    maybe_formula = maybe_formula.strip()
    # Accept as formula if it looks like a molecular formula
    # or is a charge token like 'p+1', 'H2O', etc.
    if maybe_formula and _FORMULA_RE.match(maybe_formula):
        return chem_name.strip(), maybe_formula
    # charge tokens like 'p+1' — not a formula
    return raw_name.strip(), ""


# ─────────────────────────────────────────────
# Priority 1 + 2a: metabolite annotation + formula fill
# ─────────────────────────────────────────────

def annotate_metabolites(model, chem_xref: dict, chem_prop_data: dict) -> None:
    """
    For every metabolite:
      - Extract chemical formula from name (2a)
      - Lookup MNXM via name → get cross-refs (1) + InChI/InChIKey (2a)
    Modifies model in-place.
    """
    by_source = chem_xref["by_source"]
    by_mnxid = chem_xref["by_mnxid"]
    prop = chem_prop_data["prop"]
    name_index = chem_prop_data["name_index"]

    formula_from_name = 0
    formula_from_prop = 0
    annotated = 0
    no_match = 0

    for met in model.metabolites:
        chem_name, formula_str = _parse_name_formula(met.name)

        # 2a: set formula from name if not already set
        if formula_str and not met.formula:
            met.formula = formula_str
            formula_from_name += 1

        # 1: look up MNXM
        # Strategy A: bigg.metabolite already in annotation (future-proof)
        mnxm_id = None
        bigg_raw = met.annotation.get("bigg.metabolite")
        if bigg_raw:
            bigg_id = bigg_raw[0] if isinstance(bigg_raw, list) else bigg_raw
            mnxm_id = by_source.get(f"bigg.metabolite:{bigg_id}")

        # Strategy B: match by chemical name (case-insensitive)
        if mnxm_id is None:
            mnxm_id = name_index.get(chem_name.lower())

        if mnxm_id is None:
            no_match += 1
            continue

        # Build new annotation dict
        new_ann: dict[str, list] = defaultdict(list)
        for db_prefix, db_id in by_mnxid.get(mnxm_id, []):
            new_ann[db_prefix].append(db_id)
        new_ann["metanetx.chemical"] = [mnxm_id]

        # Enrich with InChI, InChIKey, formula from chem_prop
        if mnxm_id in prop:
            p = prop[mnxm_id]
            if p["inchi"]:
                new_ann["inchi"] = p["inchi"]
            if p["inchikey"]:
                new_ann["inchikey"] = p["inchikey"]
            if p["formula"] and not met.formula and _FORMULA_RE.match(p["formula"]):
                met.formula = p["formula"]
                formula_from_prop += 1

        # Merge (keep existing, add new)
        merged = dict(met.annotation)
        merged.update(new_ann)
        met.annotation = merged
        annotated += 1

    logger.info(
        f"Metabolites: {annotated} annotated, {no_match} unmatched | "
        f"Formulas: {formula_from_name} from name, {formula_from_prop} from prop"
    )


# ─────────────────────────────────────────────
# Priority 2b: balance H+ / H2O
# ─────────────────────────────────────────────

def _get_met_by_name_fragment(model, fragment: str, compartment: str):
    """Find a metabolite matching name fragment in the given compartment."""
    frag_lower = fragment.lower()
    for met in model.metabolites:
        if compartment in met.compartment and frag_lower in met.name.lower():
            return met
    return None


def fix_proton_water_balance(model) -> None:
    """
    For each internal (non-exchange) reaction that is mass-unbalanced due
    to missing H+ or H2O, add the balancing stoichiometry.

    Only touches reactions where the only imbalance is in H and/or O/H atoms
    consistent with H+ (formula H) or H2O (formula H2O).
    """
    fixed = 0
    skipped_no_formula = 0

    for rxn in model.reactions:
        # skip exchange/demand/sink reactions
        if len(rxn.metabolites) <= 1:
            continue
        # skip if any metabolite missing formula
        if any(not met.formula for met in rxn.metabolites):
            skipped_no_formula += 1
            continue

        try:
            balance = rxn.check_mass_balance()
        except Exception:
            continue

        if not balance:
            continue

        elements = set(balance.keys())
        # Only fix if imbalance is purely H and/or O
        if not elements.issubset({"H", "O", "charge"}):
            continue

        # Determine compartment: majority compartment of metabolites
        comps = [met.compartment for met in rxn.metabolites]
        compartment = max(set(comps), key=comps.count)

        h_imb = balance.get("H", 0)
        o_imb = balance.get("O", 0)

        # Case: only H imbalance → add H+
        if h_imb != 0 and o_imb == 0:
            proton = _get_met_by_name_fragment(model, "H+", compartment)
            if proton is None:
                continue
            current = rxn.metabolites.get(proton, 0)
            rxn.add_metabolites({proton: current - h_imb})
            fixed += 1

        # Case: H and O in 2:1 ratio → add H2O
        elif h_imb != 0 and o_imb != 0 and h_imb == 2 * o_imb:
            water = _get_met_by_name_fragment(model, "H2O", compartment)
            if water is None:
                continue
            n_water = o_imb
            current = rxn.metabolites.get(water, 0)
            rxn.add_metabolites({water: current - n_water})
            fixed += 1

    logger.info(
        f"Proton/water balance: {fixed} reactions fixed, "
        f"{skipped_no_formula} skipped (missing formula)"
    )


# ─────────────────────────────────────────────
# Priority 3: diagnose & fix biomass reaction R1372
# ─────────────────────────────────────────────

def fix_biomass_reaction(model) -> None:
    """
    R1372 currently has lb=ub=0 (blocked).
    Steps:
      1. Open bounds so it can carry flux.
      2. Log stoichiometric composition for manual review.
      3. Check if all precursors can be produced (FBA check).
    """
    if "R1372" not in [r.id for r in model.reactions]:
        logger.warning("R1372 not found in model")
        return

    rxn = model.reactions.get_by_id("R1372")
    logger.info(f"R1372 current bounds: [{rxn.lower_bound}, {rxn.upper_bound}]")

    # Open bounds
    if rxn.lower_bound == 0 and rxn.upper_bound == 0:
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
        logger.info("R1372: bounds opened to [0, 1000]")

    # Log composition
    reactants = {m.name: abs(v) for m, v in rxn.metabolites.items() if v < 0}
    products  = {m.name: abs(v) for m, v in rxn.metabolites.items() if v > 0}
    logger.info(f"R1372 reactants ({len(reactants)}): {list(reactants.items())[:8]} …")
    logger.info(f"R1372 products  ({len(products)}):  {list(products.items())[:8]} …")

    # Check mass balance
    try:
        balance = rxn.check_mass_balance()
        if balance:
            logger.warning(f"R1372 mass imbalance: {balance}")
        else:
            logger.info("R1372 mass balance: OK")
    except Exception as e:
        logger.warning(f"R1372 mass balance check failed: {e}")

    # Set as objective and run FBA (use glpk to avoid Gurobi size limits)
    with model:
        model.solver = "glpk"
        model.objective = rxn
        sol = model.optimize()
        obj_val = sol.objective_value if sol.status == "optimal" else float("nan")
        logger.info(f"R1372 FBA (maximize): status={sol.status}, flux={obj_val:.4f}")

    logger.info("R1372: fix complete — review log above for manual curation needs")


# ─────────────────────────────────────────────
# Priority 4a: reaction annotation via reac_xref
# ─────────────────────────────────────────────

def annotate_reactions(model, reac_xref: dict) -> None:
    """
    Match reaction names against MetaNetX reac_xref description column.
    Falls back to bigg.reaction prefix if already annotated.
    """
    by_mnxr = reac_xref["by_mnxr"]
    desc_index = reac_xref["desc_index"]
    bigg_to_mnxr = reac_xref["bigg_to_mnxr"]

    annotated = 0
    no_match = 0

    for rxn in model.reactions:
        mnxr_id = None

        # Strategy A: existing bigg.reaction annotation
        bigg_raw = rxn.annotation.get("bigg.reaction")
        if bigg_raw:
            bid = bigg_raw[0] if isinstance(bigg_raw, list) else bigg_raw
            mnxr_id = bigg_to_mnxr.get(bid)

        # Strategy B: match reaction name to description
        if mnxr_id is None:
            mnxr_id = desc_index.get(rxn.name.lower().strip())

        if mnxr_id is None:
            no_match += 1
            continue

        new_ann: dict[str, list] = defaultdict(list)
        for db_prefix, db_id in by_mnxr.get(mnxr_id, []):
            new_ann[db_prefix].append(db_id)
        new_ann["metanetx.reaction"] = [mnxr_id]

        merged = dict(rxn.annotation)
        merged.update(new_ann)
        rxn.annotation = merged
        annotated += 1

    logger.info(f"Reactions: {annotated} annotated, {no_match} unmatched")


# ─────────────────────────────────────────────
# Priority 4b: gene annotation via UniProt REST
# ─────────────────────────────────────────────

_UNIPROT_BATCH_SIZE = 200
_UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"


def _query_uniprot_batch(gene_ids: list[str]) -> dict[str, dict]:
    """
    Query UniProt for a batch of YALI1 locus IDs.
    Returns {gene_id: {uniprot_id, gene_name, ...}}
    """
    # UniProt stores YALI1 locus IDs as YALI1_Cxxxg (with underscore).
    # The gene: field does a substring/alias match that works for both formats.
    query = " OR ".join(f'(gene:{gid})' for gid in gene_ids)
    params = {
        "query": f"({query}) AND (organism_id:4952)",  # 4952 = Y. lipolytica
        "fields": "accession,gene_names",
        "format": "tsv",
        "size": _UNIPROT_BATCH_SIZE,
    }
    try:
        resp = requests.get(_UNIPROT_URL, params=params, timeout=30)
        resp.raise_for_status()
    except Exception as e:
        logger.warning(f"UniProt request failed: {e}")
        return {}

    result: dict[str, dict] = {}
    lines = resp.text.strip().split("\n")
    if len(lines) < 2:
        return result

    # Build a normalised lookup: strip underscores, lowercase
    def _norm(s: str) -> str:
        return s.replace("_", "").lower()

    norm_to_orig = {_norm(gid): gid for gid in gene_ids}

    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        acc = parts[0].strip()
        gene_names_raw = parts[1].strip()
        for gname in gene_names_raw.split():
            orig = norm_to_orig.get(_norm(gname))
            if orig and orig not in result:
                result[orig] = {"uniprot": acc}
    return result


def annotate_genes(model) -> None:
    """
    Map YALI1* locus IDs to UniProt accessions via the UniProt REST API.
    Adds annotation["uniprot"] = [accession].
    """
    gene_ids = [g.id for g in model.genes]
    logger.info(f"Querying UniProt for {len(gene_ids)} genes in batches …")

    mapping: dict[str, dict] = {}
    for i in range(0, len(gene_ids), _UNIPROT_BATCH_SIZE):
        batch = gene_ids[i: i + _UNIPROT_BATCH_SIZE]
        result = _query_uniprot_batch(batch)
        mapping.update(result)
        time.sleep(0.5)   # be polite to UniProt

    annotated = 0
    for gene in model.genes:
        if gene.id in mapping:
            merged = dict(gene.annotation)
            merged["uniprot"] = [mapping[gene.id]["uniprot"]]
            gene.annotation = merged
            annotated += 1

    logger.info(f"Genes: {annotated}/{len(gene_ids)} mapped to UniProt")


# ─────────────────────────────────────────────
# main
# ─────────────────────────────────────────────

def main():
    if not STARTING_MODEL_PATH.exists():
        logger.error(f"Could not find starting model at {STARTING_MODEL_PATH}")
        return

    logger.info(f"Loading raw model: {STARTING_MODEL_PATH.name}")
    model = read_sbml_model(str(STARTING_MODEL_PATH))

    mnx_ok = MNX_DIR.exists() and (MNX_DIR / "chem_xref.tsv").exists()

    if mnx_ok:
        # Load MetaNetX tables once, reuse across functions
        chem_xref = load_chem_xref(MNX_DIR / "chem_xref.tsv")
        chem_prop_data = load_chem_prop(MNX_DIR / "chem_prop.tsv")
        reac_xref = load_reac_xref(MNX_DIR / "reac_xref.tsv")

        # Priority 1 + 2a
        logger.info("=== Priority 1+2a: metabolite annotation + formulas ===")
        annotate_metabolites(model, chem_xref, chem_prop_data)

        # Priority 2b
        logger.info("=== Priority 2b: H+/H2O balance ===")
        fix_proton_water_balance(model)

        # Priority 4a
        logger.info("=== Priority 4a: reaction annotation ===")
        annotate_reactions(model, reac_xref)
    else:
        logger.warning(
            f"MetaNetX files not found in {MNX_DIR}. "
            "Download chem_xref.tsv, chem_prop.tsv, reac_xref.tsv from "
            "https://www.metanetx.org/mnxdoc/mnxref.html"
        )

    # Priority 3 (always run — independent of MetaNetX)
    logger.info("=== Priority 3: biomass reaction R1372 ===")
    fix_biomass_reaction(model)

    # Priority 4b — network required, skip if offline
    logger.info("=== Priority 4b: gene annotation via UniProt ===")
    annotate_genes(model)

    logger.info(f"Saving updated model to: {OUTPUT_MODEL_PATH.name}")
    write_sbml_model(model, str(OUTPUT_MODEL_PATH))
    logger.info("Model build complete.")


if __name__ == "__main__":
    main()
