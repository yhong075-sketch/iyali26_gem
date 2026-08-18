"""
patches.py — known data-bug patches for the iYli21 model.

Each patch fixes a discrete, documented error in the source SBML.  Patches
are applied early in the pipeline, immediately after metabolite annotation,
so that downstream mass-balance, FBA, and gap analyses operate on a
chemically correct model.

Currently applied patches:

  1. NADP+ formula fix
     iYli21 stores NADP+ as C21H28N7O17P3 (one H short of the KEGG C00006
     reference C21H29N7O17P3).  This affects 6 compartmental copies and
     ~126 reactions that use NADP+/NADPH, which had cascading effects on
     ceramide synthesis (see patch 2).

  2. Ceramide formula corrections
     iYli21 stores ceramide-1-(C24/C26) with formula C63H125NO6, which is
     inconsistent with the LIPID MAPS / Yeast-GEM reference value
     C42H85NO3 (C24) / C44H89NO3 (C26).  ceramide-2/2'/3/4 are simply
     missing formulas.  These were left wrong/missing because the upstream
     NADP+ error made the relevant reactions unbalanceable — fixing patch 1
     now makes patch 2 work.

     Naming map (verified against LIPID MAPS):
       ceramide-1  = Cer(d18:0/24:0)         = C42H85NO3  / C44H89NO3
       ceramide-2  = Cer(t18:0/24:0)         = C42H85NO4  / C44H89NO4
       ceramide-2' = Cer(d18:0/24:0(2OH))    = C42H85NO4  / C44H89NO4
       ceramide-3  = Cer(t18:0/24:0(2OH))    = C42H85NO5  / C44H89NO5
       ceramide-4  = (extrapolated +1 O from ceramide-3)
                                              C42H85NO6  / C44H89NO6
"""

from __future__ import annotations

import hashlib
import json
import logging
import re
from collections import Counter
from copy import deepcopy
from math import isclose, isfinite
from numbers import Real
from pathlib import Path

logger = logging.getLogger(__name__)


# ── Patch 1: NADP+ formula fix ────────────────────────────────────────────

# iYli21 NADP+ copies (all compartments). Identified by metabolite name
# prefix "NADP(+)_" which is the iYli21 naming convention.
# Formula correction: C21H28N7O17P3 → C21H29N7O17P3 (KEGG C00006).
_NADP_PLUS_OLD_FORMULA = "C21H28N7O17P3"
_NADP_PLUS_NEW_FORMULA = "C21H29N7O17P3"


def _is_nadp_plus(met) -> bool:
    """Identify NADP+ by name (handles iYli21 'NADP(+)_C21H28N7O17P3' naming)."""
    name = (met.name or "")
    return name.startswith("NADP(+)") or name.lower().startswith("nadp(+)")


def fix_nadp_plus_formula(model) -> int:
    """
    Fix the NADP+ formula bug across all compartments.

    Returns the number of metabolites patched.
    """
    fixed = 0
    for met in model.metabolites:
        if not _is_nadp_plus(met):
            continue
        if met.formula == _NADP_PLUS_OLD_FORMULA:
            met.formula = _NADP_PLUS_NEW_FORMULA
            fixed += 1
            logger.debug(f"  NADP+ patch: {met.id} formula → {_NADP_PLUS_NEW_FORMULA}")
        elif met.formula == _NADP_PLUS_NEW_FORMULA:
            pass   # already correct
        else:
            logger.warning(
                f"  NADP+ patch: {met.id} has unexpected formula {met.formula!r} "
                f"— skipped (expected {_NADP_PLUS_OLD_FORMULA})"
            )
    return fixed


# ── Patch 2: Ceramide formula corrections ─────────────────────────────────

# Map: lowercased base name (compartment-independent) → (formula_C24, formula_C26)
# Base name is what comes BEFORE the trailing "_" or "_FORMULA" in iYli21 names.
# E.g. "ceramide-1-(C24)_" matches base "ceramide-1-(c24)".
_CERAMIDE_FORMULAS_C24 = {
    "ceramide-1-(c24)":  "C42H85NO3",   # Cer(d18:0/24:0)        LMSP02020012
    "ceramide-2-(c24)":  "C42H85NO4",   # Cer(t18:0/24:0)        LMSP02030004
    "ceramide-2'-(c24)": "C42H85NO4",   # Cer(d18:0/24:0(2OH))   (α-OH dihydroceramide)
    "ceramide-3-(c24)":  "C42H85NO5",   # Cer(t18:0/24:0(2OH))   LMSP02030002
    "ceramide-4-(c24)":  "C42H85NO6",   # (extrapolated, R208 product)
}
_CERAMIDE_FORMULAS_C26 = {
    "ceramide-1-(c26)":  "C44H89NO3",
    "ceramide-2-(c26)":  "C44H89NO4",
    "ceramide-2'-(c26)": "C44H89NO4",
    "ceramide-3-(c26)":  "C44H89NO5",
    "ceramide-4-(c26)":  "C44H89NO6",
}
_CERAMIDE_FORMULAS = {**_CERAMIDE_FORMULAS_C24, **_CERAMIDE_FORMULAS_C26}

# Old (wrong) formula on ceramide-1 — explicitly remove this value when we
# overwrite, so the log warns if it's anything else (defensive).
_CERAMIDE_1_KNOWN_BAD_FORMULA = "C63H125NO6"


def _base_name(met) -> str:
    """Return the iYli21 base name (lowercased, strip trailing _ and _FORMULA)."""
    n = (met.name or "").lower().rstrip("_").strip()
    # Strip trailing "_<FORMULA>" if present (e.g. "ceramide-1-(c24)_c63h125no6")
    if "_" in n:
        parts = n.rsplit("_", 1)
        last = parts[-1]
        # heuristic: looks like a chemical formula → strip
        if last and last[0] in "cChH" and any(c.isdigit() for c in last):
            n = parts[0].rstrip("_")
    return n


def fix_ceramide_formulas(model) -> int:
    """
    Fix ceramide-1 wrong formula + fill ceramide-2/2'/3/4 missing formulas.

    Authoritative source: LIPID MAPS LMSD (verified manually) cross-checked
    against Yeast-GEM (https://github.com/SysBioChalmers/yeast-GEM).

    Applies to all compartmental copies (C_er, C_mi, C_go, etc).
    Returns number of metabolites patched.
    """
    fixed = 0
    seen_bases = set()
    for met in model.metabolites:
        base = _base_name(met)
        if base not in _CERAMIDE_FORMULAS:
            continue
        new_formula = _CERAMIDE_FORMULAS[base]
        old_formula = met.formula or ""
        if old_formula == new_formula:
            continue   # already correct
        if old_formula and old_formula != _CERAMIDE_1_KNOWN_BAD_FORMULA:
            logger.warning(
                f"  Ceramide patch: {met.id} ({base!r}) has unexpected existing "
                f"formula {old_formula!r}, overwriting with {new_formula!r}"
            )
        met.formula = new_formula
        fixed += 1
        seen_bases.add(base)
        logger.debug(
            f"  Ceramide patch: {met.id} [{met.compartment}] {old_formula or '(none)'} "
            f"→ {new_formula}  ({base})"
        )
    # Report which base names did not appear in the model
    missing = set(_CERAMIDE_FORMULAS) - seen_bases
    if missing:
        logger.info(
            f"  Ceramide patch: {len(missing)} mapped base name(s) not found in model "
            f"(possibly C24 or C26 variants absent): {sorted(missing)[:5]}…"
        )
    return fixed


# ── Patch 3: free CoA charge fix ──────────────────────────────────────────

# iYli21 stores free coenzyme A (KEGG C00010) with charge 0, but the correct
# physiological charge is -4.  Only the 9 free-CoA copies are affected; the
# acyl-CoA species (tetracosanoyl-CoA, acetyl-CoA, …) already carry charge -4.
# Identified by: name starts with "coenzyme A", formula == C21H36N7O16P3S, charge == 0.
_COA_FREE_FORMULA   = "C21H36N7O16P3S"
_COA_CORRECT_CHARGE = -4


def fix_coa_charge(model) -> int:
    """
    Fix the free-CoA charge bug (KEGG C00010: charge should be -4, not 0).

    Patches only free coenzyme A; acyl-CoA species are left untouched because
    they already carry the correct charge.  Returns the number patched.
    """
    fixed = 0
    for met in model.metabolites:
        name = met.name or ""
        if not name.startswith("coenzyme A"):
            continue
        if met.formula == _COA_FREE_FORMULA and met.charge == 0:
            met.charge = _COA_CORRECT_CHARGE
            fixed += 1
            logger.debug(f"  CoA patch: {met.id} charge → {_COA_CORRECT_CHARGE}")
    return fixed


# ── Patch 6: cation formula/charge self-consistency ───────────────────────

# A set of protonated amines/cations carry a NEUTRAL formula but a positive
# charge — the two are mutually inconsistent (a +1 ammonium must have one more
# H than its neutral form).  This unbalances every reaction they participate
# in (e.g. sphinganine's 8 sphingolipid reactions).  Fix the formula to the
# protonated (ionic) form, matching the MetaNetX chem_prop value for the
# metabolite's own MNXM.
#
# Whitelist keyed by (current_formula, charge) → target_formula, so a same-
# formula species with a different charge is never touched.  Each target is
# the MetaNetX formula for that metabolite's MNXM (verified 2026-06).
_CATION_FORMULA_FIX: dict[tuple[str, int], str] = {
    ("C11H18N4O4", 1): "C11H19N4O4",   # MNXM735139  2-[3-carboxy-3-(methylammonio)propyl]-L-his
    ("C15H22N6O5S", 1): "C15H23N6O5S", # MNXM1363767 S-adenosyl-L-methionine
    ("C18H39NO2", 1): "C18H40NO2",     # MNXM733692  sphinganine
    ("C2H7NO", 1): "C2H8NO",           # MNXM218     ethanolamine
    ("C3H7NO", 1): "C3H8NO",           # MNXM736082  3-aminopropanal
    ("C4H12N2", 2): "C4H14N2",         # MNXM118     putrescine (+2)
    ("C4H9NO", 1): "C4H10NO",          # MNXM422     4-aminobutanal
    ("C5H12N4O", 1): "C5H13N4O",       # MNXM2617    4-guanidinobutanamide
    ("C6H11N3O", 1): "C6H12N3O",       # MNXM1281    L-histidinol
    ("C6H13NO2S", 1): "C6H14NO2S",     # MNXM681265  S-methyl-L-methionine
    ("C6H14N2O2", 1): "C6H15N2O2",     # MNXM1364268 L-lysine
    ("C6H14N4O2", 1): "C6H15N4O2",     # MNXM739527  L-arginine
    ("C7H19N3", 3): "C7H22N3",         # MNXM124     spermidine (+3)
    ("C8H12N2O2", 1): "C8H13N2O2",     # MNXM548     pyridoxamine
}


def fix_cation_formula_consistency(model) -> int:
    """
    Fix protonated cations whose formula is the neutral form but charge is
    positive.  Only metabolites matching both (formula, charge) in the
    whitelist are changed.  Returns the number of metabolites patched.
    """
    fixed = 0
    for met in model.metabolites:
        if met.formula is None or met.charge is None:
            continue
        key = (met.formula.strip(), int(met.charge))
        target = _CATION_FORMULA_FIX.get(key)
        if target:
            met.formula = target
            fixed += 1
            logger.debug(f"  cation formula patch: {met.id} {key[0]} → {target}")
    return fixed


# ── Patch 0: clean Excel-corrupted "ActiveX VT_ERROR" names ───────────────

# iYli21 was exported from Excel; some metabolite names had their trailing
# chemical-formula token replaced by the literal string "ActiveX VT_ERROR:",
# e.g. "butyrate_ActiveX VT_ERROR:".  The corrupted name prevents
# annotate_metabolites from matching the metabolite to MetaNetX (no MNXM, no
# formula).  Stripping the garbage suffix restores a clean name that the
# annotation step can match.
#
# MUST run BEFORE annotate_metabolites so the clean name is used for matching
# — so it is invoked from main() before annotation, not from apply_all_patches.
_ACTIVEX_RE = re.compile(r"_?ActiveX VT_ERROR:?\s*$")


def fix_activex_names(model) -> int:
    """
    Strip the Excel-corruption suffix '_ActiveX VT_ERROR:' from metabolite
    names.  Returns the number of names cleaned.
    """
    fixed = 0
    for met in model.metabolites:
        name = met.name or ""
        if "ActiveX" in name or "VT_ERROR" in name:
            cleaned = _ACTIVEX_RE.sub("", name).rstrip("_ ").strip()
            if cleaned and cleaned != name:
                met.name = cleaned
                fixed += 1
                logger.debug(f"  ActiveX name patch: {met.id} → {cleaned!r}")
    return fixed


# ── Patch 5: move misfiled TCDB numbers out of ec-code ────────────────────

# Some transport reactions carry a TCDB (Transporter Classification Database)
# number in their 'ec-code' annotation, e.g. "2.A.1.44.1" or "9.A.6.1.1".
# These are not EC numbers (second segment is a letter) and fail Memote's
# ec-code conformity check.  They are valid identifiers, just in the wrong
# field — move them to the identifiers.org-registered 'tcdb' annotation.
#
# NOTE: this must run after all EC codes are populated, so it is invoked from
# main() (before fix_ec_code_format), not from apply_all_patches.
_TCDB_RE = re.compile(r"^\d+\.[A-Z]\.")   # TCDB class id, e.g. 2.A.1.44.1


def move_tcdb_out_of_ec(model) -> int:
    """
    Move TCDB transporter numbers misfiled in 'ec-code' to the 'tcdb'
    annotation field.  Returns the number of TCDB ids moved.
    """
    moved = 0
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        ecs = ann.get("ec-code")
        if not ecs:
            continue
        if isinstance(ecs, str):
            ecs = [ecs]
        tcdb = [e for e in ecs if _TCDB_RE.match(e)]
        if not tcdb:
            continue
        kept = [e for e in ecs if not _TCDB_RE.match(e)]
        existing = ann.get("tcdb", [])
        if isinstance(existing, str):
            existing = [existing]
        ann["tcdb"] = sorted(set(existing) | set(tcdb))
        if kept:
            ann["ec-code"] = kept
        else:
            ann.pop("ec-code", None)
        rxn.annotation = ann
        moved += len(tcdb)
        logger.debug(f"  TCDB cleanup: {rxn.id} moved {tcdb} ec-code → tcdb")
    return moved


# ── Patch 4: EC-code format compliance ───────────────────────────────────

# Some reactions carry EC codes missing the final level (exactly three
# numeric segments, e.g. "3.1.3" or "1.2.1").  These violate the
# identifiers.org EC regex and are flagged by Memote's
# test_reaction_annotation_wrong_ids (ec-code conformity).  Padding them with
# a trailing ".-" makes them valid partial EC codes (e.g. "3.1.3.-") without
# inventing a more specific level we cannot verify.
#
# Only "exactly three pure-numeric segments" are touched.  EC strings with
# letters (preliminary IDs like "2.7.1.M29", malformed "7.4.2.i") or
# misfiled non-EC identifiers (e.g. TCDB number "2.A.29.8.3") are NOT padded
# — those need case-by-case handling and are left untouched (logged).
_EC_THREE_SEGMENT_RE = re.compile(r"^\d+\.\d+\.\d+$")


def fix_ec_code_format(model) -> int:
    """
    Pad three-segment numeric EC codes with a trailing '.-' for
    identifiers.org compliance (e.g. '3.1.3' -> '3.1.3.-').

    Only exact three-segment pure-numeric EC codes are modified.  EC strings
    containing letters or extra/fewer segments are left untouched.  Returns
    the number of individual EC strings padded.
    """
    fixed = 0
    for rxn in model.reactions:
        ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
        ecs = ann.get("ec-code")
        if not ecs:
            continue
        if isinstance(ecs, str):
            ecs = [ecs]
        new = []
        changed = False
        for ec in ecs:
            if _EC_THREE_SEGMENT_RE.match(ec):
                new.append(ec + ".-")
                fixed += 1
                changed = True
                logger.debug(f"  EC-format patch: {rxn.id} {ec} → {ec}.-")
            else:
                new.append(ec)
        if changed:
            rxn.annotation["ec-code"] = new
    return fixed


# ── Patch 7: EC-overload cleanup ──────────────────────────────────────────

# Some reactions accumulated 5+ EC numbers from MetaNetX MNXR `classifs`
# back-fill — a single, well-defined reaction tagged with EC numbers spanning
# multiple enzyme classes (e.g. R1893 mannitol dehydrogenase carried 11 EC
# including glutathione transferase 2.5.1.18).  This "EC soup" pollutes
# EC-based analysis and inflates Memote's EC inconsistency.
#
# The authoritative fix uses KEGG reaction-level ENZYME data: keep only the
# intersection of the reaction's current EC set with the EC numbers KEGG
# assigns to its kegg.reaction id.  That curation is done offline by
# scripts/audit_ec_overload.py, which writes data/ec_overload_audit.csv with
# a per-reaction keep/drop decision.  This patch *applies* the rows marked
# action=clean — it never invents EC numbers (keep_ec is always a subset of
# the reaction's current ec-code, double-checked below).
#
# Runs LAST in the pipeline (after all EC back-fill and format steps) so the
# back-fill cannot re-introduce the dropped EC numbers.
_EC_OVERLOAD_AUDIT_CSV = "data/ec_overload_audit.csv"


def clean_ec_overload(model, audit_csv: str | None = None) -> int:
    """
    Apply the KEGG-curated EC-overload cleanup from the audit CSV.

    For each audit row with action=='clean', replace the reaction's 'ec-code'
    with the curated keep_ec set — but ONLY if keep_ec is a subset of the
    reaction's current ec-code (guards against the CSV drifting out of sync
    with the model).  Returns the number of reactions cleaned.
    """
    import csv
    import os

    if audit_csv is None:
        # patches.py is at scripts/gem_annotate/patches.py -> repo root is ../../
        root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        audit_csv = os.path.join(root, _EC_OVERLOAD_AUDIT_CSV)

    if not os.path.exists(audit_csv):
        logger.warning(f"  EC-overload audit CSV not found: {audit_csv} — skipping")
        return 0

    cleaned = 0
    with open(audit_csv) as f:
        for row in csv.DictReader(f):
            if row.get("action") != "clean":
                continue
            rid = row["reaction"]
            keep = {e for e in row["keep_ec"].split(";") if e}
            if not keep:
                continue
            try:
                rxn = model.reactions.get_by_id(rid)
            except KeyError:
                logger.warning(f"  EC-overload: reaction {rid} not in model — skipping")
                continue
            ann = rxn.annotation if isinstance(rxn.annotation, dict) else {}
            cur = ann.get("ec-code")
            if not cur:
                continue
            if isinstance(cur, str):
                cur = [cur]
            cur_set = set(cur)
            # subset guard: only proceed if keep_ec really is a subset of current
            if not keep <= cur_set:
                logger.warning(
                    f"  EC-overload: {rid} keep_ec {sorted(keep)} not subset of "
                    f"current {sorted(cur_set)} — skipping (CSV out of sync)")
                continue
            if cur_set == keep:
                continue  # already clean, nothing to drop
            rxn.annotation["ec-code"] = sorted(keep)
            cleaned += 1
            logger.debug(
                f"  EC-overload: {rid} {sorted(cur_set)} → {sorted(keep)}")
    return cleaned


# ── Patch 8: isozyme GPR additions ────────────────────────────────────────

# CLIB89 expansion funnel (S2 Table -> KEGG metabolic relevance -> MetaNetX
# reaction mapping -> de-dup against model EC) identified genes whose EC is
# already carried by an existing model reaction: these are isozymes that
# should be added to that reaction's GPR (NOT new reactions).
#
# This patch applies only the SAFE subset curated in
# data/gpr_isozyme_additions.csv: each gene maps to <=3 reactions and none of
# those reactions has an 'and' (multi-subunit complex) GPR — so adding the
# gene with 'or' is unambiguous.  Genes hitting broad ECs (many reactions) or
# complex GPRs are excluded and left for manual review.
#
# Only 'or' additions are made; existing genes and reaction stoichiometry are
# never touched.  Idempotent: a gene already in the rule is skipped.
_GPR_ADDITIONS_CSV = "data/gpr_isozyme_additions.csv"


def add_isozyme_gprs(model, additions_csv: str | None = None) -> int:
    """
    Add curated isozyme genes to existing reactions' GPR via 'or'.

    Reads data/gpr_isozyme_additions.csv (columns: reaction, add_gene, ...).
    For each (reaction, gene): if the gene is not already in the reaction's
    gene_reaction_rule, append it with 'or' (or set it as the sole rule if the
    reaction had no GPR).  Returns the number of (reaction, gene) additions made.
    """
    import csv
    import os

    if additions_csv is None:
        root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        additions_csv = os.path.join(root, _GPR_ADDITIONS_CSV)

    if not os.path.exists(additions_csv):
        logger.warning(f"  GPR additions CSV not found: {additions_csv} — skipping")
        return 0

    added = 0
    with open(additions_csv) as f:
        for row in csv.DictReader(f):
            rid = row["reaction"]
            gene = row["add_gene"].strip()
            if not gene:
                continue
            try:
                rxn = model.reactions.get_by_id(rid)
            except KeyError:
                logger.warning(f"  GPR add: reaction {rid} not in model — skipping")
                continue
            rule = rxn.gene_reaction_rule.strip()
            # idempotent: skip if gene already present as a token
            existing = set(re.findall(r"[A-Za-z0-9_]+", rule))
            if gene in existing:
                continue
            new_rule = gene if not rule else f"{rule} or {gene}"
            rxn.gene_reaction_rule = new_rule
            added += 1
            logger.debug(f"  GPR add: {rid} += {gene}")
    return added


# ── Patch 9: annotate the isozyme genes added by patch 8 ───────────────────

# The genes added by add_isozyme_gprs enter the model with no annotation.
# The pipeline's main gene-annotation step (genes.annotate_genes) and SBO step
# both run BEFORE the GPR additions, so they never reach these new genes —
# leaving them with empty annotation (regresses memote gene-SBO /
# gene-product-annotation, and they re-appear empty on every full rebuild).
#
# This patch annotates exactly those genes: sbo (SBO:0000243), ncbigene +
# kegg.genes (from NCBI feature table / KEGG yli, local), and uniprot
# (best-effort network lookup via xref:geneid).  Runs right after
# add_isozyme_gprs.  Idempotent: genes already carrying an sbo are skipped.
_FEATURE_TABLE = "data/ncbi/clib89_feature_table.txt"
_KEGG_GENES = "data/kegg/yli_genes.tsv"


def _fetch_uniprot_for_geneid(geneid: str):
    """Best-effort UniProt accession from a GeneID xref. None on failure."""
    import urllib.parse
    import urllib.request
    q = urllib.parse.quote(f"xref:geneid-{geneid}")
    url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
           f"&fields=accession&format=tsv&size=1")
    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            lines = resp.read().decode().splitlines()
        if len(lines) >= 2 and lines[1].strip():
            return lines[1].strip()
    except Exception:
        pass
    return None


def annotate_isozyme_genes(model, additions_csv: str | None = None,
                           network: bool = True) -> int:
    """
    Annotate the isozyme genes added by add_isozyme_gprs.

    Adds sbo / ncbigene / kegg.genes (local) and uniprot (network, optional).
    Idempotent: a gene that already has an 'sbo' annotation is skipped.
    Returns the number of genes annotated.
    """
    import csv
    import os

    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))
    if additions_csv is None:
        additions_csv = os.path.join(root, _GPR_ADDITIONS_CSV)
    if not os.path.exists(additions_csv):
        logger.warning(f"  gene annotate: {additions_csv} not found — skipping")
        return 0

    # YALI1 (no underscore) -> GeneID from NCBI feature table
    y2g = {}
    ft = os.path.join(root, _FEATURE_TABLE)
    if os.path.exists(ft):
        with open(ft) as f:
            next(f)
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 17:
                    continue
                locus, gid = c[16].strip(), c[15].strip()
                if locus and gid:
                    y2g.setdefault(locus.replace("YALI1_", "YALI1"), gid)

    # GeneIDs present in KEGG yli
    kegg_ids = set()
    kg = os.path.join(root, _KEGG_GENES)
    if os.path.exists(kg):
        with open(kg) as f:
            for line in f:
                kegg_ids.add(line.split("\t", 1)[0].replace("yli:", ""))

    genes = sorted({r["add_gene"] for r in csv.DictReader(open(additions_csv))})
    annotated = 0
    for g in genes:
        try:
            gene = model.genes.get_by_id(g)
        except KeyError:
            continue
        ann = dict(gene.annotation) if gene.annotation else {}
        if ann.get("sbo") == "SBO:0000243":
            continue  # idempotent
        ann["sbo"] = "SBO:0000243"
        gid = y2g.get(g)
        if gid:
            ann["ncbigene"] = gid
            if gid in kegg_ids:
                ann["kegg.genes"] = f"yli:{gid}"
            if network:
                up = _fetch_uniprot_for_geneid(gid)
                if up:
                    ann["uniprot"] = up
        gene.annotation = ann
        annotated += 1
        logger.debug(f"  gene annotate: {g} ncbigene={ann.get('ncbigene','-')} "
                     f"uniprot={ann.get('uniprot','-')}")
    return annotated


# ── Patch 10: fill formulas for definite neutral metabolites ──────────────

# A set of metabolites carried no formula but have a definite chemical
# identity.  Their MetaNetX (metanetx.chemical) ids turned out to be empty
# shells with no formula in chem_prop, so the formula was looked up by name in
# PubChem instead (scripts/audit_missing_formula.py ->
# data/missing_formula_fill.csv).
#
# This patch applies ONLY the safe subset: metabolites that are neutral in
# their physiological state (terpenes, esters, amines, nitriles, alcohols,
# aldehydes, peroxide) where charge = 0 is unambiguous.  Charged species
# (carboxylates, CoA-thioesters, phosphate esters, dipeptides) are deliberately
# excluded — the model's own charge convention for those is internally
# inconsistent (e.g. UDP-glucose charge 0 vs GDP-mannose -2), so there is no
# single correct value to fill and getting it wrong would unbalance reactions.
#
# Model formula convention is neutral-H + separate charge (H_model = H_neutral),
# so the PubChem neutral formula is exactly what the model wants.  Idempotent:
# metabolites that already have a formula are skipped.
_NEUTRAL_FILL_CSV = "data/missing_formula_fill.csv"


def fill_neutral_formulas(model, fill_csv: str | None = None) -> int:
    """
    Fill formula (and charge 0) for definite-neutral metabolites listed in the
    fill CSV with status=ready and charge=0, excluding dipeptides.

    Matches by metabolite name; only fills metabolites whose formula is
    currently empty.  Returns the number of metabolite copies filled.
    """
    import csv
    import os
    import re as _re

    if fill_csv is None:
        root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        fill_csv = os.path.join(root, _NEUTRAL_FILL_CSV)
    if not os.path.exists(fill_csv):
        logger.warning(f"  neutral fill: {fill_csv} not found — skipping")
        return 0

    dipeptide = _re.compile(
        r"^(gly|ala|cys|pro|asp|glu|ser|thr|val|leu|ile|phe|tyr|trp|his|lys|"
        r"arg|asn|gln|met)[-_]", _re.I)

    targets = {}  # name -> formula
    with open(fill_csv) as f:
        for row in csv.DictReader(f):
            if row.get("status") != "ready" or row.get("charge") != "0":
                continue
            if dipeptide.match(row["name"]):
                continue
            if row["formula"] and row["formula"] != "(查不到)":
                targets[row["name"]] = row["formula"]

    filled = 0
    for met in model.metabolites:
        if met.formula:
            continue
        f = targets.get(met.name)
        if not f:
            continue
        met.formula = f
        met.charge = 0
        filled += 1
        logger.debug(f"  neutral fill: {met.id} ({met.name}) <- {f} charge=0")
    return filled


# ── Lipid chain-menu extension: add C16:1 palmitoleoyl-CoA to acyl-CoA pools ──

# Y. lipolytica W29 makes ~8% palmitoleate (C16:1) but the iYli21 acyl-CoA pools
# omit it (Carsanba 2020, Table 3: C16:1 = 8.3% on glucose, day 2 —
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7409262/, verified 2026-06).
# We add palmitoleoyl-CoA to the 3 acyl-CoA pools (xPOOL_AC_EM/LP/MM), giving it
# 8.3% of the pool and scaling the existing 6 chains by (1 - 0.083) so the substrate
# weight sum stays 0.951 (= the unchanged product coefficient). palmitoleoyl-CoA
# already exists in all three compartments, so no new metabolite is created.
# Fatty-acid pools (xPOOL_FA_*) are intentionally NOT touched here.
_AC_POOL_C161_FRACTION = 0.083                       # C16:1 share within the pool
_AC_POOLS = ("xPOOL_AC_EM", "xPOOL_AC_LP", "xPOOL_AC_MM")
# palmitoleoyl-CoA id per acyl-CoA-pool compartment (verified present in model)
_PALMITOLEOYL_COA = {"C_em": "m243[C_em]", "C_lp": "m1486[C_lp]", "C_mm": "m1624[C_mm]"}


def extend_acyl_pool_c161(model) -> int:
    """
    Add C16:1 palmitoleoyl-CoA to the 3 acyl-CoA pools and re-scale the existing
    6 chains so the substrate weight sum (= product coefficient, 0.951) is unchanged.

    Idempotent: a pool that already contains its palmitoleoyl-CoA is skipped.
    Returns the number of pools extended. Does not create metabolites and does not
    touch the fatty-acid pools.
    """
    extended = 0
    for rid in _AC_POOLS:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        comp = next(iter({m.compartment for m in rxn.metabolites}))
        c161_id = _PALMITOLEOYL_COA.get(comp)
        if c161_id is None:
            continue
        c161 = model.metabolites.get_by_id(c161_id)
        if c161 in rxn.metabolites:          # already extended → idempotent skip
            continue

        # substrate weight sum (must equal the product coefficient, preserved)
        sub_sum = sum(-c for m, c in rxn.metabolites.items() if c < 0)
        scale = 1.0 - _AC_POOL_C161_FRACTION
        # scale the 6 existing acyl-CoA substrates in place
        delta = {}
        for met, coef in rxn.metabolites.items():
            if coef < 0:
                delta[met] = coef * scale - coef     # bring coef to coef*scale
        rxn.add_metabolites(delta)
        # add palmitoleoyl-CoA at its share of the (unchanged) total
        rxn.add_metabolites({c161: -sub_sum * _AC_POOL_C161_FRACTION})
        extended += 1
        logger.info("  C16:1 pool extension: %s += %s (%.4f), 6 chains x%.3f"
                    % (rid, c161_id, sub_sum * _AC_POOL_C161_FRACTION, scale))
    return extended


# ── Lipid unlump Stage 1: split the single-acyl GPAT/DHAP reactions ──────────

# Stage 1 unlumps the 4 single-acyl acyltransferase reactions (R350/351/352/353)
# into chain-length-specific copies, one per acyl-CoA in the pool. Each copy's
# product formula is back-solved from element conservation (the only unknown).
# An "island" re-mix pool then recombines the chain-specific products back into the
# original generic product, so still-lumped downstream reactions stay connected.
# The original generic reactions are removed (fully replaced by the copies).
# Idempotent: skips if the chain-specific reactions already exist.
_STAGE1_REACTIONS = ("R350", "R351", "R352", "R353")
_REMIX_PREFIX = "xREMIX_"


def _parse_elems(f):
    """Formula -> element dict, or None if generic (empty / '*' / 'R')."""
    if not f or "*" in f or "R" in f:
        return None
    out = {}
    for sym, num in re.findall(r"([A-Z][a-z]?)(\d*)", f):
        if sym:
            out[sym] = out.get(sym, 0) + (int(num) if num else 1)
    return out


def _elems_to_formula(d):
    parts = []
    for sym in ("C", "H"):
        if d.get(sym):
            parts.append(sym + (str(d[sym]) if d[sym] != 1 else ""))
    for sym in sorted(k for k in d if k not in ("C", "H")):
        if d[sym]:
            parts.append(sym + (str(d[sym]) if d[sym] != 1 else ""))
    return "".join(parts)


def _pool_chains(model, generic_acyl_coa):
    """The concrete acyl-CoA substrates (met, weight) that define a generic acyl-CoA,
    read from the xPOOL_AC_* reaction that produces it."""
    for r in model.reactions:
        if not r.id.startswith("xPOOL_AC"):
            continue
        if generic_acyl_coa in r.metabolites and r.metabolites[generic_acyl_coa] > 0:
            return [(mt, -c) for mt, c in r.metabolites.items()
                    if c < 0 and _parse_elems(mt.formula) is not None]
    return []


def unlump_stage1(model) -> int:
    """Split R350/351/352/353 into chain-length-specific reactions (island + re-mix
    pool). Returns the number of generic reactions unlumped (0 if already done)."""
    from cobra import Metabolite, Reaction

    done = 0
    for rid in _STAGE1_REACTIONS:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        if any(r.id.startswith(rid + "_C") for r in model.reactions):
            continue  # idempotent: already unlumped

        # identify generic acyl-CoA substrate and generic product
        gen_coa = next((mt for mt in rxn.metabolites if mt.name == "acyl-CoA_"), None)
        gen_prod = next((mt for mt, c in rxn.metabolites.items()
                         if c > 0 and _parse_elems(mt.formula) is None), None)
        if gen_coa is None or gen_prod is None:
            continue
        chains = _pool_chains(model, gen_coa)
        if not chains:
            continue

        # known element sum of the non-generic, non-acyl-CoA part of the reaction
        comp = gen_prod.compartment
        remix_subs = {}
        new_rxns = []
        for coa_met, weight in chains:
            label = coa_met.name.split("_")[0].replace("-CoA", "")  # e.g. oleoyl
            # back-solve the chain-specific product formula
            known = Counter()
            for mt, c in rxn.metabolites.items():
                if mt is gen_coa:
                    e = _parse_elems(coa_met.formula)
                elif mt is gen_prod:
                    continue  # the unknown
                else:
                    e = _parse_elems(mt.formula)
                if e is None:
                    e = {}
                for s, n in e.items():
                    known[s] += c * n
            ucoef = rxn.metabolites[gen_prod]
            solved = {}
            ok = True
            for s, n in known.items():
                v = -n / ucoef
                if v < -1e-9 or abs(v - round(v)) > 1e-6:
                    ok = False
                    break
                if round(v) != 0:
                    solved[s] = int(round(v))
            if not ok:
                continue

            # new chain-specific product metabolite
            pmet = Metabolite(
                id="%s_%s" % (gen_prod.id.split("[")[0], label) + ("[%s]" % comp),
                name="%s (%s)" % (gen_prod.name, label),
                formula=_elems_to_formula(solved),
                compartment=comp,
            )
            # new chain-specific reaction: copy stoichiometry, swap generic <-> specific
            nr = Reaction(id="%s_%s" % (rid, label), name="%s (%s)" % (rxn.name, label))
            nr.bounds = rxn.bounds                      # inherit bounds (R352 stays 0,0)
            nr.gene_reaction_rule = rxn.gene_reaction_rule
            stoich = {}
            for mt, c in rxn.metabolites.items():
                if mt is gen_coa:
                    stoich[coa_met] = c
                elif mt is gen_prod:
                    stoich[pmet] = c
                else:
                    stoich[mt] = c
            nr.add_metabolites(stoich)
            new_rxns.append(nr)
            remix_subs[pmet] = weight

        if not new_rxns:
            continue
        model.add_reactions(new_rxns)

        # island re-mix pool: weighted chain-specific products -> generic product
        wsum = sum(remix_subs.values())
        remix = Reaction(id="%s%s" % (_REMIX_PREFIX, gen_prod.id.split("[")[0]),
                         name="re-mix pool: %s" % gen_prod.name)
        remix.bounds = (0.0, 1000.0)
        remix.add_metabolites({**{m: -w for m, w in remix_subs.items()}, gen_prod: wsum})
        model.add_reactions([remix])

        # remove the original generic reaction (fully replaced)
        rxn.remove_from_model()
        done += 1
        logger.info("  Stage 1 unlump: %s -> %d chain rxns + re-mix %s"
                    % (rid, len(new_rxns), remix.id))
    return done


# ── CoA / acyl-CoA biochemical-pH protonation curation (inactive) ─────────

# This is deliberately separate from apply_all_patches().  The curated core
# tuples are chemically sound, but the current model's charge/formula convention
# makes their full reaction closure reach several non-CoA acids.  The curation
# file therefore remains activation_state="blocked" until that wider migration
# is reviewed.  Keeping the machinery here provides a strict, reproducible
# preflight rather than permitting an incomplete formula-only patch.
_COA_PROTONATION_CURATION = "data/coa_protonation_curation.json"
_COA_GATE_NON_H_ELEMENTS = frozenset({"C", "N", "O", "P", "S"})
_COA_PROTONATION_SOURCE_MODEL = "../iyali26_gem/model.xml"
_COA_PROTONATION_SOURCE_STAGE = "current_main_model_xml"
_R1521_HANDOFF_CONTRACT = "2f80ce48c930c948c0c0a0d38151b8a40b49f4af7babae586b4b16e04431b5b0"
_SHA256_HEX = re.compile(r"[0-9a-f]{64}\Z")
_COA_CHEMICAL_IDENTITY_KEYS = frozenset(
    {
        "chebi", "hmdb", "inchi", "inchikey", "lipidmaps", "lipidmapsm",
        "metacyc.compound", "metacycm", "metanetx.chemical", "seed.compound", "seedm",
    }
)
_COA_ATOMIC_CLOSURE_ANNOTATIONS = {
    "tetracosanoyl_coa": {"chebi": "CHEBI:65052"},
    "3_oxohexacosanoyl_coa": {"chebi": "CHEBI:73980"},
    "er_vlcfa_3r_hydroxyhexacosanoyl_coa": {"chebi": "CHEBI:76378"},
}


class CoAProtonationCurationError(RuntimeError):
    """The curated CoA protonation patch cannot safely be applied."""


class CoAProtonationActivationBlocked(CoAProtonationCurationError):
    """The current curation intentionally has no production activation."""


def _coa_protonation_curation_path(curation_path: str | Path | None = None) -> Path:
    if curation_path is not None:
        return Path(curation_path)
    return Path(__file__).resolve().parents[2] / _COA_PROTONATION_CURATION


def _coa_protonation_repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _has_source_provenance(curation: dict) -> bool:
    fields = {"source_stage", "source_model_fingerprint"}
    return bool(fields.intersection(curation))


def _validate_curation_source_file(
    curation: dict, source_model_path: str | Path | None = None
) -> Path:
    """Verify that the manifest is tied to this exact post-annotation snapshot."""

    source_path = (
        Path(source_model_path).resolve()
        if source_model_path is not None
        else (_coa_protonation_repository_root() / curation["source_model"]).resolve()
    )
    if not source_path.is_file():
        raise CoAProtonationCurationError(
            f"CoA protonation source model is unavailable: {source_path}"
        )
    actual_sha256 = _sha256_file(source_path)
    if actual_sha256 != curation["source_sha256"]:
        raise CoAProtonationCurationError(
            "CoA protonation source model SHA-256 drifted from the curated "
            f"post-annotation snapshot ({curation['source_sha256']} -> {actual_sha256})"
        )
    return source_path


def _validate_r1521_dependency(curation: dict) -> None:
    """Keep the merged tuples tied to the separate, source-reviewed handoff."""

    if curation.get("source_model") != _COA_PROTONATION_SOURCE_MODEL:
        return
    if curation.get("r1521_current_snapshot_contract_sha256") != _R1521_HANDOFF_CONTRACT:
        raise CoAProtonationCurationError("R1521 handoff contract digest drifted")
    path = _coa_protonation_repository_root() / "data/r1521_current_snapshot_handoff.json"
    try:
        handoff = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CoAProtonationCurationError(f"cannot read R1521 handoff: {exc}") from exc
    if (
        handoff.get("target_contract_sha256") != _R1521_HANDOFF_CONTRACT
        or handoff.get("source_sha256") != curation.get("source_sha256")
        or handoff.get("source_model_fingerprint") != curation.get("source_model_fingerprint")
    ):
        raise CoAProtonationCurationError("R1521 handoff source contract drifted")
    targets = {
        entry["id"]: (entry["formula"], entry["charge"])
        for entry in handoff["local_counterfactual"]["target_tuples"]
    }
    curated = {
        met_id: _tuple_from_record(group["target_tuple"])
        for group in curation["groups"]
        for met_id in group["expected_ids"]
    }
    if any(curated.get(met_id) != target for met_id, target in targets.items()):
        raise CoAProtonationCurationError("R1521 tuple contract drifted")


def load_coa_protonation_curation(
    curation_path: str | Path | None = None,
    *,
    source_model_path: str | Path | None = None,
) -> dict:
    """Load and structurally validate the curated CoA protonation manifest."""
    path = _coa_protonation_curation_path(curation_path)
    try:
        with path.open(encoding="utf-8") as handle:
            curation = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise CoAProtonationCurationError(
            f"Cannot read CoA protonation curation {path}: {exc}"
        ) from exc
    _validate_coa_protonation_manifest(
        curation, require_source_provenance=True, require_source_file=True
    )
    _validate_curation_source_file(curation, source_model_path)
    return curation


def _validate_coa_protonation_manifest(
    curation: dict,
    *,
    require_source_provenance: bool = False,
    require_source_file: bool = False,
) -> None:
    """Reject a malformed or internally ambiguous curation manifest."""
    if not isinstance(curation, dict):
        raise CoAProtonationCurationError("CoA protonation curation must be an object")
    if type(curation.get("schema_version")) is not int or curation.get(
        "schema_version"
    ) != 1:
        raise CoAProtonationCurationError("Unsupported CoA protonation curation schema")
    if curation.get("activation_state") not in {"approved", "blocked"}:
        raise CoAProtonationCurationError(
            "CoA protonation activation_state must be 'approved' or 'blocked'"
        )

    provenance_present = _has_source_provenance(curation)
    if require_source_provenance and not provenance_present:
        raise CoAProtonationCurationError(
            "CoA protonation curation lacks required source snapshot provenance"
        )
    if provenance_present:
        required_provenance = {"source_stage", "source_model_fingerprint"}
        if not required_provenance.issubset(curation):
            raise CoAProtonationCurationError(
                "CoA protonation source provenance is incomplete"
            )
        if type(curation["source_stage"]) is not str or not curation["source_stage"]:
            raise CoAProtonationCurationError(
                "CoA protonation source_stage must be a non-empty string"
            )
        if (
            ("source_model" in curation) != ("source_sha256" in curation)
        ):
            raise CoAProtonationCurationError(
                "CoA protonation source file provenance is incomplete"
            )
        has_source_file = "source_model" in curation
        if require_source_file and not has_source_file:
            raise CoAProtonationCurationError(
                "CoA protonation curation lacks required source model file provenance"
            )
        if has_source_file:
            if curation["source_model"] != _COA_PROTONATION_SOURCE_MODEL:
                raise CoAProtonationCurationError(
                    "CoA protonation curation must reference the current model.xml snapshot"
                )
            if curation["source_stage"] != _COA_PROTONATION_SOURCE_STAGE:
                raise CoAProtonationCurationError(
                    "CoA protonation curation has an unexpected source model stage"
                )
        for field in ("source_model_fingerprint", "target_model_fingerprint"):
            if field not in curation:
                continue
            value = curation[field]
            if type(value) is not str or _SHA256_HEX.fullmatch(value) is None:
                raise CoAProtonationCurationError(
                    f"CoA protonation {field} must be a lowercase SHA-256 hex digest"
                )
        if has_source_file:
            value = curation["source_sha256"]
            if type(value) is not str or _SHA256_HEX.fullmatch(value) is None:
                raise CoAProtonationCurationError(
                    "CoA protonation source_sha256 must be a lowercase SHA-256 hex digest"
                )
    elif require_source_file:
        raise CoAProtonationCurationError(
            "CoA protonation curation lacks required source model file provenance"
        )
    if (
        curation["activation_state"] == "approved"
        and "target_model_fingerprint" not in curation
    ):
        raise CoAProtonationCurationError(
            "Approved CoA protonation curation requires target_model_fingerprint"
        )

    activation_blockers = curation.get("activation_blockers")
    if type(activation_blockers) is not list:
        raise CoAProtonationCurationError(
            "CoA protonation curation must declare activation_blockers as a list"
        )
    blocker_reaction_ids: set[str] = set()
    for blocker in activation_blockers:
        if type(blocker) is not dict or not isinstance(
            blocker.get("reaction_id"), str
        ) or not blocker["reaction_id"]:
            raise CoAProtonationCurationError(
                "CoA protonation activation blockers need non-empty reaction ids"
            )
        blocker_reaction_ids.add(blocker["reaction_id"])
    if len(blocker_reaction_ids) != len(activation_blockers):
        raise CoAProtonationCurationError(
            "CoA protonation activation blocker reaction ids must be unique"
        )
    target_balance_reaction_ids = curation.get("target_balance_reaction_ids")
    if (
        type(target_balance_reaction_ids) is not list
        or not target_balance_reaction_ids
        or len(set(target_balance_reaction_ids)) != len(target_balance_reaction_ids)
        or not all(isinstance(reaction_id, str) and reaction_id for reaction_id in target_balance_reaction_ids)
    ):
        raise CoAProtonationCurationError(
            "CoA protonation curation needs unique target_balance_reaction_ids"
        )
    if not blocker_reaction_ids.issubset(target_balance_reaction_ids):
        raise CoAProtonationCurationError(
            "Every CoA protonation activation blocker must be target-balance checked"
        )

    groups = curation.get("groups")
    if not isinstance(groups, list) or not groups:
        raise CoAProtonationCurationError("CoA protonation curation has no groups")

    seen_group_ids: set[str] = set()
    seen_metabolite_ids: set[str] = set()
    for group in groups:
        if not isinstance(group, dict):
            raise CoAProtonationCurationError("CoA protonation group must be an object")
        group_id = group.get("id")
        identity_names = group.get("identity_names")
        expected_ids = group.get("expected_ids")
        target = group.get("target_tuple")
        if not isinstance(group_id, str) or not group_id or group_id in seen_group_ids:
            raise CoAProtonationCurationError("CoA protonation group ids must be unique")
        if not isinstance(identity_names, list) or not all(
            isinstance(name, str) and name for name in identity_names
        ):
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} lacks identity_names"
            )
        if (
            not isinstance(expected_ids, list)
            or not expected_ids
            or len(set(expected_ids)) != len(expected_ids)
            or not all(isinstance(met_id, str) and met_id for met_id in expected_ids)
        ):
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} has invalid expected_ids"
            )
        if type(group.get("expected_copy_count")) is not int or group.get(
            "expected_copy_count"
        ) != len(expected_ids):
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} copy count disagrees with ids"
            )
        if (
            not isinstance(target, dict)
            or not isinstance(target.get("formula"), str)
            or type(target.get("charge")) is not int
        ):
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} lacks an exact target tuple"
            )
        target_name = group.get("target_name")
        if target_name is not None and (
            not isinstance(target_name, str) or not target_name
        ):
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} has an invalid target_name"
            )

        legacy_ids: list[str] = []
        for legacy in group.get("legacy_tuples", []):
            if (
                not isinstance(legacy, dict)
                or not (
                    legacy.get("formula") is None
                    or isinstance(legacy.get("formula"), str)
                )
                or type(legacy.get("charge")) is not int
            ):
                raise CoAProtonationCurationError(
                    f"CoA protonation group {group_id} has an invalid legacy tuple"
                )
            ids = legacy.get("ids")
            if not isinstance(ids, list) or not ids or not all(isinstance(met_id, str) for met_id in ids):
                raise CoAProtonationCurationError(
                    f"CoA protonation group {group_id} has an invalid legacy id list"
                )
            legacy_ids.extend(ids)
        if len(set(legacy_ids)) != len(legacy_ids) or not set(legacy_ids) <= set(expected_ids):
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} has ambiguous legacy ids"
            )
        if set(legacy_ids) != set(expected_ids):
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} source tuples do not cover every copy"
            )
        if group_id in _COA_ATOMIC_CLOSURE_ANNOTATIONS:
            if group.get("annotation_target") != _COA_ATOMIC_CLOSURE_ANNOTATIONS[group_id]:
                raise CoAProtonationCurationError(
                    f"CoA protonation group {group_id} has an invalid chemical annotation target"
                )
            replacement_keys = group.get("replace_annotation_keys")
            if (
                not isinstance(replacement_keys, list)
                or set(replacement_keys) != _COA_CHEMICAL_IDENTITY_KEYS
                or len(replacement_keys) != len(_COA_CHEMICAL_IDENTITY_KEYS)
            ):
                raise CoAProtonationCurationError(
                    f"CoA protonation group {group_id} has incomplete chemical annotation replacement keys"
                )
        if set(expected_ids) & seen_metabolite_ids:
            raise CoAProtonationCurationError(
                f"CoA protonation group {group_id} reuses a metabolite id"
            )
        seen_group_ids.add(group_id)
        seen_metabolite_ids.update(expected_ids)

    _validate_r1521_dependency(curation)

    reaction_corrections = curation.get("reaction_corrections", [])
    if type(reaction_corrections) is not list:
        raise CoAProtonationCurationError(
            "CoA protonation reaction_corrections must be a list"
        )
    correction_reaction_ids: set[str] = set()
    for correction in reaction_corrections:
        if not isinstance(correction, dict) or not all(
            isinstance(correction.get(key), str) and correction[key]
            for key in ("reaction_id", "metabolite_id")
        ):
            raise CoAProtonationCurationError("Invalid curated CoA reaction correction")
        for field in ("legacy_coefficient", "target_coefficient"):
            _finite_numeric(
                correction.get(field), f"CoA reaction correction {field}"
            )
        correction_reaction_ids.add(correction["reaction_id"])
    if not correction_reaction_ids.issubset(target_balance_reaction_ids):
        raise CoAProtonationCurationError(
            "Every CoA protonation reaction correction must be target-balance checked"
        )

    reaction_contracts = curation.get("reaction_contracts", [])
    if type(reaction_contracts) is not list:
        raise CoAProtonationCurationError(
            "CoA protonation reaction_contracts must be a list"
        )
    seen_contract_reactions: set[str] = set()
    for contract in reaction_contracts:
        if not isinstance(contract, dict) or set(contract) != {
            "reaction_id", "bounds", "reversible", "source_stoichiometry",
            "target_stoichiometry", "evidence",
        }:
            raise CoAProtonationCurationError("Invalid curated CoA reaction contract")
        reaction_id = contract["reaction_id"]
        if not isinstance(reaction_id, str) or not reaction_id or reaction_id in seen_contract_reactions:
            raise CoAProtonationCurationError(
                "CoA protonation reaction contract ids must be unique"
            )
        bounds = contract["bounds"]
        if (
            not isinstance(bounds, list)
            or len(bounds) != 2
            or any(isinstance(value, bool) or not isinstance(value, Real) for value in bounds)
        ):
            raise CoAProtonationCurationError(
                f"CoA reaction contract {reaction_id} has invalid bounds"
            )
        if type(contract["reversible"]) is not bool:
            raise CoAProtonationCurationError(
                f"CoA reaction contract {reaction_id} has invalid reversibility"
            )
        for field in ("source_stoichiometry", "target_stoichiometry"):
            stoichiometry = contract[field]
            if (
                not isinstance(stoichiometry, dict)
                or not stoichiometry
                or not all(
                    isinstance(metabolite_id, str)
                    and metabolite_id
                    and not isinstance(coefficient, bool)
                    and isinstance(coefficient, Real)
                    and isfinite(float(coefficient))
                    and float(coefficient) != 0.0
                    for metabolite_id, coefficient in stoichiometry.items()
                )
            ):
                raise CoAProtonationCurationError(
                    f"CoA reaction contract {reaction_id} has invalid {field}"
                )
        if not isinstance(contract["evidence"], str) or not contract["evidence"]:
            raise CoAProtonationCurationError(
                f"CoA reaction contract {reaction_id} lacks evidence"
            )
        seen_contract_reactions.add(reaction_id)
    if not seen_contract_reactions.issubset(target_balance_reaction_ids):
        raise CoAProtonationCurationError(
            "Every CoA protonation reaction contract must be target-balance checked"
        )

    identity_corrections = curation.get("reaction_identity_corrections", [])
    if not isinstance(identity_corrections, list):
        raise CoAProtonationCurationError(
            "CoA protonation reaction_identity_corrections must be a list"
        )
    seen_identity_reactions: set[str] = set()
    for correction in identity_corrections:
        if not isinstance(correction, dict) or not all(
            isinstance(correction.get(field), str) and correction[field]
            for field in ("reaction_id", "legacy_name", "target_name")
        ):
            raise CoAProtonationCurationError(
                "Invalid curated CoA reaction identity correction"
            )
        reaction_id = correction["reaction_id"]
        if reaction_id in seen_identity_reactions:
            raise CoAProtonationCurationError(
                "CoA protonation reaction identity corrections must be unique"
            )
        annotation_target = correction.get("annotation_target", {})
        if not isinstance(annotation_target, dict) or not all(
            isinstance(key, str) and key for key in annotation_target
        ):
            raise CoAProtonationCurationError(
                "CoA reaction identity annotations must be an object with string keys"
            )
        seen_identity_reactions.add(reaction_id)


def _coa_identity_key(name: str | None) -> str:
    """Normalize only the iYli21 trailing formula token for exact identity checks."""
    normalized = (name or "").strip().lower().rstrip("_")
    return re.sub(r"_[A-Z][A-Za-z0-9]*$", "", normalized, flags=re.IGNORECASE)


def _exact_live_charge(value, metabolite_id: str) -> int | None:
    """Return an integral live charge without silently coercing model state."""

    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real):
        raise CoAProtonationCurationError(
            f"{metabolite_id} has a non-numeric or boolean charge {value!r}"
        )
    numeric = float(value)
    if not isfinite(numeric) or not numeric.is_integer():
        raise CoAProtonationCurationError(
            f"{metabolite_id} has a non-finite or fractional charge {value!r}"
        )
    return int(value)


def _finite_numeric(value, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise CoAProtonationCurationError(f"{label} must be a finite numeric value")
    numeric = float(value)
    if not isfinite(numeric):
        raise CoAProtonationCurationError(f"{label} must be a finite numeric value")
    return numeric


def _tuple_of(met) -> tuple[str | None, int | None]:
    return met.formula, _exact_live_charge(met.charge, met.id)


def _tuple_from_record(record: dict) -> tuple[str, int]:
    charge = record["charge"]
    if type(charge) is not int:
        raise CoAProtonationCurationError("Curated CoA tuple charge must be an int")
    return record["formula"], charge


def _group_source_tuple_fingerprint(group: dict) -> str:
    """Hash the exact pre-migration tuple expected for every named copy."""

    legacy_by_id = {
        metabolite_id: {"formula": legacy["formula"], "charge": legacy["charge"]}
        for legacy in group.get("legacy_tuples", [])
        for metabolite_id in legacy["ids"]
    }
    payload = {
        "ids": sorted(group["expected_ids"]),
        "legacy_tuples": {metabolite_id: legacy_by_id[metabolite_id] for metabolite_id in sorted(legacy_by_id)},
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def _legacy_tuples_by_id(group: dict) -> dict[str, tuple[str, int]]:
    legacy_by_id: dict[str, tuple[str, int]] = {}
    for legacy in group.get("legacy_tuples", []):
        tuple_value = _tuple_from_record(legacy)
        for met_id in legacy["ids"]:
            legacy_by_id[met_id] = tuple_value
    return legacy_by_id


def _reaction_stoichiometry(reaction) -> dict[str, float]:
    return {
        metabolite.id: float(coefficient)
        for metabolite, coefficient in reaction.metabolites.items()
    }


def _reaction_matches_contract(reaction, contract: dict, field: str) -> bool:
    return (
        isclose(float(reaction.lower_bound), float(contract["bounds"][0]), abs_tol=1e-12)
        and isclose(float(reaction.upper_bound), float(contract["bounds"][1]), abs_tol=1e-12)
        and bool(reaction.reversibility) is contract["reversible"]
        and _reaction_stoichiometry(reaction)
        == {metabolite_id: float(coefficient) for metabolite_id, coefficient in contract[field].items()}
    )


def _model_snapshot_fingerprint(model) -> str:
    """Hash the complete in-memory model state needed to identify its source stage."""

    metabolites = []
    for met in sorted(model.metabolites, key=lambda candidate: candidate.id):
        annotation = getattr(met, "annotation", {}) or {}
        if type(annotation) is not dict:
            raise CoAProtonationCurationError(
                f"{met.id} has a non-object annotation in source snapshot"
            )
        metabolites.append(
            {
                "annotation": annotation,
                "charge": _exact_live_charge(met.charge, met.id),
                "compartment": met.compartment,
                "formula": met.formula,
                "id": met.id,
                "name": met.name,
            }
        )

    reactions = []
    objective_coefficients: list[dict[str, float | str]] = []
    for reaction in sorted(model.reactions, key=lambda candidate: candidate.id):
        coefficient = _finite_numeric(
            reaction.objective_coefficient,
            f"{reaction.id} objective coefficient",
        )
        if coefficient != 0.0:
            objective_coefficients.append({"id": reaction.id, "value": coefficient})
        reactions.append(
            {
                "gene_reaction_rule": reaction.gene_reaction_rule,
                "id": reaction.id,
                "lower_bound": _finite_numeric(
                    reaction.lower_bound, f"{reaction.id} lower bound"
                ),
                "metabolites": [
                    {
                        "id": met.id,
                        "coefficient": _finite_numeric(
                            value, f"{reaction.id}/{met.id} coefficient"
                        ),
                    }
                    for met, value in sorted(
                        reaction.metabolites.items(), key=lambda item: item[0].id
                    )
                ],
                "upper_bound": _finite_numeric(
                    reaction.upper_bound, f"{reaction.id} upper bound"
                ),
            }
        )
    payload = {
        "id": model.id,
        "metabolites": metabolites,
        "objective": {
            "direction": str(model.objective.direction),
            "coefficients": objective_coefficients,
        },
        "reactions": reactions,
    }
    try:
        serialized = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise CoAProtonationCurationError(
            "source model has a non-canonical value that cannot be fingerprinted"
        ) from exc
    return hashlib.sha256(serialized.encode("utf-8")).hexdigest()


def _source_snapshot_report(
    model, curation: dict, *, source_model_path: str | Path | None = None
) -> dict:
    """Verify both the declared source file and the supplied in-memory model."""

    if not _has_source_provenance(curation):
        return {
            "declared": False,
            "file_sha256_verified": None,
            "model_fingerprint_verified": None,
            "verified": False,
            "errors": ["CoA protonation curation lacks source provenance"],
        }

    errors: list[str] = []
    has_source_file = "source_model" in curation
    file_sha256_verified: bool | None = None
    model_fingerprint_verified = False
    target_model_fingerprint_verified: bool | None = None
    actual_fingerprint: str | None = None
    if has_source_file:
        file_sha256_verified = False
        try:
            _validate_curation_source_file(curation, source_model_path)
            file_sha256_verified = True
        except CoAProtonationCurationError as exc:
            errors.append(str(exc))
    try:
        actual_fingerprint = _model_snapshot_fingerprint(model)
        model_fingerprint_verified = (
            actual_fingerprint == curation["source_model_fingerprint"]
        )
        if "target_model_fingerprint" in curation:
            target_model_fingerprint_verified = (
                actual_fingerprint == curation["target_model_fingerprint"]
            )
        if not model_fingerprint_verified and not target_model_fingerprint_verified:
            errors.append(
                "CoA protonation model fingerprint differs from both the curated "
                "post-annotation source and declared target snapshots"
            )
    except CoAProtonationCurationError as exc:
        errors.append(str(exc))
    return {
        "declared": True,
        "source_model": curation.get("source_model"),
        "source_stage": curation["source_stage"],
        "source_sha256": curation.get("source_sha256"),
        "source_model_fingerprint": curation["source_model_fingerprint"],
        "target_model_fingerprint": curation.get("target_model_fingerprint"),
        "actual_model_fingerprint": actual_fingerprint,
        "file_sha256_verified": file_sha256_verified,
        "model_fingerprint_verified": model_fingerprint_verified,
        "target_model_fingerprint_verified": target_model_fingerprint_verified,
        "verified": (
            (model_fingerprint_verified or bool(target_model_fingerprint_verified))
            and (file_sha256_verified is None or file_sha256_verified)
        ),
        "errors": errors,
    }


def validate_coa_protonation_curation(
    model,
    curation: dict | None = None,
    *,
    source_model_path: str | Path | None = None,
) -> dict:
    """
    Check the live model against the exact curated IDs and tuple states.

    This function never mutates ``model``.  Each expected copy may be at its
    explicitly enumerated legacy tuple or target tuple; any third tuple, a
    missing copy, a new copy with the same identity, or a changed reaction
    coefficient is an error.
    """
    if curation is None:
        curation = load_coa_protonation_curation(
            source_model_path=source_model_path
        )
    else:
        _validate_coa_protonation_manifest(curation, require_source_provenance=True)

    errors: list[str] = []
    for group in curation["groups"]:
        group_id = group["id"]
        expected_ids = set(group["expected_ids"])
        identity_names = set(group["identity_names"])
        observed_ids = {
            met.id
            for met in model.metabolites
            if _coa_identity_key(met.name) in identity_names
        }
        missing = sorted(expected_ids - observed_ids)
        unexpected = sorted(observed_ids - expected_ids)
        if missing:
            errors.append(f"{group_id}: missing expected copy/copies {missing}")
        if unexpected:
            errors.append(f"{group_id}: unexpected same-identity copy/copies {unexpected}")
        if len(observed_ids) != group["expected_copy_count"]:
            errors.append(
                f"{group_id}: observed {len(observed_ids)} copies, expected "
                f"{group['expected_copy_count']}"
            )

        target = _tuple_from_record(group["target_tuple"])
        legacy_by_id = _legacy_tuples_by_id(group)
        for met_id in sorted(expected_ids):
            try:
                met = model.metabolites.get_by_id(met_id)
            except KeyError:
                continue
            actual = _tuple_of(met)
            allowed = {target}
            if met_id in legacy_by_id:
                allowed.add(legacy_by_id[met_id])
            if actual not in allowed:
                errors.append(
                    f"{group_id}: {met_id} has third tuple {actual!r}; "
                    f"allowed={sorted(allowed)!r}"
                )

            target_name = group.get("target_name")
            if (
                target_name is not None
                and met.name != target_name
                and _coa_identity_key(met.name) not in set(group["identity_names"])
            ):
                errors.append(
                    f"{group_id}: {met_id} has unexpected identity {met.name!r}"
                )

    for correction in curation.get("reaction_corrections", []):
        rid = correction["reaction_id"]
        met_id = correction["metabolite_id"]
        try:
            reaction = model.reactions.get_by_id(rid)
            met = model.metabolites.get_by_id(met_id)
        except KeyError as exc:
            errors.append(f"reaction correction {rid}: missing {exc.args[0]}")
            continue
        coefficient = reaction.metabolites.get(met, 0.0)
        allowed = {
            float(correction["legacy_coefficient"]),
            float(correction["target_coefficient"]),
        }
        if coefficient not in allowed:
            errors.append(
                f"reaction correction {rid}: {met_id} coefficient {coefficient} "
                f"is not an exact curated state {sorted(allowed)}"
            )

    for contract in curation.get("reaction_contracts", []):
        try:
            reaction = model.reactions.get_by_id(contract["reaction_id"])
        except KeyError:
            errors.append(f"reaction contract {contract['reaction_id']}: missing reaction")
            continue
        if not any(
            _reaction_matches_contract(reaction, contract, field)
            for field in ("source_stoichiometry", "target_stoichiometry")
        ):
            errors.append(
                f"reaction contract {reaction.id}: bounds, direction, or stoichiometry drifted"
            )

    for correction in curation.get("reaction_identity_corrections", []):
        rid = correction["reaction_id"]
        try:
            reaction = model.reactions.get_by_id(rid)
        except KeyError:
            errors.append(f"reaction identity correction {rid}: missing reaction")
            continue
        if reaction.name not in {correction["legacy_name"], correction["target_name"]}:
            errors.append(
                f"reaction identity correction {rid}: unexpected identity "
                f"{reaction.name!r}"
            )

    return {"valid": not errors, "errors": errors}


def _safe_mass_balance(reaction) -> dict | None:
    """Return a mass-balance residual, or None when COBRA cannot parse it."""
    try:
        return reaction.check_mass_balance()
    except (TypeError, ValueError):
        return None


def _formula_complete_balanced_reactions(model) -> set[str]:
    return {
        reaction.id
        for reaction in model.reactions
        if all(met.formula for met in reaction.metabolites)
        and _safe_mass_balance(reaction) == {}
    }


def _all_mass_balances(model) -> dict[str, dict | None]:
    return {reaction.id: _safe_mass_balance(reaction) for reaction in model.reactions}


def _target_reaction_residuals(model, curation: dict) -> dict[str, dict | None]:
    """Require every explicitly declared target reaction to be exactly balanced."""

    residuals: dict[str, dict | None] = {}
    for reaction_id in curation["target_balance_reaction_ids"]:
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            residuals[reaction_id] = None
            continue
        balance = _safe_mass_balance(reaction)
        if balance != {}:
            residuals[reaction_id] = balance
    return residuals


def _balance_delta(before: dict | None, after: dict | None) -> dict[str, float]:
    if before is None or after is None:
        return {}
    keys = set(before) | set(after)
    return {
        key: float(after.get(key, 0.0) - before.get(key, 0.0))
        for key in keys
        if not isclose(float(after.get(key, 0.0)), float(before.get(key, 0.0)), abs_tol=1e-12)
    }


def _objective_value(model) -> float | None:
    """Return a finite FBA optimum when one is available; otherwise None."""
    try:
        value = model.slim_optimize(error_value=None)
    except Exception:  # solver availability is outside patch correctness
        return None
    if value is None or not isfinite(float(value)):
        return None
    return float(value)


def _apply_coa_protonation_curation_unsafe(model, curation: dict) -> dict[str, int]:
    """Apply a preflighted curation to an in-memory model only."""
    counts = {
        "metabolites": 0,
        "annotations": 0,
        "reaction_corrections": 0,
        "identity_corrections": 0,
    }
    for group in curation["groups"]:
        target = _tuple_from_record(group["target_tuple"])
        legacy_by_id = _legacy_tuples_by_id(group)
        for met_id in group["expected_ids"]:
            met = model.metabolites.get_by_id(met_id)
            actual = _tuple_of(met)
            if actual != target:
                if legacy_by_id.get(met_id) != actual:
                    raise CoAProtonationCurationError(
                        f"Refusing third tuple for {met_id}: {actual!r}"
                    )
                met.formula, met.charge = target
                counts["metabolites"] += 1

            target_name = group.get("target_name")
            if target_name is not None and met.name != target_name:
                met.name = target_name
                counts["identity_corrections"] += 1

            annotation_target = group.get("annotation_target")
            if annotation_target:
                annotation = dict(met.annotation or {})
                before = dict(annotation)
                for key in group.get("replace_annotation_keys", []):
                    annotation.pop(key, None)
                annotation.update(annotation_target)
                if annotation != before:
                    met.annotation = annotation
                    counts["annotations"] += 1

    for correction in curation.get("reaction_corrections", []):
        reaction = model.reactions.get_by_id(correction["reaction_id"])
        met = model.metabolites.get_by_id(correction["metabolite_id"])
        legacy = float(correction["legacy_coefficient"])
        target = float(correction["target_coefficient"])
        actual = float(reaction.metabolites.get(met, 0.0))
        if actual == target:
            continue
        if actual != legacy:
            raise CoAProtonationCurationError(
                f"Refusing third coefficient for {reaction.id}/{met.id}: {actual}"
            )
        reaction.add_metabolites({met: target - actual})
        counts["reaction_corrections"] += 1

    for correction in curation.get("reaction_identity_corrections", []):
        reaction = model.reactions.get_by_id(correction["reaction_id"])
        if reaction.name != correction["target_name"]:
            if reaction.name != correction["legacy_name"]:
                raise CoAProtonationCurationError(
                    f"Refusing unexpected reaction identity for {reaction.id}: "
                    f"{reaction.name!r}"
                )
            reaction.name = correction["target_name"]
            counts["identity_corrections"] += 1
        annotation_target = correction.get("annotation_target", {})
        if annotation_target:
            annotation = dict(reaction.annotation or {})
            before = dict(annotation)
            annotation.update(annotation_target)
            if annotation != before:
                reaction.annotation = annotation
                counts["annotations"] += 1

    for contract in curation.get("reaction_contracts", []):
        reaction = model.reactions.get_by_id(contract["reaction_id"])
        if not _reaction_matches_contract(reaction, contract, "target_stoichiometry"):
            raise CoAProtonationCurationError(
                f"Refusing incomplete CoA reaction contract for {reaction.id}"
            )
    return counts


def _apply_coa_protonation_curation(model, curation: dict) -> dict[str, int]:
    """Apply the staged curation transactionally after a caller's preflight."""

    metabolite_ids = {
        metabolite_id
        for group in curation["groups"]
        for metabolite_id in group["expected_ids"]
    }
    reaction_ids = {
        correction["reaction_id"]
        for correction in curation.get("reaction_corrections", [])
    } | {
        correction["reaction_id"]
        for correction in curation.get("reaction_identity_corrections", [])
    } | {
        contract["reaction_id"]
        for contract in curation.get("reaction_contracts", [])
    }
    metabolites = [model.metabolites.get_by_id(metabolite_id) for metabolite_id in metabolite_ids]
    reactions = [model.reactions.get_by_id(reaction_id) for reaction_id in reaction_ids]
    metabolite_snapshot = [
        (metabolite, metabolite.name, metabolite.formula, metabolite.charge, deepcopy(metabolite.annotation or {}))
        for metabolite in metabolites
    ]
    reaction_snapshot = [
        (reaction, reaction.name, deepcopy(reaction.annotation or {}), dict(reaction.metabolites))
        for reaction in reactions
    ]
    try:
        return _apply_coa_protonation_curation_unsafe(model, curation)
    except Exception:
        for metabolite, name, formula, charge, annotation in metabolite_snapshot:
            metabolite.name, metabolite.formula, metabolite.charge, metabolite.annotation = (
                name, formula, charge, annotation,
            )
        for reaction, name, annotation, stoichiometry in reaction_snapshot:
            current = dict(reaction.metabolites)
            changes = {
                metabolite: stoichiometry.get(metabolite, 0.0) - current.get(metabolite, 0.0)
                for metabolite in set(current) | set(stoichiometry)
                if not isclose(
                    float(stoichiometry.get(metabolite, 0.0)),
                    float(current.get(metabolite, 0.0)),
                    abs_tol=1e-12,
                )
            }
            if changes:
                reaction.add_metabolites(changes)
            reaction.name, reaction.annotation = name, annotation
        raise


def _evaluate_coa_protonation_gate(
    model, curation: dict, *, source_model_path: str | Path | None = None
) -> dict:
    """Evaluate all source, target-balance, and objective gates on a model copy."""

    validation = validate_coa_protonation_curation(
        model, curation, source_model_path=source_model_path
    )
    source_snapshot = _source_snapshot_report(
        model, curation, source_model_path=source_model_path
    )
    activation_blockers_clear = not curation["activation_blockers"]
    report = {
        "validation": validation,
        "source_snapshot": source_snapshot,
        "activation_blockers_clear": activation_blockers_clear,
    }
    if not validation["valid"]:
        report.update(
            {
                "gate_passed": False,
                "newly_unbalanced": {},
                "non_h_charge_deltas": {},
                "target_reactions_unbalanced": {},
                "balanced_count_before": None,
                "balanced_count_after": None,
                "objective_before": None,
                "objective_after": None,
                "objective_unchanged": False,
            }
        )
        return report

    candidate = model.copy()
    balanced_before = _formula_complete_balanced_reactions(candidate)
    formula_complete_before = {
        reaction.id
        for reaction in candidate.reactions
        if all(met.formula for met in reaction.metabolites)
    }
    balances_before = _all_mass_balances(candidate)
    objective_before = _objective_value(candidate)
    counts = _apply_coa_protonation_curation(candidate, curation)
    balances_after = _all_mass_balances(candidate)
    balanced_after = _formula_complete_balanced_reactions(candidate)
    objective_after = _objective_value(candidate)
    target_reactions_unbalanced = _target_reaction_residuals(candidate, curation)

    newly_unbalanced = {
        reaction_id: balances_after[reaction_id]
        for reaction_id in sorted(balanced_before)
        if balances_after[reaction_id] not in ({}, None)
    }
    changed_residuals = {
        reaction_id: _balance_delta(balances_before[reaction_id], balances_after[reaction_id])
        for reaction_id in balances_before
        if _balance_delta(balances_before[reaction_id], balances_after[reaction_id])
    }
    non_h_charge_deltas = {
        reaction_id: delta
        for reaction_id, delta in changed_residuals.items()
        if (
            reaction_id in formula_complete_before
            and _COA_GATE_NON_H_ELEMENTS.intersection(delta)
        )
    }
    objective_unchanged = (
        objective_before is not None
        and objective_after is not None
        and isclose(objective_before, objective_after, rel_tol=0.0, abs_tol=1e-9)
    )
    report.update(
        {
            "counts": counts,
            "newly_unbalanced": newly_unbalanced,
            "non_h_charge_deltas": non_h_charge_deltas,
            "target_reactions_unbalanced": target_reactions_unbalanced,
            "balanced_count_before": len(balanced_before),
            "balanced_count_after": len(balanced_after),
            "objective_before": None if objective_before is None else round(objective_before, 12),
            "objective_after": None if objective_after is None else round(objective_after, 12),
            "objective_unchanged": objective_unchanged,
        }
    )
    report["gate_passed"] = (
        source_snapshot["verified"]
        and activation_blockers_clear
        and not newly_unbalanced
        and not non_h_charge_deltas
        and not target_reactions_unbalanced
        and len(balanced_after) >= len(balanced_before)
        and objective_unchanged
    )
    return report


def audit_coa_protonation_curation(
    model,
    curation: dict | None = None,
    curation_path: str | Path | None = None,
    *,
    source_model_path: str | Path | None = None,
) -> dict:
    """
    Read-only audit of the curated CoA protonation proposal.

    It validates exact identity/copy/tuple states and evaluates all safety
    gates on a copy.  ``activation_blockers`` are reported verbatim so callers
    cannot mistake a chemically incomplete closure for an applicable patch.
    """
    if curation is None:
        curation = load_coa_protonation_curation(
            curation_path, source_model_path=source_model_path
        )
    else:
        _validate_coa_protonation_manifest(curation, require_source_provenance=True)
    report = _evaluate_coa_protonation_gate(
        model, curation, source_model_path=source_model_path
    )
    report["activation_state"] = curation["activation_state"]
    report["activation_reason"] = curation.get("activation_reason", "")
    report["activation_blockers"] = curation["activation_blockers"]
    report["ready_for_activation"] = (
        curation["activation_state"] == "approved"
        and report["activation_blockers_clear"]
        and report["gate_passed"]
    )
    return report


def normalize_coa_protonation(
    model,
    curation: dict | None = None,
    curation_path: str | Path | None = None,
    *,
    source_model_path: str | Path | None = None,
) -> dict[str, int]:
    """
    Apply the CoA curation only after every exact preflight gate passes.

    The repository curation is intentionally blocked and this function raises
    before mutating the model.  It is not wired into apply_all_patches() or the
    pipeline until the documented external protonation closure is approved.
    """
    if curation is None:
        curation = load_coa_protonation_curation(
            curation_path, source_model_path=source_model_path
        )
    else:
        _validate_coa_protonation_manifest(curation, require_source_provenance=True)
    if curation["activation_state"] != "approved":
        raise CoAProtonationActivationBlocked(
            "CoA protonation curation is blocked: "
            f"{curation.get('activation_reason', 'no approval recorded')}"
        )

    report = _evaluate_coa_protonation_gate(
        model, curation, source_model_path=source_model_path
    )
    if not report["gate_passed"]:
        raise CoAProtonationCurationError(
            "CoA protonation safety gate failed: "
            f"newly_unbalanced={report['newly_unbalanced']}; "
            f"non_h_charge_deltas={report['non_h_charge_deltas']}; "
            f"target_reactions_unbalanced={report['target_reactions_unbalanced']}; "
            f"source_snapshot_verified={report['source_snapshot']['verified']}; "
            f"balanced={report['balanced_count_before']}->{report['balanced_count_after']}; "
            f"objective={report['objective_before']}->{report['objective_after']}"
        )
    counts = _apply_coa_protonation_curation(model, curation)
    logger.info(
        "CoA protonation curation applied: metabolites=%d annotations=%d reactions=%d",
        counts["metabolites"],
        counts["annotations"],
        counts["reaction_corrections"],
    )
    return counts


# ── Top-level driver ──────────────────────────────────────────────────────

def apply_all_patches(model) -> dict:
    """
    Apply all known model patches. Returns a dict with per-patch counts.
    Safe to call multiple times (idempotent — already-correct values are skipped).
    """
    logger.info("Applying iYli21 known-bug patches …")
    # NOTE: fix_ec_code_format is intentionally NOT called here.  Reaction EC
    # codes are populated later in the pipeline (gene EC enrichment, EC
    # backfill, reaction xref backfill), so EC formatting must run after those
    # steps — it is invoked separately near the end of main().
    # NOTE: fix_coa_charge is intentionally NOT called.  Setting free CoA to
    # charge -4 in isolation unbalances every reaction that produces/consumes
    # it (the other side is not adjusted), which regressed Memote's
    # reaction_charge_balance (+10 offenders).  The function is kept for
    # reference but disabled until a charge fix is applied consistently across
    # whole reactions rather than to the metabolite alone.
    counts = {
        "nadp_plus_fixed":   fix_nadp_plus_formula(model),
        "ceramide_fixed":    fix_ceramide_formulas(model),
        "cation_formula_fixed": fix_cation_formula_consistency(model),
    }
    logger.info(
        f"  patches applied: NADP+={counts['nadp_plus_fixed']} copies, "
        f"ceramide={counts['ceramide_fixed']} copies, "
        f"cation-formula={counts['cation_formula_fixed']} copies"
    )
    return counts
