"""
patches.py — known data-bug patches for the iYali26 model.

Each patch fixes a discrete, documented error in the source SBML.  Patches
are applied early in the pipeline, immediately after metabolite annotation,
so that downstream mass-balance, FBA, and gap analyses operate on a
chemically correct model.

Currently applied patches:

  1. NADP+ formula fix
     iYali26 stores NADP+ as C21H28N7O17P3 (one H short of the KEGG C00006
     reference C21H29N7O17P3).  This affects 6 compartmental copies and
     ~126 reactions that use NADP+/NADPH, which had cascading effects on
     ceramide synthesis (see patch 2).

  2. Ceramide formula corrections
     iYali26 stores ceramide-1-(C24/C26) with formula C63H125NO6, which is
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

import copy
import hashlib
import json
from collections import Counter
from copy import deepcopy
from math import isclose, isfinite
from numbers import Real
from pathlib import Path
import logging
import math
import re

from .config import ESSENTIALITY_DIR, load_project_paths

logger = logging.getLogger(__name__)


# ── Patch 1: NADP+ formula fix ────────────────────────────────────────────

# iYali26 NADP+ copies (all compartments). Identified by metabolite name
# prefix "NADP(+)_" which is the iYali26 naming convention.
# Formula correction: C21H28N7O17P3 → C21H29N7O17P3 (KEGG C00006).
_NADP_PLUS_OLD_FORMULA = "C21H28N7O17P3"
_NADP_PLUS_NEW_FORMULA = "C21H29N7O17P3"


def _is_nadp_plus(met) -> bool:
    """Identify NADP+ by name (handles iYali26 'NADP(+)_C21H28N7O17P3' naming)."""
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
# Base name is what comes BEFORE the trailing "_" or "_FORMULA" in iYali26 names.
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
    """Return the iYali26 base name (lowercased, strip trailing _ and _FORMULA)."""
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

# iYali26 was exported from Excel; some metabolite names had their trailing
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


# ── Patch 0b: D-arabinokinase direction and proton bookkeeping ──────────

# R2041 was inherited as a reversible ATP:D-arabinose 5-phosphotransferase.
# The verified EC 2.7.1.54 equation is:
#
#   ATP + D-arabinose -> ADP + D-arabinose 5-phosphate + H+
#
# ENZYME: https://enzyme.expasy.org/EC/2.7.1.54
#
# Allowing the reverse direction lets the D-arabinose-phosphate isomerase /
# phosphatase loop synthesize ATP.  Directionality and the missing product
# proton are reaction-level facts and can be corrected independently of the
# deferred ATP/ADP connected-component microspecies migration.  With the
# current legacy neutral-H nucleotide formulas the corrected reaction retains
# an H=-1 formula residual; it becomes fully balanced only when ATP and ADP are
# migrated atomically to the curated pH-7.3 microspecies.
_D_ARABINOKINASE_REACTION_ID = "R2041"
_D_ARABINOKINASE_PROTON_ID = "m10[C_cy]"
_D_ARABINOKINASE_SOURCE = "https://enzyme.expasy.org/EC/2.7.1.54"


def fix_d_arabinokinase_direction_and_proton(model) -> int:
    """Make R2041 forward-only and add its verified product proton.

    The patch is idempotent and fails closed if R2041 already contains an
    unexpected proton coefficient.  It returns one when the reaction changed
    and zero when it was already canonical.
    """

    try:
        reaction = model.reactions.get_by_id(_D_ARABINOKINASE_REACTION_ID)
        proton = model.metabolites.get_by_id(_D_ARABINOKINASE_PROTON_ID)
    except KeyError:
        return 0

    proton_coefficient = float(reaction.metabolites.get(proton, 0.0))
    if proton_coefficient not in {0.0, 1.0}:
        raise ValueError(
            f"{_D_ARABINOKINASE_REACTION_ID} has unexpected product-proton "
            f"coefficient {proton_coefficient}"
        )

    changed = False
    if reaction.bounds != (0.0, 1000.0):
        reaction.bounds = (0.0, 1000.0)
        changed = True
    if proton_coefficient == 0.0:
        reaction.add_metabolites({proton: 1.0})
        changed = True

    notes = dict(reaction.notes or {})
    expected_note = {
        "source": _D_ARABINOKINASE_SOURCE,
        "equation": (
            "ATP + D-arabinose -> ADP + D-arabinose 5-phosphate + H+"
        ),
        "status": "active",
        "remaining_gate": "ATP/ADP connected-component microspecies migration",
        "lock_proton_water_stoichiometry": True,
    }
    if notes.get("curated_reaction_correction") != expected_note:
        notes["curated_reaction_correction"] = expected_note
        reaction.notes = notes
        changed = True

    return int(changed)


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
        audit_csv = str(
            load_project_paths().resolve_legacy_path(_EC_OVERLOAD_AUDIT_CSV)
        )

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


# ── Curated EC cleanup: ADP/ATP transporters ────────────────────────────────

# R815 is the mitochondrial adenine-nucleotide antiporter and R816 is a
# peroxisomal adenine-nucleotide transporter. Both inherited EC 2.7.4.6 from
# a broad MetaNetX/xref annotation. That EC describes a nucleoside-diphosphate
# kinase, not either transporter. Do not replace it with a guessed transport
# EC: current Y. lipolytica evidence supports carrier identity but not a
# reaction-specific EC assignment. This patch is metadata only: equation,
# bounds, GPR and non-EC cross-references are preserved.
#
# It runs after generic EC backfills, preventing a later annotation stage from
# silently reintroducing the demonstrated stale value.
_ADP_ATP_TRANSPORT_REACTION_IDS = ("R815", "R816")
_STALE_ADP_ATP_TRANSPORT_EC = "2.7.4.6"


def remove_stale_adp_atp_transporter_ec_codes(model) -> int:
    """Remove EC 2.7.4.6 from the two curated ADP/ATP transporters.

    Returns the number of reactions whose annotation changed. The operation is
    idempotent and removes only the demonstrated stale EC; an independently
    curated future transport EC would be retained.
    """

    changed = 0
    for reaction_id in _ADP_ATP_TRANSPORT_REACTION_IDS:
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            logger.warning(
                "  ADP/ATP transporter EC cleanup: %s not in model — skipping",
                reaction_id,
            )
            continue
        annotation = reaction.annotation if isinstance(reaction.annotation, dict) else {}
        raw_ec_codes = annotation.get("ec-code", [])
        ec_codes = [raw_ec_codes] if isinstance(raw_ec_codes, str) else list(raw_ec_codes)
        retained = [
            str(ec_code).strip()
            for ec_code in ec_codes
            if str(ec_code).strip()
            and str(ec_code).strip() != _STALE_ADP_ATP_TRANSPORT_EC
        ]
        if len(retained) == len(ec_codes):
            continue
        if retained:
            annotation["ec-code"] = sorted(set(retained))
        else:
            annotation.pop("ec-code", None)
        reaction.annotation = annotation
        changed += 1
        logger.info(
            "  ADP/ATP transporter EC cleanup: %s removed stale EC %s",
            reaction_id,
            _STALE_ADP_ATP_TRANSPORT_EC,
        )
    return changed


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
        additions_csv = str(
            load_project_paths().resolve_legacy_path(_GPR_ADDITIONS_CSV)
        )

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


# ── Patch 8b: remove mis-annotated genes from reaction GPRs ────────────────
#
# Two genes were assigned to reactions whose enzyme they do NOT catalyse — a
# reaction-mis-annotation (the gene's real protein is a different enzyme). Both
# are OR-isozyme partners, so removing them leaves the true catalyst in place and
# does NOT change any FBA/growth result — it only corrects GPR biology and removes
# a false-negative contaminant from essential-gene screening.
#
# Verified this session against UniProt (opened):
#   YALI1E07744g -> UniProt Q6C6P1 = glycoside hydrolase family 65 (trehalase),
#     NOT transketolase (EC 2.2.1.1). Wrongly placed in R765/R766; real
#     transketolase is the partner YALI1D02625g.
#     https://rest.uniprot.org/uniprotkb/Q6C6P1.json
#   YALI1E11370g -> similar to PET112/GatB, UniProt P33893 = glutamyl-tRNA(Gln)
#     amidotransferase subunit B (EC 6.3.5.-), NOT prephenate dehydrogenase.
#     Wrongly placed in R671; real prephenate dehydrogenase is the partner
#     YALI1F23441g.  https://rest.uniprot.org/uniprotkb/P33893.json
#
# (R671 also carries a self-contradictory EC 6.3.5.7 vs its "prephenate
# dehydrogenase" name — flagged for separate curation, NOT changed here.)
_GPR_MISANNOTATION_REMOVALS = {
    "R765": "YALI1E07744g",
    "R766": "YALI1E07744g",
    "R671": "YALI1E11370g",
}


def remove_misannotated_gprs(model) -> int:
    """Remove mis-annotated isozyme genes from reaction GPRs (safe: all 'or' cases,
    true partner retained). Returns the number of (reaction, gene) removals made.
    Idempotent; refuses to empty a GPR."""
    removed = 0
    for rid, gene in _GPR_MISANNOTATION_REMOVALS.items():
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            logger.warning(f"  GPR removal: reaction {rid} not in model — skipping")
            continue
        toks = set(re.findall(r"YALI1[A-Za-z0-9]+", rxn.gene_reaction_rule))
        if gene not in toks:
            continue  # idempotent: already removed
        remaining = sorted(t for t in toks if t != gene)
        if not remaining:
            logger.warning(f"  GPR removal: {rid} -= {gene} would empty GPR — skipping")
            continue
        rxn.gene_reaction_rule = " or ".join(remaining) if len(remaining) > 1 else remaining[0]
        removed += 1
        logger.info(f"  GPR removal (mis-annotation): {rid} -= {gene} -> '{rxn.gene_reaction_rule}'")
    return removed


# ── Patch 8c: reviewed R612 / external-NDH2 GPR corrections ─────────────
#
# These corrections have direct, reaction-specific evidence and are deliberately
# separate from the automated isozyme-expansion table above:
#
# * R612: YALI1E31685g / URA3 encodes orotidine-5'-phosphate decarboxylase.
#   The PO1f ura3-302 phenotype remains represented by the runtime strain
#   overlay, which disables R612 only in that strain context.
# * R570: YALI1F32476g / NDH2 is the external alternative NADH dehydrogenase.
#   R2063 is an exact, same-bounds duplicate of R570 and is removed.  The
#   mitochondrial complex-I reaction R1889 is intentionally *not* assigned an
#   AND GPR here; its subunit inventory is retained in
#   docs/curation/complex_i_gpr_evidence.csv pending identity/requiredness review.

_R612_URA3_GENE = "YALI1E31685g"
_R570_NDH2_GENE = "YALI1F32476g"


def _upsert_gene_annotation(gene, *, name: str, annotation: dict) -> None:
    """Set curated identity fields without discarding prior identifier mapping."""
    gene.name = name
    current = dict(gene.annotation) if isinstance(gene.annotation, dict) else {}
    current.update(annotation)
    gene.annotation = current


def add_r612_ura3_gpr(model) -> int:
    """Assign the experimentally supported URA3 GPR to R612.

    Returns one when the GPR changes and zero when it was already correct.
    Refuses an unexpected pre-existing rule so a future source-model change
    cannot be silently overwritten.
    """
    try:
        reaction = model.reactions.get_by_id("R612")
    except KeyError as exc:
        raise ValueError("R612 is required for the curated URA3 GPR correction") from exc

    target_rule = _R612_URA3_GENE
    old_rule = reaction.gene_reaction_rule.strip()
    if old_rule not in {"", target_rule}:
        raise ValueError(
            "R612 has an unexpected GPR; refusing to replace it: "
            f"{old_rule!r}"
        )

    changed = int(old_rule != target_rule)
    reaction.gene_reaction_rule = target_rule
    gene = model.genes.get_by_id(_R612_URA3_GENE)
    _upsert_gene_annotation(
        gene,
        name="URA3",
        annotation={
            "sbo": "SBO:0000243",
            "uniprot": "A0A1H6PUU4",
            "ec-code": "4.1.1.23",
        },
    )
    notes = dict(reaction.notes) if isinstance(reaction.notes, dict) else {}
    notes["curated_gpr_correction"] = (
        "YALI1E31685g (URA3), orotidine-5'-phosphate decarboxylase; "
        "direct Yarrowia lipolytica ura3 genetic evidence: "
        "https://onlinelibrary.wiley.com/doi/full/10.1002/yea.3673"
    )
    notes["strain_context"] = (
        "PO1f ura3-302 is represented by the runtime PO1f overlay, which "
        "sets R612 bounds to zero without changing the reference-model GPR."
    )
    notes["chemistry_status"] = (
        "UMP microspecies remains component_review pending a family-wide "
        "pyrimidine/nucleotide balance repair."
    )
    reaction.notes = notes
    if changed:
        logger.info("  Reviewed GPR correction: R612 = %s (URA3)", target_rule)
    return changed


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
    if old_rule not in {"", target_rule}:
        legacy_tokens = set(re.findall(r"[A-Za-z0-9_]+", old_rule))
        required_legacy_markers = {
            "YALI1A21711g",
            "YALI1F32476g",
            "YALIfMp29",
        }
        if not required_legacy_markers <= legacy_tokens:
            raise ValueError(
                "R570 has an unexpected GPR; refusing to replace it: "
                f"{old_rule!r}"
            )
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


# ── Patch 8d: direct enzyme-like GPR assignments ─────────────────────────
#
# These are reaction-specific assignments that are supported either by direct
# Yarrowia experiments (LIP2) or by a shared, already-modelled enzyme activity
# with the same EC and substrate family (the three vitamin-B6 kinases).  They
# are intentionally not placed in the broad automatic isozyme-expansion table.
#
# * R2274: YALI1A21372g / LIP2 is the secreted extracellular lipase.  lip2Δ
#   eliminates extracellular lipase activity, and Lip2 has direct activity on
#   long-chain triacylglycerols including triolein.
# * R1302/R1303: YALI1A08512g has curated pyridoxal-kinase EC 2.7.1.35
#   annotation and already catalyses the pyridoxine congener reaction R1306.
#   EC 2.7.1.35 covers pyridoxal, pyridoxamine and pyridoxine substrate forms.
#   This is a curated annotation assignment, not direct knockout evidence.

_DIRECT_ENZYME_LIKE_GPRS = {
    "R2274": "YALI1A21372g",
    "R1302": "YALI1A08512g",
    "R1303": "YALI1A08512g",
}


def add_direct_enzyme_like_gprs(model) -> int:
    """Assign three reviewed direct enzyme-like GPRs.

    The source rules must be empty (or already equal to the curated rule), so
    a future source-model change cannot be silently replaced.  Returns the
    number of reaction GPRs changed and is idempotent.
    """
    changed = 0
    for reaction_id, target_rule in _DIRECT_ENZYME_LIKE_GPRS.items():
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError as exc:
            raise ValueError(
                f"{reaction_id} is required for the direct enzyme-like GPR correction"
            ) from exc
        old_rule = reaction.gene_reaction_rule.strip()
        if old_rule not in {"", target_rule}:
            raise ValueError(
                f"{reaction_id} has an unexpected GPR; refusing to replace it: "
                f"{old_rule!r}"
            )
        if old_rule != target_rule:
            reaction.gene_reaction_rule = target_rule
            changed += 1
            logger.info(
                "  Direct enzyme-like GPR correction: %s = %s",
                reaction_id,
                target_rule,
            )

        notes = dict(reaction.notes) if isinstance(reaction.notes, dict) else {}
        if reaction_id == "R2274":
            notes["curated_gpr_correction"] = (
                "YALI1A21372g (LIP2), secreted extracellular triacylglycerol "
                "lipase; direct Yarrowia lipolytica lip2Δ and Lip2 triolein "
                "activity evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC101989/"
            )
            notes["gpr_evidence_status"] = "experimentally_verified"
        else:
            notes["curated_gpr_correction"] = (
                "YALI1A08512g, pyridoxal kinase EC 2.7.1.35; curated UniProt "
                "automatic annotation, with the same gene already assigned to "
                "the pyridoxine-congener reaction R1306. EC 2.7.1.35 covers "
                "pyridoxal, pyridoxamine and pyridoxine."
            )
            notes["gpr_evidence_status"] = "curated_annotation"
        reaction.notes = notes

    lip2 = model.genes.get_by_id("YALI1A21372g")
    _upsert_gene_annotation(
        lip2,
        name="LIP2",
        annotation={"sbo": "SBO:0000243", "ec-code": "3.1.1.3"},
    )
    pyridoxal_kinase = model.genes.get_by_id("YALI1A08512g")
    current = (
        dict(pyridoxal_kinase.annotation)
        if isinstance(pyridoxal_kinase.annotation, dict)
        else {}
    )
    current.update(
        {
            "sbo": "SBO:0000243",
            "uniprot": "A0A1D8N468",
            "ec-code": "2.7.1.35",
        }
    )
    pyridoxal_kinase.annotation = current
    return changed


# ── Patch 8e: remove the spurious cytosolic/nuclear quinone branches ────
#
# The source model contains two disconnected alternatives to the mitochondrial
# CoQ pathway:
#
# * R189 has 4-hydroxybenzoate/nonaprenyl-diphosphate chemistry, but is named
#   CAAX farnesyltransferase and assigned the RAM1/RAM2-like protein
#   farnesyltransferase subunits. Its branch-specific substrate and product are
#   degree-one endpoints, so steady state forces its flux to zero.
# * R2250 -> R2247 -> R2248 -> R2249 -> R2242 is a bacterial-style
#   octaprenyl/2-polyprenyl-6-hydroxyphenol fragment split between cytosol and
#   nucleus. R2250 joins an octaprenyl-labelled substrate to a chemically
#   hexaprenyl product, while R2242 is assigned DIM1 (EC 2.1.1.183), an 18S-rRNA
#   dimethyltransferase. Eukaryotic COQ3 instead performs mitochondrial
#   O-methylations on different, carboxylated CoQ intermediates.
#
# A separate mitochondrial COQ2 reaction, R407, already represents the retained
# pathway. Removing these six reactions therefore deletes an inert duplicate/
# mis-annotation; it does not open a pathway, add a demand, or tune essentiality
# recall. Native Y. lipolytica CoQ9 chain-length repair and decomposition of the
# broad COQ-synthome GPR are handled by separate, independently gated patches.
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

# Exact reviewed signatures prevent this cleanup from deleting a future,
# legitimately reconstructed CoQ9 reaction that happens to reuse a legacy ID.
# R189 has two explicitly accepted forms: the raw source contains two product
# protons, while the charge-aware pipeline removes them. The other five
# reactions are identical in the raw and pre-cleanup canonical models.
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


# ── Patch 8f: replace the legacy CoQ6 main chain with balanced CoQ9 ──

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

# R2063 is still present at the pre-FVA pipeline stage and is removed later as
# the reviewed duplicate of R570.  Supporting exactly that optional member lets
# FVA operate on Q9 without weakening the connected-component precondition.
_COQ9_OPTIONAL_CONNECTED_REACTION = "R2063"

# name, legacy Q6 formula, Q9 formula, charge, exact Q9 annotations
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


# ── Patch 8g: decompose the inherited quinone-synthome GPR ───────────────

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


# ── Patch 8c: remove spurious carrier-free CoA-thioester transport ─────────
#
# R1172 models 3-hydroxy-3-methylglutaryl-CoA (HMG-CoA) crossing the inner
# mitochondrial membrane directly and reversibly, with no carrier and no GPR:
#     m646[C_cy] <=> m648[C_mi]   (HMG-CoA cyt <-> HMG-CoA mito)
# HMG-CoA is a CoA-thioester, and the inner mitochondrial membrane is impermeable
# to acyl-CoA / CoA-thioesters — acyl groups cross only as carnitine esters via the
# carnitine shuttle, and HMG-CoA has no such carrier. Verified this session (opened):
#   "Fatty acyl CoA is impermeable to the inner mitochondrial membrane, so it is
#    carried in the form of fatty acyl carnitine."
#     https://library.med.utah.edu/NetBiochem/FattyAcids/8_4.html
#   "Since the mitochondrial inner membrane is not permeable to acyl-CoAs, acyl
#    groups are transferred from CoA to carnitine..."
#     https://pmc.ncbi.nlm.nih.gov/articles/PMC8066319/
# Both HMG-CoA pools are synthesized independently in their own compartment
# (cytosol R411; mito R412/R1973), so the transport is not needed for connectivity.
# Under SD-Leu- it carries zero flux and has no GPR, so removing it is WT-safe and
# changes no gene's essentiality — this is a model-correctness fix, not a recall fix.
# User-requested retention (2026-09-09); historical concern is preserved above
# and in data/reference_build/retained_reactions.json. No carrier GPR is invented.
_SPURIOUS_TRANSPORT_REMOVALS = []


def remove_spurious_transport_reactions(model) -> int:
    """Remove biochemically impossible transport reactions (carrier-free CoA-thioester
    crossing the inner mitochondrial membrane). Returns the number removed. Idempotent;
    keeps the shared metabolites (remove_orphans=False) since other reactions use them."""
    removed = 0
    for rid in _SPURIOUS_TRANSPORT_REMOVALS:
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue  # idempotent: already removed
        model.remove_reactions([rxn], remove_orphans=False)
        removed += 1
        logger.info(f"  Spurious transport removal: {rid} (carrier-free CoA-thioester transport)")
    return removed


# ── Evidence-gated essentiality curation patches ──────────────────────────

_ESSENTIALITY_PATCHES_CSV = ESSENTIALITY_DIR / "curated_model_patches.csv"
_ESSENTIALITY_OPERATIONS = {
    "set_gpr",
    "set_bounds",
    "remove_reaction",
    "couple_trna_biomass",
    "partition_cpa_ura2",
}
_ESSENTIALITY_VALUE_OPERATIONS = {
    "set_gpr",
    "set_bounds",
    "couple_trna_biomass",
}


def partition_cpa_ura2_pathways(model) -> dict[str, str]:
    """Separate the arginine CPA1/CPA2 pool from the channelled URA2 pathway.

    YALI1E11768g (URA2) is a multifunctional pyrimidine enzyme whose carbamoyl
    phosphate intermediate is channelled directly from its CPS domains to its
    ATCase domain.  The source model instead exposes that intermediate through
    the same reaction and metabolite pool used by the arginine-specific CPA1/
    CPA2 complex.  That representation makes all three genes interchangeable.

    This patch keeps the existing reaction and metabolite counts unchanged:

    * R159 becomes the net, channelled URA2 CPS + ATCase reaction;
    * R190 becomes the arginine-specific CPA1 AND CPA2 reaction;
    * m325 remains the explicit carbamoyl-phosphate pool used by R190/R607.

    The physical mitochondrial relocation is deliberately deferred because the
    current GEM lacks an evidence-backed glutamine/citrulline transport module.
    The function is idempotent and returns an auditable before/after snapshot.
    """
    required_reactions = ("R159", "R190", "R607")
    required_metabolites = (
        "m10[C_cy]",
        "m32[C_cy]",
        "m35[C_cy]",
        "m50[C_cy]",
        "m130[C_cy]",
        "m141[C_cy]",
        "m143[C_cy]",
        "m199[C_cy]",
        "m267[C_cy]",
        "m325[C_cy]",
        "m326[C_cy]",
    )
    try:
        reactions = {rid: model.reactions.get_by_id(rid) for rid in required_reactions}
        metabolites = {
            mid: model.metabolites.get_by_id(mid) for mid in required_metabolites
        }
    except KeyError as exc:
        raise ValueError(f"CPA/URA2 partition requires {exc.args[0]}") from exc

    r159 = reactions["R159"]
    r190 = reactions["R190"]
    r607 = reactions["R607"]
    carbamoyl_phosphate = metabolites["m325[C_cy]"]

    def assign_subsystem(reaction, subsystem_name: str) -> None:
        """Keep both the reaction attribute and SBML Groups package in sync."""
        subsystem_groups = [
            group
            for group in model.groups
            if group.kind == "partonomy"
            and group.annotation.get("sbo") == "SBO:0000633"
        ]
        target_groups = [
            group for group in subsystem_groups if group.name == subsystem_name
        ]
        if len(target_groups) != 1:
            raise ValueError(
                f"CPA/URA2 partition requires one subsystem group named {subsystem_name!r}"
            )
        target_group = target_groups[0]
        for group in subsystem_groups:
            if reaction in group.members and group is not target_group:
                group.remove_members([reaction])
        if reaction not in target_group.members:
            target_group.add_members([reaction])
        reaction.subsystem = subsystem_name

    if (
        carbamoyl_phosphate not in r190.metabolites
        or carbamoyl_phosphate not in r607.metabolites
    ):
        raise ValueError("CPA/URA2 partition requires m325 in both R190 and R607")

    before = (
        f"R159={r159.reaction} [{r159.gene_reaction_rule}]; "
        f"R190={r190.reaction} [{r190.gene_reaction_rule}]"
    )

    # Net stoichiometry of the source R190 CPS reaction plus R159 ATCase
    # reaction, with the channelled carbamoyl-phosphate intermediate cancelled.
    target_r159 = {
        "m10[C_cy]": -1.0,
        "m32[C_cy]": -1.0,
        "m35[C_cy]": 2.0,
        "m50[C_cy]": 1.0,
        "m130[C_cy]": -1.0,
        "m141[C_cy]": -2.0,
        "m143[C_cy]": 2.0,
        "m199[C_cy]": -1.0,
        "m267[C_cy]": -1.0,
        "m326[C_cy]": 1.0,
    }
    if r159.metabolites:
        r159.add_metabolites(
            {metabolite: -coefficient for metabolite, coefficient in r159.metabolites.items()}
        )
    r159.add_metabolites(
        {metabolites[mid]: coefficient for mid, coefficient in target_r159.items()}
    )
    r159.name = (
        "pyrimidine-specific carbamoyl-phosphate synthase/aspartate "
        "carbamoyltransferase (channelled)"
    )
    assign_subsystem(r159, "Pyrimidine metabolism")
    r159.gene_reaction_rule = "YALI1E11768g"
    r159.annotation["ec-code"] = ["2.1.3.2", "3.5.1.2", "6.3.4.16", "6.3.5.5"]
    r159.notes["essentiality_curation"] = (
        "URA2-only net reaction; its carbamoyl-phosphate intermediate is "
        "channelled and is not shared with arginine biosynthesis."
    )

    r190.name = "arginine-specific carbamoyl-phosphate synthase (glutamine-hydrolysing)"
    assign_subsystem(r190, "Arginine and proline metabolism")
    r190.gene_reaction_rule = "YALI1C33005g and YALI1D09420g"
    r190.annotation["ec-code"] = ["3.5.1.2", "6.3.4.16", "6.3.5.5"]
    r190.notes["essentiality_curation"] = (
        "Arginine-specific CPA2/CPA1 heterodimer; no longer interchangeable with URA2."
    )
    carbamoyl_phosphate.notes["pathway_pool"] = (
        "Explicit arginine-specific pool produced by R190 and consumed by R607; "
        "the URA2 intermediate is channelled inside R159."
    )

    after = (
        f"R159={r159.reaction} [{r159.gene_reaction_rule}]; "
        f"R190={r190.reaction} [{r190.gene_reaction_rule}]"
    )
    return {"before": before, "after": after}


def couple_trna_charging_to_biomass(
    model,
    biomass_id: str = "biomass_C",
    template_id: str = "R1387",
) -> list[dict[str, str | float]]:
    """Replace free amino-acid biomass drains with charged-tRNA carriers.

    The source model already contains 20 cytosolic aminoacyl-tRNA synthetase
    reactions and a locked, carrier-balanced tRNA biomass template. This patch
    uses the active biomass amino-acid coefficients but requires each amino acid
    to pass through its charging reaction: biomass consumes charged tRNA and
    returns the corresponding uncharged tRNA with the same coefficient.

    This preserves each tRNA carrier exactly and is idempotent. It refuses to
    make a partial change unless all 20 amino-acid/tRNA mappings are complete.
    """
    candidate_ids = {"R1387", "R1710"}
    try:
        biomass = model.reactions.get_by_id(biomass_id)
        template = model.reactions.get_by_id(template_id)
    except KeyError as exc:
        raise ValueError(
            f"tRNA coupling requires reactions {biomass_id} and {template_id}"
        ) from exc

    mappings: list[tuple[object, object, object, object, float]] = []
    already_coupled = 0
    charged_metabolites = sorted(
        [
            metabolite
            for metabolite, coefficient in template.metabolites.items()
            if coefficient < 0 and "trna" in (metabolite.name or "").lower()
        ],
        key=lambda metabolite: metabolite.id,
    )
    for charged in charged_metabolites:
        charging_reactions = sorted(
            [reaction for reaction in charged.reactions if reaction.id not in candidate_ids],
            key=lambda reaction: reaction.id,
        )
        if len(charging_reactions) != 1:
            raise ValueError(
                f"Expected one charging reaction for {charged.id}; found "
                f"{[reaction.id for reaction in charging_reactions]}"
            )
        charging = charging_reactions[0]
        uncharged = [
            metabolite
            for metabolite, coefficient in charging.metabolites.items()
            if coefficient < 0 and "trna" in (metabolite.name or "").lower()
        ]
        amino_acids = [
            metabolite
            for metabolite, coefficient in charging.metabolites.items()
            if coefficient < 0
            and metabolite in biomass.metabolites
            and biomass.metabolites[metabolite] < 0
            and (metabolite.name or "").lower().startswith("l-")
        ]

        if len(uncharged) != 1:
            raise ValueError(
                f"Expected one uncharged tRNA for {charged.id}; found "
                f"{[metabolite.id for metabolite in uncharged]}"
            )
        uncharged_metabolite = uncharged[0]
        if uncharged_metabolite not in template.metabolites or not math.isclose(
            abs(float(template.metabolites[charged])),
            abs(float(template.metabolites[uncharged_metabolite])),
            rel_tol=1e-9,
            abs_tol=1e-12,
        ):
            raise ValueError(f"Template does not conserve tRNA carrier for {charged.id}")

        if not amino_acids:
            charged_coefficient = biomass.metabolites.get(charged, 0.0)
            uncharged_coefficient = biomass.metabolites.get(uncharged_metabolite, 0.0)
            if charged_coefficient < 0 and math.isclose(
                abs(float(charged_coefficient)),
                abs(float(uncharged_coefficient)),
                rel_tol=1e-9,
                abs_tol=1e-12,
            ):
                already_coupled += 1
                continue
            raise ValueError(f"No biomass amino-acid substrate found for {charged.id}")
        if len(amino_acids) != 1:
            raise ValueError(
                f"Ambiguous biomass amino acid for {charged.id}: "
                f"{[metabolite.id for metabolite in amino_acids]}"
            )
        amino_acid = amino_acids[0]
        amount = -float(biomass.metabolites[amino_acid])
        mappings.append((amino_acid, charged, uncharged_metabolite, charging, amount))

    if already_coupled == 20 and not mappings:
        return []
    if already_coupled or len(mappings) != 20:
        raise ValueError(
            f"Refusing partial tRNA coupling: {len(mappings)} new, {already_coupled} existing"
        )

    audit: list[dict[str, str | float]] = []
    for amino_acid, charged, uncharged, charging, amount in mappings:
        biomass.add_metabolites(
            {
                amino_acid: amount,
                charged: -amount,
                uncharged: amount,
            }
        )
        audit.append(
            {
                "amino_acid_id": amino_acid.id,
                "charged_trna_id": charged.id,
                "uncharged_trna_id": uncharged.id,
                "charging_reaction_id": charging.id,
                "coefficient": amount,
            }
        )
    biomass.notes["essentiality_trna_coupling"] = (
        "Free amino-acid biomass drains replaced by carrier-balanced charged-tRNA "
        "consumption and uncharged-tRNA recycling."
    )
    return audit


_SPLIT_TRNA_BIOMASS_MODE = "split_v1"
_SPLIT_TRNA_REACTION_PREFIX = "TRNA_BIOMASS_"
_SPLIT_TRNA_RESIDUE_PREFIX = "trna_biomass_residue_"


def split_trna_charging_from_biomass(
    model,
    biomass_id: str = "biomass_C",
    template_id: str = "R1387",
) -> list[dict[str, str | float]]:
    """Build the fully split, experimental B-group translation layer.

    For each of the 20 cytosolic amino acids, this overlay replaces the free
    amino-acid coefficient in ``biomass_id`` with a private protein-residue
    intermediate and adds one carrier-conserving reaction::

        AA-tRNA(i) -> tRNA(i) + protein_residue(i)

    The biomass reaction consumes ``a_i protein_residue(i)``. Consequently the
    corresponding aminoacyl-tRNA synthetase must carry ``a_i`` flux per unit of
    biomass, while every tRNA carrier is returned one-for-one. The 20 separate
    reactions make each amino-acid requirement independently auditable.

    This is an experimental overlay, not an evidence-approved canonical patch.
    It performs a complete preflight and refuses partial or mixed states.
    """
    from cobra import Metabolite, Reaction

    try:
        biomass = model.reactions.get_by_id(biomass_id)
    except KeyError as exc:
        raise ValueError(f"Split tRNA coupling requires reaction {biomass_id}") from exc

    split_reactions = sorted(
        (
            reaction
            for reaction in model.reactions
            if reaction.id.startswith(_SPLIT_TRNA_REACTION_PREFIX)
        ),
        key=lambda reaction: reaction.id,
    )
    split_residues = sorted(
        (
            metabolite
            for metabolite in model.metabolites
            if metabolite.id.startswith(_SPLIT_TRNA_RESIDUE_PREFIX)
        ),
        key=lambda metabolite: metabolite.id,
    )
    mode = str(biomass.notes.get("experimental_trna_biomass_mode", ""))
    if mode == _SPLIT_TRNA_BIOMASS_MODE:
        if len(split_reactions) != 20 or len(split_residues) != 20:
            raise ValueError(
                "Partial B-group tRNA biomass state: expected 20 split reactions "
                f"and residues, found {len(split_reactions)} and {len(split_residues)}"
            )
        for reaction in split_reactions:
            residue_id = str(reaction.notes.get("protein_residue_id", ""))
            charged_id = str(reaction.notes.get("charged_trna_id", ""))
            uncharged_id = str(reaction.notes.get("uncharged_trna_id", ""))
            amount = float(reaction.notes.get("biomass_coefficient", 0.0))
            required_ids = (residue_id, charged_id, uncharged_id)
            if not all(required_ids) or amount <= 0:
                raise ValueError(
                    f"Partial B-group metadata on split reaction {reaction.id}"
                )
            residue = model.metabolites.get_by_id(residue_id)
            charged = model.metabolites.get_by_id(charged_id)
            uncharged = model.metabolites.get_by_id(uncharged_id)
            expected = {
                charged: -1.0,
                uncharged: 1.0,
                residue: 1.0,
            }
            if dict(reaction.metabolites) != expected or not math.isclose(
                float(biomass.metabolites.get(residue, 0.0)),
                -amount,
                rel_tol=1e-9,
                abs_tol=1e-12,
            ):
                raise ValueError(
                    f"Partial B-group stoichiometry on split reaction {reaction.id}"
                )
        return []
    if mode or split_reactions or split_residues:
        raise ValueError(
            "Refusing partial or mixed B-group tRNA biomass state: "
            f"mode={mode!r}, reactions={len(split_reactions)}, "
            f"residues={len(split_residues)}"
        )

    # Reuse the direct-coupling preflight on a copy. It is the established
    # complete 20-pair mapper, but no mutation is allowed to escape the probe.
    probe = model.copy()
    mapping = couple_trna_charging_to_biomass(
        probe,
        biomass_id=biomass_id,
        template_id=template_id,
    )
    if len(mapping) != 20:
        raise ValueError(
            f"B-group tRNA biomass requires 20 complete pairs; found {len(mapping)}"
        )

    planned_ids: set[str] = set()
    for row in mapping:
        charging_id = str(row["charging_reaction_id"])
        planned_ids.update(
            {
                f"{_SPLIT_TRNA_REACTION_PREFIX}{charging_id}",
                f"{_SPLIT_TRNA_RESIDUE_PREFIX}{charging_id}",
            }
        )
    existing_ids = {
        item_id
        for item_id in planned_ids
        if item_id in model.reactions or item_id in model.metabolites
    }
    if existing_ids:
        raise ValueError(
            "Refusing B-group tRNA biomass ID collision: "
            + ", ".join(sorted(existing_ids))
        )

    audit: list[dict[str, str | float]] = []
    for row in mapping:
        amino_acid = model.metabolites.get_by_id(str(row["amino_acid_id"]))
        charged = model.metabolites.get_by_id(str(row["charged_trna_id"]))
        uncharged = model.metabolites.get_by_id(str(row["uncharged_trna_id"]))
        charging_id = str(row["charging_reaction_id"])
        amount = float(row["coefficient"])
        residue_id = f"{_SPLIT_TRNA_RESIDUE_PREFIX}{charging_id}"
        reaction_id = f"{_SPLIT_TRNA_REACTION_PREFIX}{charging_id}"

        residue = Metabolite(
            residue_id,
            name=(
                f"Protein-incorporated {amino_acid.name or amino_acid.id} "
                "requirement (experimental B split)"
            ),
            compartment=amino_acid.compartment,
        )
        residue.annotation = {"sbo": "SBO:0000247"}
        residue.notes = {
            "experimental_role": "split_trna_biomass_protein_residue",
            "source_amino_acid_id": amino_acid.id,
        }

        translation = Reaction(reaction_id)
        translation.name = (
            f"Biomass incorporation of {amino_acid.name or amino_acid.id} "
            "through charged tRNA"
        )
        translation.bounds = (0.0, 1000.0)
        translation.subsystem = "Protein synthesis"
        translation.annotation = {"sbo": "SBO:0000176"}
        translation.notes = {
            "experimental_trna_biomass_mode": _SPLIT_TRNA_BIOMASS_MODE,
            "amino_acid_id": amino_acid.id,
            "charged_trna_id": charged.id,
            "uncharged_trna_id": uncharged.id,
            "charging_reaction_id": charging_id,
            "protein_residue_id": residue_id,
            "biomass_coefficient": amount,
            "carrier_conservation": "1 AA-tRNA consumed; 1 tRNA returned",
        }
        translation.add_metabolites(
            {
                charged: -1.0,
                uncharged: 1.0,
                residue: 1.0,
            }
        )
        model.add_reactions([translation])
        biomass.add_metabolites({amino_acid: amount, residue: -amount})
        audit.append(
            {
                **row,
                "split_reaction_id": reaction_id,
                "protein_residue_id": residue_id,
            }
        )

    biomass.notes["experimental_trna_biomass_mode"] = _SPLIT_TRNA_BIOMASS_MODE
    biomass.notes["experimental_trna_biomass_template"] = template_id
    biomass.notes["experimental_trna_biomass_design"] = (
        "B group: 20 independent AA-tRNA -> tRNA + protein-residue reactions"
    )
    return audit


def _validate_schema_v2_patch_gate(
    model,
    row,
    repo_root: str,
    essentiality_dir: str,
    patch_id: str,
) -> dict:
    """Validate evidence, human approval and the live target fingerprint."""
    import json
    import os

    from .essentiality_evidence import (
        SIMULATION_CONTEXT_FIELDS,
        chemistry_fingerprint,
        read_ledger,
        require_valid_evidence_dossier,
        sha256_file,
        target_fingerprint,
    )

    required = (
        "case_id",
        "evidence_path",
        "approved_by",
        "approved_at",
        "target_fingerprint",
    )
    missing = [field for field in required if not (row.get(field) or "").strip()]
    if missing:
        raise ValueError(
            f"Essentiality patch {patch_id} lacks schema-v2 gate fields: {missing}"
        )
    if row["approved_by"].strip() != "human_user":
        raise ValueError(
            f"Essentiality patch {patch_id} approved_by must be human_user"
        )

    def resolve_essentiality_path(raw_path: str) -> str:
        if os.path.isabs(raw_path):
            return os.path.realpath(raw_path)
        normalized = raw_path.replace("\\", "/")
        legacy_prefix = "data/essentiality/"
        if normalized.startswith(legacy_prefix):
            return os.path.realpath(
                os.path.join(
                    essentiality_dir,
                    normalized[len(legacy_prefix) :],
                )
            )
        return os.path.realpath(os.path.join(repo_root, raw_path))

    evidence_path = resolve_essentiality_path(row["evidence_path"].strip())
    evidence_path = os.path.realpath(evidence_path)
    evidence_root = os.path.realpath(os.path.join(essentiality_dir, "evidence"))
    if os.path.commonpath([evidence_path, evidence_root]) != evidence_root:
        raise ValueError(
            f"Essentiality patch {patch_id} evidence must be inside data/essentiality/evidence"
        )
    if not os.path.exists(evidence_path):
        raise ValueError(
            f"Essentiality patch {patch_id} evidence dossier is missing: {evidence_path}"
        )
    with open(evidence_path, encoding="utf-8") as handle:
        dossier = json.load(handle)
    require_valid_evidence_dossier(dossier, require_human_approval=True)

    case_id = row["case_id"].strip()
    expected_fingerprint = row["target_fingerprint"].strip()
    human = dossier["human_decision"]
    proposal = dossier["proposed_operation"]
    if dossier.get("case_id") != case_id:
        raise ValueError(f"Essentiality patch {patch_id} case_id does not match dossier")
    if dossier.get("target_fingerprint") != expected_fingerprint:
        raise ValueError(
            f"Essentiality patch {patch_id} target fingerprint does not match dossier"
        )
    if human.get("approved_by") != row["approved_by"].strip():
        raise ValueError(f"Essentiality patch {patch_id} approver does not match dossier")
    if human.get("approved_at") != row["approved_at"].strip():
        raise ValueError(
            f"Essentiality patch {patch_id} approval timestamp does not match dossier"
        )
    if proposal.get("operation") != row["operation"].strip():
        raise ValueError(
            f"Essentiality patch {patch_id} operation does not match evidence proposal"
        )
    if proposal.get("target_id") != row["target_id"].strip():
        raise ValueError(
            f"Essentiality patch {patch_id} target_id does not match evidence proposal"
        )
    operation = row["operation"].strip()
    if operation in _ESSENTIALITY_VALUE_OPERATIONS:
        value = row["value"].strip()
        if not value:
            raise ValueError(
                f"Essentiality patch {patch_id} operation {operation} requires a value"
            )
        if proposal.get("value") != value:
            raise ValueError(
                f"Essentiality patch {patch_id} value does not match evidence proposal"
            )

    chemistry = dossier["chemistry_review"]
    chemistry_audit_path = resolve_essentiality_path(
        chemistry["audit_path"].strip()
    )
    if os.path.commonpath([chemistry_audit_path, evidence_root]) != evidence_root:
        raise ValueError(
            f"Essentiality patch {patch_id} chemistry audit must be inside "
            "data/essentiality/evidence"
        )
    if not os.path.exists(chemistry_audit_path):
        raise ValueError(
            f"Essentiality patch {patch_id} chemistry audit is missing: "
            f"{chemistry_audit_path}"
        )
    chemistry_audit_sha256 = sha256_file(chemistry_audit_path)
    if chemistry_audit_sha256 != chemistry["audit_sha256"].strip():
        raise ValueError(
            f"Essentiality patch {patch_id} chemistry audit SHA does not match dossier"
        )

    ledger_path = os.path.join(essentiality_dir, "curation_cases.csv")
    ledger_matches = [
        item for item in read_ledger(ledger_path) if item["case_id"] == case_id
    ]
    if len(ledger_matches) != 1:
        raise ValueError(
            f"Essentiality patch {patch_id} requires one durable ledger row; "
            f"found {len(ledger_matches)}"
        )
    ledger = ledger_matches[0]
    if ledger["status"] not in {"accepted", "implemented", "regression_passed"}:
        raise ValueError(
            f"Essentiality patch {patch_id} ledger status is {ledger['status']!r}, "
            "not accepted"
        )
    for field in ("target_fingerprint", "approved_by", "approved_at"):
        if ledger.get(field, "") != row[field].strip():
            raise ValueError(
                f"Essentiality patch {patch_id} {field} does not match durable ledger"
            )
    expected_chemistry_fingerprint = dossier["chemistry_fingerprint"]
    if ledger.get("chemistry_fingerprint", "") != expected_chemistry_fingerprint:
        raise ValueError(
            f"Essentiality patch {patch_id} chemistry fingerprint does not match "
            "durable ledger"
        )
    if str(dossier.get("schema_version", "2.0")) >= "2.1":
        for field in SIMULATION_CONTEXT_FIELDS:
            dossier_value = dossier.get(field)
            ledger_value = ledger.get(field, "")
            if field == "strain_overlay_enabled":
                matches = str(ledger_value).strip().lower() == str(bool(dossier_value)).lower()
            elif dossier_value is None:
                matches = str(ledger_value).strip() == ""
            else:
                matches = str(ledger_value) == str(dossier_value)
            if not matches:
                raise ValueError(
                    f"Essentiality patch {patch_id} {field} does not match durable ledger"
                )

    reaction_ids = sorted(
        {
            str(context.get("reaction_id", "")).strip()
            for context in dossier.get("model_context", {}).get("reactions", [])
            if str(context.get("reaction_id", "")).strip()
        }
    )
    if not reaction_ids:
        reaction_ids = [row["target_id"].strip()]
    live_contexts = []
    for reaction_id in reaction_ids:
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError as exc:
            raise ValueError(
                f"Essentiality patch {patch_id} evidence targets missing reaction {reaction_id}"
            ) from exc
        live_contexts.append(
            {
                "reaction_id": reaction.id,
                "stoichiometry": {
                    metabolite.id: float(coefficient)
                    for metabolite, coefficient in sorted(
                        reaction.metabolites.items(), key=lambda item: item[0].id
                    )
                },
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
                "gpr": reaction.gene_reaction_rule,
                "metabolite_chemistry": {
                    metabolite.id: {
                        "formula": metabolite.formula,
                        "charge": metabolite.charge,
                        "compartment": metabolite.compartment,
                    }
                    for metabolite in sorted(
                        reaction.metabolites, key=lambda item: item.id
                    )
                },
            }
        )
    live_fingerprint = target_fingerprint(live_contexts)
    if live_fingerprint != expected_fingerprint:
        raise ValueError(
            f"Essentiality patch {patch_id} is stale: target fingerprint changed "
            f"({expected_fingerprint} -> {live_fingerprint})"
        )
    live_chemistry_fingerprint = chemistry_fingerprint(live_contexts)
    if live_chemistry_fingerprint != expected_chemistry_fingerprint:
        raise ValueError(
            f"Essentiality patch {patch_id} is stale: chemistry fingerprint changed "
            f"({expected_chemistry_fingerprint} -> {live_chemistry_fingerprint})"
        )
    return {
        "case_id": case_id,
        "evidence_path": evidence_path,
        "approved_by": row["approved_by"].strip(),
        "approved_at": row["approved_at"].strip(),
        "target_fingerprint": expected_fingerprint,
        "chemistry_fingerprint": expected_chemistry_fingerprint,
        "chemistry_audit_path": chemistry_audit_path,
        "chemistry_audit_sha256": chemistry_audit_sha256,
        "audited_reaction_ids": ";".join(chemistry["audited_reaction_ids"]),
    }


def _assert_schema_v2_post_patch_balance(
    model, audited_reaction_ids: str, patch_id: str
) -> None:
    """Stop a schema-v2 build if an audited target is no longer balanced."""
    failures: dict[str, object] = {}
    for reaction_id in (
        item.strip() for item in audited_reaction_ids.split(";") if item.strip()
    ):
        try:
            reaction = model.reactions.get_by_id(reaction_id)
        except KeyError:
            failures[reaction_id] = "missing_after_patch"
            continue
        try:
            residual = reaction.check_mass_balance()
        except (TypeError, ValueError) as exc:
            failures[reaction_id] = f"uncheckable: {exc}"
            continue
        if residual:
            failures[reaction_id] = residual
    if failures:
        raise ValueError(
            f"Essentiality patch {patch_id} failed the post-patch mass/charge "
            f"gate: {failures}"
        )


def apply_curated_essentiality_patches(
    model,
    patches_csv: str | None = None,
) -> list[dict[str, str]]:
    """Apply reviewed essentiality patches from the audit table.

    Only ``status=accepted`` rows are eligible. Schema-v2 rows must also match
    a direct-evidence dossier, independent skeptic pass, durable human approval
    and the current target fingerprint. ``EG-GPR-001`` is the sole legacy-v1
    exception. Eligible rows may use one of five deliberately narrow operations:

    - ``set_gpr``: replace a reaction's gene-reaction rule with ``value``;
    - ``set_bounds``: set ``lower;upper`` bounds from ``value``;
    - ``remove_reaction``: remove the target reaction without removing orphans;
    - ``couple_trna_biomass``: use ``target_id`` as biomass and ``value`` as
      the locked tRNA template reaction.
    - ``partition_cpa_ura2``: replace the shared CPA1/CPA2/URA2 pool with a
      channelled URA2 reaction and an arginine-specific CPA1 AND CPA2 reaction.

    The function returns an audit list describing applied changes. Empty and
    review-only tables are valid and leave the model unchanged.
    """
    import csv
    import os

    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    legacy_essentiality_dir = os.path.join(root, "data", "essentiality")
    if patches_csv is None:
        patches_csv = os.fspath(_ESSENTIALITY_PATCHES_CSV)
        essentiality_dir = os.fspath(ESSENTIALITY_DIR)
    else:
        essentiality_dir = (
            legacy_essentiality_dir
            if os.path.isdir(legacy_essentiality_dir)
            else os.path.dirname(os.path.abspath(patches_csv))
        )
    if not os.path.exists(patches_csv):
        logger.warning("  Essentiality patch table not found: %s", patches_csv)
        return []

    required = {
        "patch_id",
        "status",
        "operation",
        "target_id",
        "value",
        "evidence_url",
        "rationale",
    }
    applied: list[dict[str, str]] = []
    accepted_rows: list[tuple[int, dict[str, str]]] = []
    with open(patches_csv, newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Essentiality patch table missing columns: {sorted(missing)}"
            )
        for line_number, row in enumerate(reader, start=2):
            if row["status"].strip().lower() != "accepted":
                continue
            accepted_rows.append((line_number, row))

    # Validate every accepted row against its evidence dossier before applying
    # any of them. This prevents a later stale or mismatched schema-v2 row from
    # leaving the caller's model partially patched.
    preflighted: list[dict[str, object]] = []
    for line_number, row in accepted_rows:
        patch_id = row["patch_id"].strip()
        operation = row["operation"].strip()
        target_id = row["target_id"].strip()
        value = row["value"].strip()
        evidence_url = row["evidence_url"].strip()
        rationale = row["rationale"].strip()
        schema_version_text = (row.get("schema_version") or "1").strip() or "1"
        try:
            schema_version = int(float(schema_version_text))
        except ValueError as exc:
            raise ValueError(
                f"Essentiality patch {patch_id} has invalid schema_version "
                f"{schema_version_text!r}"
            ) from exc
        if not patch_id or not target_id:
            raise ValueError(f"Essentiality patch row {line_number} lacks ID/target")
        if operation not in _ESSENTIALITY_OPERATIONS:
            raise ValueError(
                f"Essentiality patch {patch_id} has unsupported operation {operation!r}"
            )
        if not evidence_url.startswith(("https://", "http://")) or not rationale:
            raise ValueError(
                f"Essentiality patch {patch_id} requires evidence_url and rationale"
            )
        gate_audit: dict[str, str] = {}
        if schema_version >= 2:
            gate_audit = _validate_schema_v2_patch_gate(
                model,
                row,
                root,
                essentiality_dir,
                patch_id,
            )
        else:
            if patch_id != "EG-GPR-001":
                raise ValueError(
                    f"Essentiality patch {patch_id} cannot use the legacy "
                    "schema-v1 gate"
                )
            logger.warning(
                "  Essentiality patch %s uses legacy schema-v1 approval; "
                "evidence backfill is still required",
                patch_id,
            )
        try:
            reaction = model.reactions.get_by_id(target_id)
        except KeyError as exc:
            raise ValueError(
                f"Essentiality patch {patch_id} targets missing reaction {target_id}"
            ) from exc

        bounds: tuple[float, float] | None = None
        if operation == "set_gpr" and not value:
            raise ValueError(f"Essentiality patch {patch_id} has an empty GPR")
        if operation == "set_bounds":
            try:
                lower_text, upper_text = value.split(";", 1)
                bounds = (float(lower_text), float(upper_text))
            except (ValueError, TypeError) as exc:
                raise ValueError(
                    f"Essentiality patch {patch_id} bounds must be 'lower;upper'"
                ) from exc
            if bounds[0] > bounds[1]:
                raise ValueError(f"Essentiality patch {patch_id} has inverted bounds")

        preflighted.append(
            {
                "patch_id": patch_id,
                "operation": operation,
                "target_id": target_id,
                "value": value,
                "evidence_url": evidence_url,
                "rationale": rationale,
                "schema_version": schema_version,
                "gate_audit": gate_audit,
                "reaction": reaction,
                "bounds": bounds,
            }
        )

    for patch in preflighted:
        patch_id = str(patch["patch_id"])
        operation = str(patch["operation"])
        target_id = str(patch["target_id"])
        value = str(patch["value"])
        evidence_url = str(patch["evidence_url"])
        rationale = str(patch["rationale"])
        schema_version = int(patch["schema_version"])
        gate_audit = patch["gate_audit"]
        reaction = patch["reaction"]
        before = ""
        after = ""
        if operation == "set_gpr":
            before = reaction.gene_reaction_rule
            reaction.gene_reaction_rule = value
            after = reaction.gene_reaction_rule
        elif operation == "set_bounds":
            bounds = patch["bounds"]
            before = f"{reaction.lower_bound};{reaction.upper_bound}"
            reaction.bounds = bounds
            after = f"{reaction.lower_bound};{reaction.upper_bound}"
        elif operation == "remove_reaction":
            before = reaction.reaction
            model.remove_reactions([reaction], remove_orphans=False)
            after = "removed"
        elif operation == "couple_trna_biomass":
            template_id = value or "R1387"
            before = reaction.reaction
            mapping = couple_trna_charging_to_biomass(
                model,
                biomass_id=target_id,
                template_id=template_id,
            )
            after = reaction.reaction
            if not mapping and before == after:
                logger.info("  Essentiality patch %s already applied", patch_id)
        else:
            partition_audit = partition_cpa_ura2_pathways(model)
            before = partition_audit["before"]
            after = partition_audit["after"]

        if schema_version >= 2 and operation != "remove_reaction":
            _assert_schema_v2_post_patch_balance(
                model,
                str(gate_audit["audited_reaction_ids"]),
                patch_id,
            )

        audit_row = {
            "patch_id": patch_id,
            "schema_version": str(schema_version),
            "operation": operation,
            "target_id": target_id,
            "before": before,
            "after": after,
            "evidence_url": evidence_url,
            "rationale": rationale,
        }
        audit_row.update(gate_audit)
        applied.append(audit_row)
    return applied


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

    project_paths = load_project_paths()
    if additions_csv is None:
        additions_csv = str(
            project_paths.resolve_legacy_path(_GPR_ADDITIONS_CSV)
        )
    if not os.path.exists(additions_csv):
        logger.warning(f"  gene annotate: {additions_csv} not found — skipping")
        return 0

    # YALI1 (no underscore) -> GeneID from NCBI feature table
    y2g = {}
    ft = str(project_paths.resolve_legacy_path(_FEATURE_TABLE))
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
    kg = str(project_paths.resolve_legacy_path(_KEGG_GENES))
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
        fill_csv = str(
            load_project_paths().resolve_legacy_path(_NEUTRAL_FILL_CSV)
        )
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

# Y. lipolytica W29 makes ~8% palmitoleate (C16:1) but the iYali26 acyl-CoA pools
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


# ── Charge fix, Stage 1: free monocarboxylate/sulfonate anions left at charge 0 ──
#
# These 4 metabolites store the DEPROTONATED (anion) formula but were left at
# charge 0 — a missing-charge bug, not a wrong formula. Setting charge = -1 makes
# each metabolite formula/charge self-consistent and balances the reactions they
# touch, WITHOUT any formula edit. This is the "provably-safe" subset from the
# charge-convention strategy review: each touches only reactions that become
# balanced (or were already mass-imbalanced for other reasons), breaking none.
#
# Verified anchors (pH ~7.3 major microspecies):
#   m965 taurocholic acid C26H44NO7S -> -1: ChEBI:36257 taurocholate(1-), formula
#        matches the model exactly (verified this session)
#        https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:36257
#   m1403 kynurenic acid C10H7NO3 -> -1: monocarboxylic acid, KEGG C01717 (verified)
#        https://www.kegg.jp/entry/C01717
#   m286 ureidoglycolic acid C3H6N2O4 -> -1: carboxylate, KEGG C00603 (verified)
#        https://www.kegg.jp/entry/C00603
#   m753 1-pyrroline-3-hydroxy-5-carboxylate C5H6NO3 -> -1: name is the carboxylate;
#        formula already deprotonated; in-model R1194 balances and R500 improves
#        (charge gap -3 -> -2), consistent with a -1 anion (inferred, model-self-consistent)
#
# Note: m1403 = -1 also adds a charge term to R1375 (spontaneous kynurenic ->
# quinaldic), but R1375 is ALREADY mass-imbalanced (drops C10H7NO3 from nothing) —
# a pre-existing modeling error, flagged separately for curation. The -1 charge is
# the chemically correct value; we do not freeze it to flatter a broken reaction.
_CHARGE_STAGE1 = {
    "m286[C_cy]":  -1,   # ureidoglycolic acid
    "m753[C_mi]":  -1,   # 1-pyrroline-3-hydroxy-5-carboxylate
    "m965[C_cy]":  -1,   # taurocholic acid
    "m1403[C_cy]": -1,   # kynurenic acid
}


def fix_charge_stage1(model) -> int:
    """Set 4 free anionic metabolites from charge 0 to -1 (formula unchanged).

    Idempotent: a metabolite already at the target charge is skipped. Returns the
    number of metabolites changed. See _CHARGE_STAGE1 for per-metabolite sources.
    """
    changed = 0
    for mid, target in _CHARGE_STAGE1.items():
        try:
            met = model.metabolites.get_by_id(mid)
        except KeyError:
            continue
        if met.charge == target:
            continue
        old = met.charge
        met.charge = target
        changed += 1
        logger.info("  Charge Stage 1: %s (%s) charge %s -> %d"
                    % (mid, met.name, old, target))
    return changed


# ── Charge fix, Stage 2: anion-stored metabolites identified by InChI dH ──────
#
# Found by the InChI-dH discriminator (full-model scan): a metabolite whose stored
# formula has FEWER H than its own embedded-InChI neutral formula (dH < 0) is already
# the deprotonated anion, so setting charge = dH leaves the (formula, charge) pair a
# real protonation microspecies WITHOUT any formula edit. This avoids the FMN/FAD
# trap (those are dH == 0 = neutral-stored, and are deliberately NOT touched here).
#
# All 7 below were reaction-safety-gated on a cold model reload: balanced 914 -> 921
# (+7: R313/R322/R344/R567/R764/R768/R769), 0 reactions broken (full mass+charge).
#
# The 3 farnesyl-diphosphate copies MUST be changed together: R1116 is an FPP
# cross-membrane transport (m512[C_cy] <=> m610[C_mi]) that only balances when both
# ends carry the same charge.
#
# Verified charges (pH ~7.3 major microspecies, opened this session):
#   glyceraldehyde-3-phosphate(2-) C3H5O6P  -> -2  ChEBI:59776 (verified, formula matches)
#       https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:59776
#   farnesyl diphosphate(3-)       C15H25O7P2 -> -3  ChEBI:175763 (verified, formula matches)
#       https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:175763
#   ADP-ribose                     C15H21N5O14P2 -> -2  (InChI dH = -2, model-embedded InChI)
#   1-pyrroline-3-OH-5-carboxylate C5H6NO3       -> -1  (InChI dH = -1; sibling m753 already -1)
#
# NOT included: m966[C_va] vacuolar taurocholate — chemically -1 but R746 needs an
# H+ inserted first (see the r746-taurocholate-m966-needs-hplus note); deferred.
_CHARGE_STAGE2 = {
    "m536[C_cy]":  -2,   # D-glyceraldehyde 3-phosphate
    "m818[C_nu]":  -2,   # ADP-ribose (nucleus)
    "m1725[C_cy]": -2,   # ADP-ribose (cytosol)
    "m512[C_cy]":  -3,   # farnesyl diphosphate  ┐
    "m610[C_mi]":  -3,   # farnesyl diphosphate  ├ change together (R1116 symmetry)
    "m411[C_lp]":  -3,   # farnesyl diphosphate  ┘
    "m764[C_cy]":  -1,   # 1-pyrroline-3-hydroxy-5-carboxylate
}


def fix_charge_stage2(model) -> int:
    """Set 7 anion-stored metabolites (InChI dH < 0) to their dH charge (formula unchanged).

    Idempotent. Returns the number of metabolites changed. See _CHARGE_STAGE2 for sources.
    """
    changed = 0
    for mid, target in _CHARGE_STAGE2.items():
        try:
            met = model.metabolites.get_by_id(mid)
        except KeyError:
            continue
        if met.charge == target:
            continue
        old = met.charge
        met.charge = target
        changed += 1
        logger.info("  Charge Stage 2: %s (%s) charge %s -> %d"
                    % (mid, met.name, old, target))
    return changed


# ── Top-level driver ──────────────────────────────────────────────────────

def apply_all_patches(model) -> dict:
    """
    Apply all known model patches. Returns a dict with per-patch counts.
    Safe to call multiple times (idempotent — already-correct values are skipped).
    """
    logger.info("Applying iYali26 known-bug patches …")
    # NOTE: fix_ec_code_format is intentionally NOT called here.  Reaction EC
    # codes are populated later in the pipeline (gene EC enrichment, EC
    # backfill, reaction xref backfill), so EC formatting must run after those
    # steps — it is invoked separately near the end of main().
    counts = {
        "nadp_plus_fixed":   fix_nadp_plus_formula(model),
        "ceramide_fixed":    fix_ceramide_formulas(model),
        "cation_formula_fixed": fix_cation_formula_consistency(model),
        "d_arabinokinase_fixed": fix_d_arabinokinase_direction_and_proton(model),
    }
    logger.info(
        f"  patches applied: NADP+={counts['nadp_plus_fixed']} copies, "
        f"ceramide={counts['ceramide_fixed']} copies, "
        f"cation-formula={counts['cation_formula_fixed']} copies, "
        f"D-arabinokinase={counts['d_arabinokinase_fixed']} reaction"
    )
    return counts


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
_R1521_HANDOFF_CONTRACT = "18d60ef4e583311fd95581dcce2a96f4af6ac6c1d0dc19bc36799c163334d678"
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
        else (
            _coa_protonation_repository_root()
            / ("model.xml" if curation["source_model"] == _COA_PROTONATION_SOURCE_MODEL else curation["source_model"])
        ).resolve()
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
    contract = dict(handoff)
    contract.pop("target_contract_sha256", None)
    actual_contract = hashlib.sha256(
        json.dumps(
            contract,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode()
    ).hexdigest()
    if (
        handoff.get("target_contract_sha256") != _R1521_HANDOFF_CONTRACT
        or actual_contract != _R1521_HANDOFF_CONTRACT
        or handoff.get("source_sha256") != curation.get("source_sha256")
        or handoff.get("source_model_fingerprint") != curation.get("source_model_fingerprint")
    ):
        raise CoAProtonationCurationError("R1521 handoff source contract drifted")
    dependencies = handoff.get("evidence_dependencies")
    if not isinstance(dependencies, dict):
        raise CoAProtonationCurationError("R1521 evidence dependency contract drifted")
    root = _coa_protonation_repository_root()
    for name, dependency in dependencies.items():
        try:
            evidence = json.loads(
                (root / dependency["path"]).read_text(encoding="utf-8")
            )
            actual_digest = hashlib.sha256(
                json.dumps(
                    evidence,
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=True,
                ).encode()
            ).hexdigest()
        except (KeyError, OSError, TypeError, json.JSONDecodeError) as exc:
            raise CoAProtonationCurationError(
                f"cannot read R1521 evidence dependency {name}: {exc}"
            ) from exc
        if actual_digest != dependency.get("canonical_contract_sha256"):
            raise CoAProtonationCurationError(
                f"R1521 evidence dependency drifted: {name}"
            )
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
