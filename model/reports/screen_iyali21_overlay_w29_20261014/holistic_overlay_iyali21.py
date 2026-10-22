"""Task A's diagnostic 'holistic package' adapted to iYali21 (iyli21.xml), in memory only.

Copied from reports/r132_mito_adh_20261012/holistic_overlay.py (SHA-256 bf6c7ce5...), itself a verbatim
copy of Task A's apply_overlay. Two adaptations, approved by the user on 2026-10-14:
1. R2069 starts from iYali21's form (H+ + NADH + D-glyceraldehyde <=> NAD+ + m277 + H2O) instead of
   iYali26's; it is rewritten to the same KEGG R01752 end equation.
2. R_MDHy and R_ALCDH_nadp_hi do not exist in iYali21, so there is nothing to close (asserted absent).
Everything else is unchanged. R1384 <= 8 and R860 <= 1 are ASSUMPTION (diagnostic) caps.
Task A audit F3 applies: the R1864/R1866 rewrite leaves +1 H and +1 charge; m1949 has no formula.
"""

FAO = ("R1882", "R1883", "R1884", "R2255", "R2259", "R2263")
ABSENT_IN_IYALI21 = ("R_MDHy", "R_ALCDH_nadp_hi")


def apply_overlay(model):
    """Diagnostic 'holistic package' in memory only (not curated, not in the XML)."""
    met = model.metabolites.get_by_id
    record = {}
    for rid, alcohol, aldehyde in (("R1864", "m1789[C_cy]", "m1790[C_cy]"), ("R1866", "m1792[C_cy]", "m1793[C_cy]")):
        r = model.reactions.get_by_id(rid)
        assert {m.id: c for m, c in r.metabolites.items()} == {
            "m122[C_cy]": -2.0, alcohol: -1.0, "m32[C_cy]": -1.0, "m10[C_cy]": 2.0, "m123[C_cy]": 2.0, aldehyde: 1.0}
        r.subtract_metabolites(dict(r.metabolites))
        r.add_metabolites({met("m122[C_cy]"): -1, met(alcohol): -1, met("m10[C_cy]"): 1,
                           met("m123[C_cy]"): 1, met(aldehyde): 1})
        record[rid] = r.reaction
    r = model.reactions.R2069
    # Adaptation 1: iYali21's starting form.
    assert {m.id: c for m, c in r.metabolites.items()} == {
        "m10[C_cy]": -1.0, "m123[C_cy]": -1.0, "m276[C_cy]": -1.0, "m122[C_cy]": 1.0, "m277[C_cy]": 1.0, "m32[C_cy]": 1.0}
    record["R2069_before"] = r.reaction
    r.subtract_metabolites(dict(r.metabolites))
    # KEGG R01752: D-glyceraldehyde + NAD+ + H2O <=> D-glycerate + NADH + H+ (bounds kept).
    r.add_metabolites({met("m276[C_cy]"): -1, met("m122[C_cy]"): -1, met("m32[C_cy]"): -1,
                       met("m1949[C_cy]"): 1, met("m123[C_cy]"): 1, met("m10[C_cy]"): 1})
    record["R2069"] = r.reaction
    record["m1949_other_reactions"] = sorted(x.id for x in met("m1949[C_cy]").reactions if x.id != "R2069")
    for rid in FAO:
        model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    # Adaptation 2: these iYali26 reactions do not exist here.
    assert not any(rid in model.reactions for rid in ABSENT_IN_IYALI21)
    model.reactions.R1113.upper_bound = 0.0
    r171 = model.reactions.R171
    core = r171.gene_reaction_rule.split(") or ")[0].lstrip("(")
    assert r171.gene_reaction_rule == f"({core}) or YALI1B28025g or YALI1B15733g"
    r171.gene_reaction_rule = core
    model.reactions.R1384.upper_bound = 8.0  # ASSUMPTION (diagnostic): proton-leak cap
    model.reactions.R860.upper_bound = 1.0   # ASSUMPTION (diagnostic): Ggc1 GTP/GDP carrier cap
    record.update(closed=list(FAO), absent_not_closed=list(ABSENT_IN_IYALI21), R1113_upper=0.0,
                  R171_gpr=core, R1384_upper=8.0, R860_upper=1.0)
    return record
