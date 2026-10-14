"""Diagnostic 'holistic package' overlay from the prompt header, in memory only (not curated, not in any XML).

Copied verbatim (function body unchanged) from Task A's
model/reports/w29_glucose_uptake_20261012/run_task_a.py, SHA-256 f28cef1670fc07ad8e6a2a66616c4db693fe1635cd660ba5a308ae5f367a6392,
uncommitted in worktree competent-turing-bc6e9b. R1384 <= 8 and R860 <= 1 are ASSUMPTION (diagnostic) caps.
Task A audit F3: with this model's NAD/NADH formulas the R1864/R1866 rewrite leaves +1 H and +1 charge,
and D-glycerate m1949 has no formula, so overlay numbers are diagnostic only.
"""

FAO = ("R1882", "R1883", "R1884", "R2255", "R2259", "R2263")


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
    assert {m.id: c for m, c in r.metabolites.items()} == {
        "m123[C_cy]": -1.0, "m276[C_cy]": -1.0, "m122[C_cy]": 1.0, "m277[C_cy]": 1.0}
    r.subtract_metabolites(dict(r.metabolites))
    # KEGG R01752: D-glyceraldehyde + NAD+ + H2O <=> D-glycerate + NADH + H+ (bounds kept).
    r.add_metabolites({met("m276[C_cy]"): -1, met("m122[C_cy]"): -1, met("m32[C_cy]"): -1,
                       met("m1949[C_cy]"): 1, met("m123[C_cy]"): 1, met("m10[C_cy]"): 1})
    record["R2069"] = r.reaction
    record["m1949_other_reactions"] = sorted(x.id for x in met("m1949[C_cy]").reactions if x.id != "R2069")
    for rid in FAO + ("R_MDHy", "R_ALCDH_nadp_hi"):
        model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    model.reactions.R1113.upper_bound = 0.0
    r171 = model.reactions.R171
    core = r171.gene_reaction_rule.split(") or ")[0].lstrip("(")
    assert r171.gene_reaction_rule == f"({core}) or YALI1B28025g or YALI1B15733g"
    r171.gene_reaction_rule = core
    model.reactions.R1384.upper_bound = 8.0  # ASSUMPTION (diagnostic): proton-leak cap
    model.reactions.R860.upper_bound = 1.0   # ASSUMPTION (diagnostic): Ggc1 GTP/GDP carrier cap
    record.update(closed=list(FAO) + ["R_MDHy", "R_ALCDH_nadp_hi"], R1113_upper=0.0,
                  R171_gpr=core, R1384_upper=8.0, R860_upper=1.0)
    return record
