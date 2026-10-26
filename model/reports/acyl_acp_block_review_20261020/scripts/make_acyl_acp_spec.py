"""Write model/curation/chemistry/acyl_acp_chemistry.json from the approved P1/P2/P3 design.

Run from the repository root with the lipoylation v2 candidate as input. Records the before state of
every touched species and reaction, applies the change in memory, and records the after state and
the residual of every block reaction. Asserts the intended residual pattern before writing.
"""
import json

import cobra

MODEL = "model/candidates/E5_coq9_resp_pkg_r132_combined_lipoylation_v2_20261020.xml"
OUT = "model/curation/chemistry/acyl_acp_chemistry.json"
M = cobra.io.read_sbml_model(MODEL)
g, rx = M.metabolites.get_by_id, M.reactions.get_by_id

FILLS = {  # acyl group + S + R, as for acetyl-ACP C2H3OSR, malonyl-ACP C3H3O3SR, octanoyl-ACP C8H15OSR, ACP1 HSR
    "m1410[C_mi]": "C4H5O2SR", "m1415[C_mi]": "C4H7O2SR", "m1418[C_mi]": "C4H5OSR", "m1411[C_mi]": "C4H7OSR",
    "m1412[C_mi]": "C6H9O2SR", "m1416[C_mi]": "C6H11O2SR", "m1419[C_mi]": "C6H9OSR", "m1413[C_mi]": "C6H11OSR",
    "m1414[C_mi]": "C8H13O2SR", "m1417[C_mi]": "C8H15O2SR", "m1420[C_mi]": "C8H13OSR",
    "m1901[C_mi]": "C10H19OSR", "m1910[C_cy]": "C10H19OSR",
    "m1902[C_mi]": "C12H23OSR", "m1911[C_cy]": "C12H23OSR",
    "m1903[C_mi]": "C14H27OSR", "m1912[C_cy]": "C14H27OSR",
    "m1904[C_mi]": "C14H25OSR", "m1913[C_cy]": "C14H25OSR",
    "m1905[C_mi]": "C16H31OSR", "m1914[C_cy]": "C16H31OSR",
    "m1906[C_mi]": "C16H29OSR", "m1915[C_cy]": "C16H29OSR",
    "m1907[C_mi]": "C18H35OSR", "m1916[C_cy]": "C18H35OSR",
    "m1908[C_mi]": "C18H33OSR", "m1917[C_cy]": "C18H33OSR",
    "m1909[C_mi]": "C18H31OSR", "m1918[C_cy]": "C18H31OSR",
    "m1919[C_cy]": "C8H15OSR", "m1920[C_cy]": "HSR", "m1922[C_cy]": "C18H32O2",
}
KS = ["R1394", "R1395", "R1396"]
ELONGATION = ["R2218", "R2221", "R2224", "R2230"]
DESATURATING = ["R2227", "R2233"]
HYDROLASES = ["R2217", "R2220", "R2223", "R2226", "R2229", "R2232", "R2235", "R2237", "R2239", "R2241"]
REVERSIBLE_HYDROLASES = ["R2226", "R2229", "R2232", "R2235", "R2237", "R2239", "R2241"]
H_MITO, H_CYTO = "m28[C_mi]", "m10[C_cy]"
BLOCK = ([f"R{i}" for i in range(1394, 1406)] + [f"R{i}" for i in range(2216, 2242)])

assert len(FILLS) == 32
assert len(KS + ELONGATION + DESATURATING + HYDROLASES) == 19


def state(r):
    return {"stoichiometry": {m.id: float(c) for m, c in r.metabolites.items()}, "bounds": [float(b) for b in r.bounds],
            "gpr": r.gene_reaction_rule}


def imbalance(r):
    try:
        b = r.check_mass_balance()
    except ValueError as e:
        return None
    return {k: v for k, v in b.items() if abs(v) > 1e-9}


fill_users = {mid: sorted(r.id for r in g(mid).reactions) for mid in FILLS}
assert all(set(u) <= set(BLOCK) for u in fill_users.values()), "a filled species is used outside the block"
participants = sorted({m.id for rid in BLOCK for m in rx(rid).metabolites} - set(FILLS))
species = {mid: {"formula": g(mid).formula, "charge": g(mid).charge, "compartment": g(mid).compartment} for mid in participants}
fills = {mid: {"before": {"formula": g(mid).formula, "charge": g(mid).charge},
               "after": {"formula": f, "charge": g(mid).charge}} for mid, f in FILLS.items()}
before = {rid: state(rx(rid)) for rid in BLOCK}
residual_before = {rid: imbalance(rx(rid)) for rid in BLOCK}

# apply in memory
for mid, f in FILLS.items():
    g(mid).formula = f
for rid in KS + ELONGATION + DESATURATING:
    rx(rid).add_metabolites({g(H_MITO): 1.0})        # one H+ substrate fewer
for rid in HYDROLASES:
    rx(rid).add_metabolites({g(H_CYTO): -1.0})       # one H+ product fewer
for rid in REVERSIBLE_HYDROLASES:
    rx(rid).lower_bound = 0.0
after = {rid: state(rx(rid)) for rid in BLOCK}
residual_after = {rid: imbalance(rx(rid)) for rid in BLOCK}

expected_pattern = {rid: {} for rid in BLOCK}
for rid in ELONGATION:
    expected_pattern[rid] = {"charge": -2.0}
for rid in DESATURATING:
    expected_pattern[rid] = {"charge": -3.0}
for rid in ["R1397", "R1398", "R1399", "R1403", "R1404", "R1405"]:
    expected_pattern[rid] = {"charge": -1.0}
assert all(residual_after[r] is not None and {k: round(v, 9) for k, v in residual_after[r].items()} == expected_pattern[r]
           for r in BLOCK), {r: residual_after[r] for r in BLOCK if residual_after[r] != expected_pattern[r]}

changed = sorted(r for r in BLOCK if before[r] != after[r])
spec = {
    "schema_version": 1, "change_id": "ACYL-ACP-CHEMISTRY-20261020", "status": "user_authorized_candidate",
    "reviewer": "human_user",
    "authorization": "2026-10-20: user approved P1 (32 formula fills), P2 (remove the 19 extra H+), P3 (7 hydrolases "
                     "forward-only after a thermodynamic citation) and P5 (lipoylation guard compatibility) as an opt-in "
                     "change; no transport-direction change, no deletions, no commit.",
    "report": "model/reports/acyl_acp_block_review_20261020/REPORT.md",
    "requires": {"reaction": "R2286", "species_formula": {"m179[C_mi]": "HSR"},
                 "reason": "the fills balance only with ACP1 = HSR, set by LIPOYLATION-MODULE-20261020"},
    "formula_fills": fills, "fill_users": fill_users, "species": species,
    "block_reactions": BLOCK, "changed_reactions": changed,
    "reactions_before": before, "reactions_after": after,
    "residuals_before": residual_before, "residuals_after": residual_after,
    "h_removed": {"mitochondrial_substrate": KS + ELONGATION + DESATURATING, "cytosolic_product": HYDROLASES},
    "made_irreversible": REVERSIBLE_HYDROLASES,
    "evidence": {
        "formulas": "Acyl group + S + R, the convention of acetyl-ACP C2H3OSR, malonyl-ACP C3H3O3SR, octanoyl-ACP C8H15OSR "
                    "and ACP1 HSR. Replaces placeholders 'ACP' (read as elements A, C, P), 'C141'/'C161' (read as 141/161 "
                    "carbons), formulas with '*', or none.",
        "protons": "The block copies BiGG iMM904/iND750 stoichiometry written for ionized species; with the model's neutral "
                   "acyl and fatty-acid species one H+ is surplus in each condensation, hydrolysis and elongation lump. "
                   "The remaining -1 charge per NADPH is the model-wide NADP(+) (m178) convention.",
        "hydrolase_direction": "Thioester hydrolysis is strongly exergonic: eQuilibrator estimates dG'0 = -29.6 +/- 12.6 kJ/mol "
                               "(K'eq 1.6e5) and dG'm = -46.8 kJ/mol for palmitoyl-CoA + H2O = palmitate + CoA (pH 7.5, "
                               "I 0.25 M), the analogous phosphopantetheine thioester. The reverse direction would make an "
                               "acyl thioester from a free fatty acid without ATP. Making the 7 reversible cytosolic "
                               "hydrolases forward-only (as R2217, R2220, R2223 already are) removes the C14->C16 cycle "
                               "observed through R2226(reverse)/R2225/R2230/R2231/R2232.",
        "sources": ["https://equilibrator.weizmann.ac.il/search?query=palmitoyl-CoA+%2B+H2O+%3D+palmitate+%2B+CoA",
                    "https://bigg.ucsd.edu/models/iMM904/reactions/FA180ACPH"],
    },
}
json.dump(spec, open(OUT, "w"), indent=1)
print("written", OUT, "| changed reactions:", len(changed), changed)
