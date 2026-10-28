"""Write model/curation/biomass/cofactor_biomass_terms.json from before.json (the a1237821 before-state).

Usage (from the worktree root): python -I make_curation_file.py before.json OUT.json
The species rows, the biomass_C before-state (stoichiometry, bounds, note keys) and the baseline SHA come
from before.json; the coefficients, sources and authorization are the TASK.md values.
"""
import json
import sys

before_path, out = sys.argv[1:3]
before = json.load(open(before_path))
SOURCE_URL = "https://github.com/SysBioChalmers/yeast-GEM/blob/main/model/yeast-GEM.yml"
SOURCE_SHA = "2ddfaaa73ead41ab243127d1b34fa56db2c86052f6eabcd12bb4eea93947dc17"
TERMS = {
    "m479[C_cy]": ("THF", 6.34e-5, "s_1487", "-6.34000025456771e-05",
                   "CHEBI:166569 and MNXM1128450 resolve to 6R-tetrahydrofolic acid, C19H23N7O6 (read 2026-10-27); KEGG C00101 formula C19H23N7O6"),
    "m528[C_cy]": ("FAD", 1e-5, "s_0687", "-9.99999974737875e-06",
                   "KEGG C00016 (FAD) formula C27H33N9O15P2 matches the model (read 2026-10-27)"),
    "m926[C_cy]": ("riboflavin", 9.9e-4, "s_1405", "-0.000989999971352518",
                   "KEGG C00255 (riboflavin) formula C17H20N4O6 matches the model (read 2026-10-27)"),
    "m122[C_cy]": ("NAD", 2.65e-3, "s_1198", "-0.00264999992214143",
                   "CHEBI:13389 and MNXM64096 are the NAD class entry with unspecified oxidation state (no structure); the species is NAD+ by its reactions (R534 malate dehydrogenase and the other NAD-reducing steps) and its neutral formula C21H27N7O14P2 is KEGG C00003 (C21H28N7O14P2, the cation) less one proton; cross-reference left unchanged"),
    "m41[C_cy]": ("NADP", 5.7e-4, "s_1207", "-0.000569999974686652",
                  "KEGG C00006 (NADP+) formula C21H29N7O17P3 matches the model formula field; the model name says C21H28 (finding, not changed)"),
    "m69[C_cy]": ("CoA", 1.9e-4, "s_0529", "-0.000190000006114133",
                  "CHEBI:23354 and MNXM1094981 are the generic class coenzyme, not coenzyme A (finding, not changed); the species is coenzyme A by its name, its formula C21H36N7O16P3S (KEGG C00010) and its reactions (dephospho-CoA kinase R1092 product, acyl-CoA chemistry)"),
}
spec = {
    "schema_version": 1,
    "curation_id": "COFACTOR-BIOMASS-TERMS-20261027",
    "status": "user_authorized_biomass_terms",
    "evidence_tier": "consensus_model_value",
    "authorization": ("2026-10-27 user: Approved: implement the cofactor biomass-terms plan as written: --cofactor-biomass-terms "
                      "on top of a1237821 (running after --acyl-acp-chemistry), THF m479, FAD m528, riboflavin m926, NAD m122, "
                      "NADP m41 and CoA m69 in biomass_C at the yeast-GEM 9.1.1 r_4598 values 6.34e-5, 1e-5, 9.9e-4, 2.65e-3, 5.7e"),
    "authorization_note": ("The message as received ends mid-sentence at the app reply-button limit of 300 characters; the plan it "
                           "approves as written is model/reports/cofactor_biomass_terms_20261027/TASK.md (CoA 1.9e-4 mmol/gDW, "
                           "phase-1 checks, builds, diff, knockout sweep, re-census, audit, no commit or push)."),
    "baseline_model_sha256": before["sha256"],
    "report": "model/reports/cofactor_biomass_terms_20261027/",
    "reviewer": "Claude (this session, 2026-10-27); independent audit in the report",
    "biomass_reaction": "biomass_C",
    "requires_objective": "biomass_C",
    "source": {"model": "yeast-GEM", "version": "9.1.1", "reaction": "r_4598 cofactor pseudoreaction", "url": SOURCE_URL,
               "raw_file_sha256": SOURCE_SHA, "read": "2026-10-27 (same bytes as the 2026-10-22 read for TPP-REQUIREMENT-20261022)"},
    "organism_caveat": "Saccharomyces cerevisiae consensus-model values, not Yarrowia lipolytica measurements",
    "terms": {},
    "accepted_before_state": {
        "id": "a1237821_biomass_C",
        "description": "biomass_C of candidate a1237821 (2026-10-23): base plus the reviewed TPP, heme a / iron and lipoyl terms",
        "biomass_stoichiometry": before["biomass_C"]["stoichiometry"],
        "biomass_bounds": before["biomass_C"]["bounds"],
        "biomass_note_keys": sorted(before["biomass_C"]["notes"]),
    },
}
for mid, (label, coefficient, species, exact, identity) in TERMS.items():
    row = before["species"][mid]
    assert row["compartment"] == "C_cy", mid
    spec["terms"][mid] = {
        "label": label, "name": row["name"], "formula": row["formula"], "charge": row["charge"], "compartment": row["compartment"],
        "reactions_before": [r for r in row["reactions"] if r != "biomass_C"],
        "coefficient_mmol_per_gDW": coefficient,
        "source": f"yeast-GEM v9.1.1 r_4598 cofactor pseudoreaction, {label} [c] {species} = {exact} mmol/gDW",
        "identity_check": identity,
    }
assert not set(spec["terms"]) & set(spec["accepted_before_state"]["biomass_stoichiometry"])
with open(out, "w") as handle:
    json.dump(spec, handle, indent=1)
    handle.write("\n")
print(out, len(spec["terms"]), "terms; before-state metabolites", len(spec["accepted_before_state"]["biomass_stoichiometry"]),
      "note keys", spec["accepted_before_state"]["biomass_note_keys"])
