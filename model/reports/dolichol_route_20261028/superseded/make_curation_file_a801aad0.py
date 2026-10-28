"""Write model/curation/dolichol/dolichol_route.json from before.json (the ce9e1112 before-state).

Usage (from the worktree root): python -I make_curation_file.py before.json OUT.json
Every pinned field (R2178 and R287 definitions, species identities and reaction contexts, the biomass_C
before-state) comes from before.json; the decision, the new reaction, the coefficient, the evidence and the
notes are the TASK.md values (Option B, approved 2026-10-28). Note texts avoid <, > and & (SBML notes).
"""
import json
import sys

before_path, out = sys.argv[1:3]
before = json.load(open(before_path))
T, S = before["targets"], before["species"]
assert before["sha256"] == "ce9e1112ee940531469457009234a2fd11de7f589c5948a7dec6b91118005294", before["sha256"]
assert before["new_reaction_id_absent"] and before["objective"] == {"biomass_C": 1.0}

CURATION_ID = "DOLICHOL-ROUTE-20261028"
NEW_ID = "R_DOLSYN"
COEFFICIENT = 1e-4
# Net of: polyprenyl-PP + H2O = polyprenol + PPi (uncharacterised phosphatase); polyprenol + NAD = polyprenal + NADH
# (Tda5); polyprenal + NADPH + H+ = dolichal + NADP (Dfg10); dolichal + NADPH + H+ = dolichol + NADP (Env9). Written in
# the model's own convention (NAD C21H27 / NADH C21H29 dehydrogenases carry no H+, NADPH reductases one H+ each, e.g.
# R534, R16; scripts/proton_convention_scan.py), hence 2 H+; the textbook notation with charged NAD+ / NADP+ has 1 H+ (TASK.md section 10).
NEW_STOICHIOMETRY = {"m2019[C_cy]": -1.0, "m32[C_cy]": -1.0, "m122[C_cy]": -1.0, "m40[C_cy]": -2.0, "m10[C_cy]": -2.0,
                     "m499[C_cy]": 1.0, "m203[C_cy]": 1.0, "m123[C_cy]": 1.0, "m41[C_cy]": 2.0}
SOURCES = {
    "sato2001": "https://europepmc.org/abstract/MED/11442630",
    "hoffmann2017": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5714454",
    "sato1999": "https://pmc.ncbi.nlm.nih.gov/articles/PMC83905",
    "cantagrel2010": "https://europepmc.org/abstract/MED/20637498",
    "hanaoka2026": "https://europepmc.org/abstract/MED/42201967",
    "gryz2019": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6628320",
    "kegg_sec59": "https://www.kegg.jp/entry/yli:2911738",
    "kegg_dfg10": "https://www.kegg.jp/entry/yli:2912457",
    "kegg_ec_1.3.1.94": "https://www.kegg.jp/entry/ec:1.3.1.94",
    "kegg_rer2": "https://www.kegg.jp/entry/yli:2909965",
    "sgd_sec59": "https://www.yeastgenome.org/locus/YMR013C",
}


def identity(mid):
    row = S[mid]
    return {key: row[key] for key in ("name", "formula", "charge", "compartment")}


def pinned(rid):
    row = T[rid]
    return {"reaction_name": row["name"], "gpr": row["gpr"], "bounds": row["bounds"],
            "stoichiometry": row["stoichiometry"], "annotation": row["annotation"], "groups": row["groups"]}


species = {}
for mid in sorted(set(NEW_STOICHIOMETRY) | {"m503[C_er]"}):
    species[mid] = identity(mid)
# The route species whose reaction context is part of the decision: dehydrodolichol-PP (made by R2180, used by the
# R2178 being closed), dolichol (R287 its only consumer) and Dol-P in the ER (the biomass term's species).
for mid in ("m2019[C_cy]", "m499[C_cy]", "m503[C_er]"):
    species[mid]["reactions_before"] = S[mid]["reactions"]
assert S["m2019[C_cy]"]["reactions"] == ["R2178", "R2180"] and S["m499[C_cy]"]["reactions"] == ["R287"]
assert "biomass_C" not in S["m503[C_er]"]["reactions"]
assert T["R2178"]["bounds"] == [0.0, 1000.0] and T["R2178"]["gpr"] == "" and T["R2178"]["notes"] == {}
assert "m503[C_er]" not in T["biomass_C"]["stoichiometry"]

authorization = ("2026-10-28 user, after reading the Stage 1 draft (model/reports/dolichol_route_20261028/TASK.md): "
                 "\"Approve Option B\" = close R2178, add the lumped gene-less dehydrodolichol-PP to dolichol step with the "
                 "three-step-detour net stoichiometry, and add dolichyl phosphate m503[C_er] to biomass_C at 1e-4 mmol/gDW, "
                 "through a guarded --dolichol-route step on ce9e1112. No commit or push.")
spec = {
    "schema_version": 1,
    "curation_id": CURATION_ID,
    "status": "user_authorized",
    "authorization": authorization,
    "baseline_model_sha256": before["sha256"],
    "report": "model/reports/dolichol_route_20261028/",
    "reviewer": "Claude Code session (Stage 1 investigation, 2026-10-28); decision by the user; independent audit in the report",
    "requires_objective": "biomass_C",
    "species": species,
    "closed_reaction": {
        "reaction_id": "R2178", **pinned("R2178"), "after_bounds": [0.0, 0.0],
        "reasons": [
            "No gene, no name (NO NAME), no annotation beyond SBO; added with the iYli21 R2130-R2200 block and un-annotated in the raw iyli21.xml.",
            "Its equation (6.5 O2 + dehydrodolichol-PP = dolichyl-PP + 5 CO2 + 3 H2O) matches no enzyme of the yeast dolichol route: polyprenyl-PP is dephosphorylated and reduced to dolichol, which Sec59 phosphorylates (TASK.md B7-B12).",
            "Kept as a closed reaction rather than deleted: deletion needs evidence that it does not exist, and closing records the decision without losing the inherited entry.",
            "With a Dol-P drain and R2178 open, R2178 is obligatory (FVA min = max) and the CWH8 ortholog becomes essential while the SEC59 ortholog stays dispensable, the reverse of S. cerevisiae genetics (TASK.md section 4, Option A).",
        ],
        "notes": {
            "dolichol_route_curation_id": CURATION_ID,
            "dolichol_route_previous_bounds": "[0, 1000]",
            "dolichol_route_r2178_reason": "Closed (not deleted): gene-less, un-annotated iYli21 placeholder whose O2-consuming, CO2-releasing equation matches no enzyme of the dolichol route; replaced by R_DOLSYN (dehydrodolichol-PP to dolichol, then Sec59 R287 to Dol-P).",
            "dolichol_route_authorization": "User authorized Option B on 2026-10-28 (model/reports/dolichol_route_20261028/TASK.md).",
        },
    },
    "kept_reaction": {"reaction_id": "R287", **pinned("R287"),
                      "role": "dolichol kinase (SEC59 ortholog YALI1E22718g, KEGG K00902); the only consumer of dolichol m499 and the step that the new reaction feeds"},
    "new_reaction": {
        "reaction_id": NEW_ID,
        "reaction_name": "dehydrodolichyl diphosphate to dolichol (lumped dephosphorylation, polyprenol oxidation, polyprenal and dolichal reduction)",
        "bounds": [0.0, 1000.0],
        "gpr": "",
        "group": "N-Glycan biosynthesis",
        "stoichiometry": NEW_STOICHIOMETRY,
        "annotation": {"sbo": "SBO:0000176"},
        "component_steps": [
            {"step": "polyprenyl-PP + H2O = polyprenol + PPi", "enzyme": "not characterised in the sources read (TASK.md B12)"},
            {"step": "polyprenol + NAD = polyprenal + NADH (model convention; NAD+ + H+ release in textbook notation)", "enzyme": "S. cerevisiae Tda5 (Hanaoka et al. 2026); Y. lipolytica ortholog not identified"},
            {"step": "polyprenal + NADPH + H+ = dolichal + NADP+", "enzyme": "S. cerevisiae Dfg10, EC 1.3.1.94, KEGG R13444; Y. lipolytica candidate ortholog YALI1E09357g (YALI0E07788g, K12345), not in the model"},
            {"step": "dolichal + NADPH + H+ = dolichol + NADP+", "enzyme": "S. cerevisiae Env9 (Hanaoka et al. 2026); Y. lipolytica ortholog not identified"},
        ],
        "why_gene_less": "Three enzymes, one of them uncharacterised and two without an identified Y. lipolytica ortholog; a single-gene GPR on the lumped step would not meet AGENTS.md section 5, and dfg10 deletion is viable in S. cerevisiae (alternative route, Cantagrel 2010).",
        "balance_note": ("Cannot be balance-checked as built: dehydrodolichol-PP m2019 has no formula and dolichol m499 carries a farnesol formula C15H28O. With a Dol-16 formula set filled in memory (m2019 C80H132O7P2, m499 C80H132O) C, H, N, O and P balance and the charge is -2, the model's NADP charge convention (NADP C21H29 charge 0, -1 per NADPH; the same pattern as R16 and R1393). The dolichol-family formula set is a separate chemistry task."),
        "notes": {
            "dolichol_route_curation_id": CURATION_ID,
            "dolichol_route_status": "user_authorized",
            "dolichol_route_evidence": "S. cerevisiae: polyprenol is reduced to dolichol through Tda5 (NAD+), Dfg10 and Env9 (NADPH) (Hanaoka et al. 2026, PMID 42201967); SRD5A3/Dfg10 needed for dolichol (PMID 20637498); Sec59 phosphorylates dolichol. Y. lipolytica: DFG10 candidate ortholog YALI1E09357g by KEGG orthology only.",
            "dolichol_route_evidence_tier": "S. cerevisiae experiment plus Y. lipolytica orthology; no Y. lipolytica dolichol data",
            "dolichol_route_limits": "Lumped and gene-less; not mass-balanced (m2019 has no formula, m499 has a farnesol formula); the polyprenyl-PP phosphatase is uncharacterised; Tda5 and Env9 orthologs not identified; NUS1 ortholog YALI1F26726g not in the model.",
            "dolichol_route_sources": "; ".join(SOURCES[k] for k in ("hanaoka2026", "cantagrel2010", "kegg_ec_1.3.1.94", "kegg_dfg10", "kegg_sec59")),
            "dolichol_route_authorization": "User authorized Option B on 2026-10-28 (model/reports/dolichol_route_20261028/TASK.md).",
        },
    },
    "biomass_term": {
        "biomass_reaction": "biomass_C",
        "metabolite_id": "m503[C_er]",
        "coefficient_mmol_per_gDW": COEFFICIENT,
        "derivation": ("Order-of-magnitude value. S. cerevisiae dolichol 32 (late log) to 48 (early stationary) ug/g, basis not "
                       "stated and probably wet weight (Gryz et al. 2019, citing Surmacz and Swiezewska 2011, unread); Dol-16 "
                       "C80H132O 1109.9 g/mol; with dry matter 20-25 % of wet weight this is 1.2e-4 to 2.2e-4 mmol/gDW, and "
                       "2.9e-5 to 4.3e-5 if the figure were dry weight. The term stands for the dolichol-derived carrier pool, "
                       "since free dolichol has no consumer other than R287. Stage 1: growth differs by less than 1e-4 relative "
                       "across 4e-5 to 1.6e-4 and no targeted knockout changes with the value."),
        "source": SOURCES["gryz2019"],
        "evidence_tier": "literature value for S. cerevisiae, unit basis assumed; not a Y. lipolytica measurement",
        "accepted_before_state": {
            "id": "ce9e1112_biomass_C",
            "description": "biomass_C of candidate ce9e1112 (2026-10-28): a1237821 plus the six cofactor terms",
            "biomass_stoichiometry": T["biomass_C"]["stoichiometry"],
            "biomass_bounds": T["biomass_C"]["bounds"],
            "biomass_note_keys": sorted(T["biomass_C"]["notes"]),
        },
        "notes": {
            "dolichol_route_curation_id": CURATION_ID,
            "dolichol_route_biomass_term": f"m503[C_er] {COEFFICIENT:g} mmol/gDW",
            "dolichol_route_biomass_evidence": "S. cerevisiae dolichol 32-48 ug/g (Gryz et al. 2019, PMC6628320; basis probably wet weight); order-of-magnitude value",
            "dolichol_route_biomass_caveat": "Saccharomyces cerevisiae literature value with an assumed wet-to-dry conversion, not a Yarrowia lipolytica measurement",
        },
    },
    "evidence": [
        {"url": SOURCES["kegg_sec59"], "claim": "YALI1E22718g = YALI0E18942g, KEGG yli:2911738, dolichol kinase K00902 (EC 2.7.1.108); on the consensus essential list.", "tier": "verified", "scope": "Orthology, not experiment.", "accessed": "2026-10-28"},
        {"url": SOURCES["kegg_rer2"], "claim": "YALI1C26036g = YALI0C18799g, RER2-type dehydrodolichyl diphosphate synthase K11778.", "tier": "verified", "scope": "Orthology.", "accessed": "2026-10-28"},
        {"url": SOURCES["kegg_dfg10"], "claim": "YALI1E09357g = YALI0E07788g, polyprenol reductase K12345 (DFG10 / SRD5A3), not in the model.", "tier": "verified", "scope": "Orthology; recorded as a note, not a GPR.", "accessed": "2026-10-28"},
        {"url": SOURCES["kegg_ec_1.3.1.94"], "claim": "EC 1.3.1.94 polyprenal reductase, KEGG R13444 (dolichal + NADP+ = polyprenal + NADPH + H+); lists sce YIL049W and yli 2912457.", "tier": "verified", "scope": "Database entry.", "accessed": "2026-10-28"},
        {"url": SOURCES["sato2001"], "claim": "Srt1 products are longer, mainly on lipid particles, and not sufficiently converted to dolichol and Dol-P.", "tier": "verified", "scope": "Abstract read; S. cerevisiae; why UM6 is left as it is.", "accessed": "2026-10-28"},
        {"url": SOURCES["hoffmann2017"], "claim": "Rer2 is also on lipid droplets; long Srt1 polyprenols are not efficiently converted to dolichol.", "tier": "verified", "scope": "Article read; S. cerevisiae.", "accessed": "2026-10-28"},
        {"url": SOURCES["sato1999"], "claim": "rer2 deletion grows very slowly with about 2.7 % of wild-type cis-prenyltransferase activity; rer2 srt1 double disruption is lethal.", "tier": "verified", "scope": "Article read; S. cerevisiae.", "accessed": "2026-10-28"},
        {"url": SOURCES["sgd_sec59"], "claim": "sec59 deletion is inviable in S. cerevisiae.", "tier": "verified", "scope": "SGD phenotype record; S. cerevisiae.", "accessed": "2026-10-28"},
        {"url": SOURCES["cantagrel2010"], "claim": "SRD5A3 / Dfg10 are required to reduce polyprenol to dolichol; depleted cells still contain some dolichol.", "tier": "verified", "scope": "Abstract read.", "accessed": "2026-10-28"},
        {"url": SOURCES["hanaoka2026"], "claim": "In yeast, Tda5 oxidises polyprenol to polyprenal, Dfg10 reduces polyprenal to dolichal and Env9 reduces dolichal to dolichol (NADPH).", "tier": "verified", "scope": "Abstract read; S. cerevisiae; Y. lipolytica TDA5 / ENV9 orthologs not identified.", "accessed": "2026-10-28"},
        {"url": SOURCES["gryz2019"], "claim": "S. cerevisiae dolichol 32 ug/g (late log) and 48 ug/g (early stationary), citing Surmacz and Swiezewska 2011 (unread).", "tier": "verified", "scope": "Passage read; basis not stated in the sentence; primary measurement unread.", "accessed": "2026-10-28"},
    ],
}
with open(out, "w") as handle:
    json.dump(spec, handle, indent=1)
    handle.write("\n")
print(out, "species", len(species), "before-state biomass species", len(T["biomass_C"]["stoichiometry"]),
      "note keys", len(T["biomass_C"]["notes"]))
