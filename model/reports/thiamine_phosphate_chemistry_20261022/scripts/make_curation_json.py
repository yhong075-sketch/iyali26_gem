"""Write model/curation/chemistry/thiamine_phosphate_chemistry.json from before.json plus the reviewed values.

Before-values (formula, charge, annotation, notes, reaction lists, stoichiometry) are copied from before.json
so the guard compares exact values; the after-values and evidence are written here by hand.
Usage (from the worktree root): python make_curation_json.py BEFORE.json OUT.json
"""
import json
import sys

before_path, out = sys.argv[1:3]
before = json.load(open(before_path))
assert before["sha256"] == "8ae6c93b0ec41cde3db206f50ee9ad2549b1f96a3114d4f9df96a0d75bf01a76"

KEGG = "KEGG REST get/cpd:{kegg} (read 2026-10-22) gives {formula}"
# One row per molecule: KEGG fully protonated formula, ChEBI charge of that form, MetaNetX 4.5 ID, BiGG ID.
MOLECULES = {
    "thiamine monophosphate": {
        "species": ["m977[C_cy]", "m973[C_ex]"], "formula": "C12H18N4O4PS", "charge": 1,
        "kegg": "C01081", "chebi": "CHEBI:9533", "mnx": "MNXM662", "bigg": "thmmp",
        "chebi_name": "thiamine(1+) monophosphate", "ph73": "MNXM662 C12H16N4O4PS -1 (CHEBI:37575)",
        "balance": "R1337, R751, R754 and R1342 balance as written",
        "previous": ("formula C6H4Cl3NO, charge 0; annotated as 3,5,6-trichloro-2-methoxypyridine (MetaNetX "
                     "MNXM1137754, envipath, InChIKey RLIVUWLXZBDMBL), matched by the name TMP"),
    },
    "thiamine diphosphate": {
        "species": ["m975[C_cy]", "m971[C_ex]", "m1116[C_mi]"], "formula": "C12H19N4O7P2S", "charge": 1,
        "kegg": "C00068", "chebi": "CHEBI:9532", "mnx": "MNXM256", "bigg": "thmpp",
        "chebi_name": "thiamine(1+) diphosphate", "ph73": "MNXM256 C12H16N4O7P2S -2 (CHEBI:58937)",
        "balance": ("R1342, R753 and R933 balance as written; R752 and R750 keep the one H+ of their pH 7.3 "
                    "equations"),
        "previous": "no formula, charge 0, no cross-reference",
    },
    "thiamine triphosphate": {
        "species": ["m976[C_cy]"], "formula": "C12H20N4O10P3S", "charge": 1,
        "kegg": "C03028", "chebi": "CHEBI:9534", "mnx": "MNXM1341", "bigg": "thmtp",
        "chebi_name": "thiamine(1+) triphosphate", "ph73": "MNXM1341 C12H16N4O10P3S -3 (CHEBI:58938)",
        "balance": "R753 balances as written",
        "previous": "no formula, charge 0; annotated as BiGG ttp and MetaNetX MNXM148527 (Ttp, no formula)",
    },
    "4-methyl-5-(2-phosphonooxyethyl)thiazole": {
        "species": ["m643[C_cy]"], "formula": "C6H10NO4PS", "charge": 0,
        "kegg": "C04327", "chebi": "CHEBI:17857", "mnx": "MNXM960", "bigg": "4mpetz",
        "chebi_name": "4-methyl-5-(2-phosphonooxyethyl)thiazole", "ph73": "MNXM960 C6H8NO4PS -2 (CHEBI:58296)",
        "balance": ("R754 balances as written; R408 keeps the one H+ of its pH 7.3 equation; R1343 and R1344 "
                    "stay unbalanced in every protonation state tried, a 2-H discrepancy in their lumped "
                    "equations"),
        "previous": "no formula, charge 0, no cross-reference",
    },
}
CONVENTION = ("KEGG fully protonated form, matching the still-neutral ATP, ADP, AMP, phosphate, diphosphate "
              "and HMP-PP of this model: {balance}. With every partner at its pH 7.3 form all nine thiamine "
              "reactions balance as written, so this species must move to its pH 7.3 form together with them.")

species = {}
for molecule, row in MOLECULES.items():
    for mid in row["species"]:
        b = before["species"][mid]
        species[mid] = {
            "molecule": molecule, "name": b["name"], "compartment": b["compartment"],
            "reactions_before": b["reactions"], "notes_before": b["notes"],
            "before": {"formula": b["formula"], "charge": b["charge"], "annotation": b["annotation"]},
            "after": {"formula": row["formula"], "charge": row["charge"],
                      "annotation": {"sbo": "SBO:0000247", "kegg.compound": row["kegg"], "chebi": row["chebi"],
                                     "metanetx.chemical": row["mnx"], "bigg.metabolite": row["bigg"]}},
            "evidence": (KEGG.format(kegg=row["kegg"], formula=row["formula"])
                         + f"; ChEBI OLS4 {row['chebi']} {row['chebi_name']}, charge {row['charge']}; "
                         + f"MetaNetX 4.5 maps KEGG {row['kegg']} to {row['mnx']} (pH 7.3 form {row['ph73']})"),
            "evidence_tier": ("formula, charge and the KEGG, ChEBI and MetaNetX IDs verified on 2026-10-22; the "
                              "BiGG ID only through MetaNetX cross-references. The identity is inferred from the "
                              "reactions that use the species, since the source model gives no formula or ID: "
                              "R752 is KEGG R00619, R753 R00616, R1342 R00617, R1337 and R751 R02135, R408 "
                              "R04448 (equations read from KEGG on 2026-10-22)."),
            "rationale": CONVENTION.format(balance=row["balance"]),
            "previous_state": row["previous"],
        }

r1342 = before["reactions"]["R1342"]
spec = {
    "schema_version": 1,
    "curation_id": "THIAMINE-P-CHEMISTRY-20261022",
    "status": "user_authorized_formula_charge_annotation_correction",
    "authorization": ("2026-10-22 user: Approved: implement the thiamine phosphate chemistry plan as written: KEGG "
                      "fully protonated forms, m971 and m643 included, R1342 EC 2.7.4.9 removed, names unchanged, "
                      "on top of 8ae6c93b with --thiamine-phosphate-chemistry requiring --tpp-requirement, no H+ "
                      "edits, WT spot check only. No commit or push"),
    "baseline_model_sha256": before["sha256"],
    "report": "model/reports/thiamine_phosphate_chemistry_20261022/",
    "species": species,
    # Every reaction on these species, pinned so the balance review stays valid; biomass_C only by its TPP term.
    "context_reactions": {rid: row["stoichiometry"] for rid, row in before["reactions"].items()},
    # Formula and charge of every other species on these reactions; the balance review assumed them.
    "partner_species": before["partners"],
    "expected_balanced_after": ["R1337", "R1342", "R751", "R753", "R754", "R933"],
    "reactions": {
        "R1342": {
            "name": r1342["name"], "notes_before": r1342["notes"],
            "annotation_before": r1342["annotation"],
            "annotation_after": {**r1342["annotation"], "ec-code": "2.7.4.16"},
            "evidence": ("MetaNetX 4.5 reac_prop lists MNXR104907 with EC 2.7.4.16;2.7.4.9, while all of its "
                         "cross-references are thiamine-phosphate kinase (Rhea 15913, KEGG R00617, MetaCyc "
                         "THI-P-KIN-RXN, BiGG TMPK, SEED rxn00438 and rxn38718, SABIO 1007). ExPASy enzyme.dat "
                         "(release of 28-Jan-2026): 2.7.4.9 is dTMP kinase, 2.7.4.16 is thiamine-phosphate kinase."),
            "evidence_tier": "verified on 2026-10-22 from the local MetaNetX 4.5 and ExPASy files",
            "rationale": ("EC 2.7.4.9 belongs to thymidylate (dTMP) kinase; R1342 phosphorylates thiamine "
                          "monophosphate, so only EC 2.7.4.16 fits. Every other cross-reference is kept."),
            "previous_state": "ec-code 2.7.4.16 and 2.7.4.9",
        },
    },
}
# libSBML drops a note holding < or > and escapes ' and " on the next read, so new note texts avoid them.
texts = [spec["authorization"]] + [row[k] for row in [*species.values(), *spec["reactions"].values()]
                                   for k in ("evidence", "evidence_tier", "rationale", "previous_state")]
for text in texts:
    assert not set(text) & set("<>&'\""), text
with open(out, "w") as handle:
    json.dump(spec, handle, indent=1)
    handle.write("\n")
print("wrote", out, "species", len(species), "context reactions", len(spec["context_reactions"]))
