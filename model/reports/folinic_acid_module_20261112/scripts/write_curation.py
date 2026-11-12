"""Write model/curation/folate/folinic_acid_module.json from before.json (the e584151f state) and the decisions."""
import copy
import json
import sys
from pathlib import Path

before = json.loads(Path(sys.argv[1]).read_text())
out_path = Path(sys.argv[2])
R = before["reactions"]
SHMT = "YALI1D28569g or YALI1E19557g"
CID = "FOLINIC-ACID-MODULE-20261112"
TIER = ("Chemistry: Rhea 34767 (reviewed enzymes are SHMTs only) and Stover and Schirch 1990 (PMID 2201683, "
        "abstract read). Yeast source: Holmes and Appling 2002 (PMID 11923304, abstract read), 5-formyl-THF "
        "undetectable without SHMT. Y. lipolytica genes: KEGG orthology K00600 only; paralog localisation not checked; "
        "the side reaction takes the same rule as the SHMT step in its compartment.")


def pinned(row):
    return {k: copy.deepcopy(row[k]) for k in ("name", "bounds", "gpr", "stoichiometry", "species", "annotation",
                                               "notes", "groups")}


r748_after = pinned(R["R748"])
r748_after.update(
    name="5,10-methenyl-THF hydrolysis (SHMT side activity)",
    gpr=SHMT,
    annotation={**R["R748"]["annotation"], "ec-code": "2.1.2.1"},
    notes={"PROTEIN_CLASS": "2.1.2.1",
           "folinic_acid_module_curation": f"{CID}: rule moved from the glycine-cleavage T protein YALI1F04082g "
                                           "(K00605, EC 2.1.2.10) to the SHMT genes, as on R359; name and EC follow.",
           "folinic_acid_module_evidence": "KEGG R02300 is filed under EC 2.1.2.10, which ExPASy ENZYME describes "
                                           "only as the lipoyl-aminomethyl to methylene-THF step; the reaction is "
                                           "the SHMT side activity (Rhea 34767). Yeast-GEM r_1030 and iMM904 THFATm "
                                           "also carry glycine-cleavage-complex rules (different subunit sets); "
                                           "that iYali26 inherited its rule from them is inferred.",
           "folinic_acid_module_evidence_tier": TIER})

r990_after = pinned(R["R990"])
r990_after.update(
    bounds=[0.0, 0.0],
    notes={"folinic_acid_module_curation": f"{CID}: closed, kept in the model. The stoichiometry is exactly R64 plus "
                                           "R543 (5-formyl-THF cyclo-ligase then methenyl-THF cyclohydrolase), "
                                           "which the model already has with genes.",
           "folinic_acid_module_evidence": "No EC, KEGG or Rhea entry describes a 5-formyl-THF isomerase (searched "
                                           "2026-11-12); BiGG FTHFI is in iMM904 only, gene-less. Inferred: open "
                                           "and gene-less it would bypass a knockout of the FAU1-type gene "
                                           "YALI1D01489g once cytosolic 5-formyl-THF has a source (R_THFAT).",
           "folinic_acid_module_evidence_tier": "Lumped duplicate: verified by stoichiometry on e584151f. Absence of a "
                                                "distinct isomerase: no source found (KEGG, Rhea, BiGG checked)."})

new = {
    "reaction_id": "R_THFAT",
    "name": "5,10-methenyl-THF hydrolysis (SHMT side activity), cytosol",
    "subsystem": "One carbon pool by folate",
    "bounds": [0.0, 1000.0],
    "gpr": SHMT,
    "stoichiometry": {"m160[C_cy]": -1.0, "m32[C_cy]": -1.0, "m159[C_cy]": 1.0, "m10[C_cy]": 1.0},
    "annotation": {"sbo": "SBO:0000176", "bigg.reaction": "THFAT", "ec-code": "2.1.2.1", "rhea": "34767",
                   "kegg.reaction": "R02300", "metanetx.reaction": "MNXR153595"},
    "notes": {"PROTEIN_CLASS": "2.1.2.1",
              "folinic_acid_module_curation": f"{CID}: cytosolic source of 5-formyl-THF m159, which had none, so the "
                                              "FAU1-type R64 could never carry flux. Rule as R358 (cytosolic SHMT).",
              "folinic_acid_module_evidence": "Stover and Schirch 1990: cytosolic SHMT hydrolyses methenyl-THF to "
                                              "5-formyl-THF. Holmes and Appling 2002: SHMT is the in vivo source in "
                                              "S. cerevisiae. Charge imbalance -1 as in R748: methenyl-THF m160 is "
                                              "stored at charge 0 with the formula of the -1 species (ChEBI 57455).",
              "folinic_acid_module_evidence_tier": TIER},
    "group": "One carbon pool by folate",
}

species = {mid: {k: before["metabolites"][mid][k] for k in ("name", "formula", "charge", "compartment")}
           for mid in ("m159[C_cy]", "m160[C_cy]", "m32[C_cy]", "m10[C_cy]")}
species["m159[C_cy]"]["reactions_before"] = before["metabolites"]["m159[C_cy]"]["reactions"]
species["m160[C_cy]"]["reactions_before"] = before["metabolites"]["m160[C_cy]"]["reactions"]

genes = {}
for gid, row in before["genes"].items():
    if gid == "YALI1D01489g":
        continue
    rb = row["reactions"]
    ra = [r for r in rb if r != "R748"] if gid == "YALI1F04082g" else sorted(rb + ["R748", "R_THFAT"])
    genes[gid] = {"annotation": row["annotation"], "reactions_before": rb, "reactions_after": ra}

spec = {
    "schema_version": 1,
    "curation_id": CID,
    "status": "user_authorized",
    "authorization": "User, 2026-11-12: 'I think there is a gap to fix', 'Authorize F1+F2+F4', 'Approve plan'.",
    "reviewed_on": "e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2",
    "reviewer": "Claude Code session 2026-11-12 (folinic acid review); independent audit in the report",
    "reactions": {"R748": {"before": pinned(R["R748"]), "after": r748_after},
                  "R990": {"before": pinned(R["R990"]), "after": r990_after}},
    "new_reaction": new,
    "species": species,
    "genes": genes,
    "evidence": [
        {"claim": "Rhea 34767 is methenyl-THF + H2O = 5-formyl-THF + H+; reviewed enzymes are SHMTs",
         "url": "https://www.rhea-db.org/rhea/34767", "tier": "verified 2026-11-12 (Rhea and UniProt REST)"},
        {"claim": "SHMT hydrolyses methenyl-THF to 5-formyl-THF (cytosolic and mitochondrial rabbit, E. coli)",
         "url": "https://pubmed.ncbi.nlm.nih.gov/2201683/", "tier": "verified 2026-11-12, abstract only"},
        {"claim": "SHMT is the in vivo 5-formyl-THF source in S. cerevisiae; FAU1 is the only MTHFS",
         "url": "https://pubmed.ncbi.nlm.nih.gov/11923304/", "tier": "verified 2026-11-12, abstract only"},
        {"claim": "KEGG R02300 is filed under EC 2.1.2.10 / K00605",
         "url": "https://www.kegg.jp/entry/R02300", "tier": "verified 2026-11-12 (KEGG REST)"},
        {"claim": "EC 2.1.2.10 is the lipoyl-aminomethyl to methylene-THF step only",
         "url": "https://enzyme.expasy.org/EC/2.1.2.10", "tier": "verified 2026-11-12"},
        {"claim": "YALI1D28569g (yli:2910775) and YALI1E19557g (yli:2911591) are K00600 SHMT; YALI1F04082g "
                  "(yli:2907779) is K00605", "url": "https://www.kegg.jp/entry/yli:2911591",
         "tier": "verified 2026-11-12 (KEGG REST); orthology only"},
        {"claim": "Methenyl-THF C20H20N7O6 is the -1 species", "url":
            "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:57455", "tier": "verified 2026-11-12 (ChEBI API)"},
        {"claim": "iMM904 THFATm carries the glycine-cleavage rule; FTHFI is iMM904-only and gene-less",
         "url": "https://bigg.ucsd.edu/models/iMM904/reactions/THFATm", "tier": "verified 2026-11-12 (BiGG API)"},
    ],
}
out_path.write_text(json.dumps(spec, indent=1, sort_keys=False) + "\n")
print("wrote", out_path)
