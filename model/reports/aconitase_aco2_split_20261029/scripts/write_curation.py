"""Write model/curation/gpr/aconitase_aco2_split.json from before.json and the decision.

Usage (from the repo root): python3 model/reports/aconitase_aco2_split_20261029/scripts/write_curation.py
Note text avoids apostrophes, quotes, ampersands and angle brackets: libSBML writes them as entities, which breaks the
SBML roundtrip test.
"""
import json
from pathlib import Path

REPORT = Path("model/reports/aconitase_aco2_split_20261029")
OUT = Path("model/curation/gpr/aconitase_aco2_split.json")
before = json.loads((REPORT / "before.json").read_text())
ACO1, ACO2 = "YALI1D11984g", "YALI1E17935g"
TARGETS = ("R1556", "R216", "R238")
AUTH = ("2026-10-29 user: Approve the aconitase plan in TASK.md section 7 on baseline 059316ef "
        "(reports/aconitase_aco2_split_20261029: R1556, R216, R238 to YALI1D11984g, YALI1E17935g kept unlinked, flag "
        "--aconitase-aco2-split). No commit or push.")
TIER = ("homology (PANTHER subfamily, pairwise identity to the S. cerevisiae paralogs) plus S. cerevisiae and "
        "filamentous-fungal biochemistry; no Y. lipolytica enzymology, localisation or knockout for YALI1E17935g")
EVIDENCE = (
    f"{ACO1} (YALI0D09361g, W29 XP_502616.2, UniProt Q6C9P6) is ACO1-type: PANTHER PTHR43160:SF3 aconitate "
    "hydratase, mitochondrial; 77.7 percent identity to S. cerevisiae Aco1p P19414 and 52.3 percent to Aco2p P39533. "
    f"{ACO2} (YALI0E14949g, XP_503960.1, Q6C5V2) is ACO2-type: PANTHER PTHR43160:SF2 homocitrate dehydratase, "
    "mitochondrial; 70.3 percent identity to Aco2p and 54.9 percent to Aco1p (global alignment, BLOSUM62). "
    "S. cerevisiae Aco1p serves the citric acid cycle while Aco2p contributes only to lysine biosynthesis, and the "
    "Aco2 homologues of filamentous fungi were enzymatically inactive (Fazius 2012, PMID 23106124, abstract). "
    "Multicopy Y. lipolytica ACO1 shifts citric towards isocitric acid (Holz 2009, PMID 18850095, abstract; the "
    f"locus is not named there). {ACO2} is therefore taken as homocitrate-specific or inactive, not as a citrate "
    "aconitase; it is not assigned to R24 (separate review). UniProt automatic names and GO terms call both "
    "proteins aconitate hydratase. Report reports/aconitase_aco2_split_20261029.")

notes = {
    "aconitase_aco2_split_curation_id": "ACONITASE-ACO2-SPLIT-20261029",
    "aconitase_aco2_split_status": "user_authorized_gpr_correction",
    "aconitase_aco2_split_authorization": AUTH,
    "aconitase_aco2_split_previous_gpr": f"{ACO1} or {ACO2}; {ACO2} is the ACO2-type paralog and is now unlinked",
    "aconitase_aco2_split_evidence": EVIDENCE,
    "aconitase_aco2_split_evidence_tier": TIER,
}
reactions = {}
for rid in TARGETS:
    r = before["reactions"][rid]
    reactions[rid] = {
        "reaction_name": r["name"], "before_gpr": r["gpr"], "after_gpr": ACO1, "bounds": r["bounds"],
        "stoichiometry": r["stoichiometry"],
        "species": {m: {k: v[k] for k in ("formula", "charge", "compartment")} for m, v in r["species"].items()},
        "annotation": r["annotation"], "before_notes": r["notes"], "notes": notes,
    }
spec = {
    "schema_version": 1,
    "assignment_id": "ACONITASE-ACO2-SPLIT-20261029",
    "status": "user_authorized_gpr_correction",
    "authorization": AUTH,
    "reactions": reactions,
    "genes": {g: {"annotation": before["genes"][g]["annotation"], "reactions_before": before["genes"][g]["reactions"],
                  "reactions_after": [] if g == ACO2 else before["genes"][g]["reactions"]} for g in (ACO1, ACO2)},
    "unlinked_gene": {"id": ACO2, "after": "kept in the model with no reaction (as YALI1B07957g after the aaRS split "
                                          "and YALI1D19110g after the R673 split)",
                      "identity": "YALI0E14949g, W29 XP_503960.1, UniProt Q6C5V2 (CLIB122) / A0A1H6Q7Y5 (W29), 790 aa, "
                                  "ACO2-type aconitase-family protein; homocitrate-specific or inactive by homology"},
    "reviewed_model": {"path": "model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_"
                               "thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_notes_v3_20261029.xml",
                       "sha256": before["baseline_model"]["sha256"]},
    "task": "model/reports/aconitase_aco2_split_20261029/TASK.md",
    "report": "model/reports/aconitase_aco2_split_20261029/REPORT.md",
    "reviewer": "Claude Code session 2026-10-29 (evidence and implementation); decision by the user; independent "
                "audit in the report",
    "evidence_tier": TIER,
    "evidence": [
        {"url": "https://pubmed.ncbi.nlm.nih.gov/23106124/",
         "claim": "Fazius 2012 abstract: yeast Aco1p is essential for the citric acid cycle; Aco2p specifically and "
                  "exclusively contributes to lysine biosynthesis; Aco2p homologues of filamentous fungi were "
                  "enzymatically inactive and did not complement yeast aconitase mutants",
         "tier": "S. cerevisiae and fungal biochemistry (abstract)", "raw": f"{REPORT}/evidence/pmid23106124.txt"},
        {"url": "https://pubmed.ncbi.nlm.nih.gov/18850095/",
         "claim": "Holz 2009 abstract: multicopy Y. lipolytica ACO1 shifts the citric/isocitric acid ratio towards "
                  "isocitrate; locus not named (YALI0D09361g inferred)",
         "tier": "Y. lipolytica phenotype (abstract)", "raw": f"{REPORT}/evidence/yl_aconitase_abstracts.txt"},
        {"url": "https://rest.uniprot.org/uniprotkb/Q6C9P6.json",
         "claim": f"{ACO1}: RefSeq XP_502616.2, ORF YALI0_D09361g, PANTHER PTHR43160:SF3, automatic name Aconitate "
                  "hydratase, mitochondrial (entry version 113)", "tier": "homology", "raw": f"{REPORT}/evidence/Q6C9P6.json"},
        {"url": "https://rest.uniprot.org/uniprotkb/Q6C5V2.json",
         "claim": f"{ACO2}: RefSeq XP_503960.1, ORF YALI0_E14949g, PANTHER PTHR43160:SF2 homocitrate dehydratase, "
                  "automatic name Aconitate hydratase, mitochondrial, GO aconitate hydratase activity IBA (entry "
                  "version 111)", "tier": "homology", "raw": f"{REPORT}/evidence/Q6C5V2.json"},
        {"url": "https://rest.uniprot.org/uniprotkb/P39533.fasta",
         "claim": "global alignment (Biopython 1.85, BLOSUM62, open -10, extend -0.5, end gaps free): identity over "
                  f"the aligned span {ACO1} vs Aco1p 77.7, vs Aco2p 52.3; {ACO2} vs Aco2p 70.3, vs Aco1p 54.9; "
                  "Aco1p vs Aco2p 54.7", "tier": "homology", "raw": f"{REPORT}/evidence/pairwise_alignment.tsv"},
    ],
}
assert not any(c in "".join(notes.values()) for c in "'\"&<>"), "note text must not need XML entities"
OUT.write_text(json.dumps(spec, indent=1) + "\n")
print(OUT, sorted(reactions))
