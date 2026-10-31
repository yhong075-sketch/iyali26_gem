"""Write model/curation/gpr/r24_homocitrate_aconitase.json from before.json and the decision (option C).

Usage (from the repo root): python3 model/reports/r24_homocitrate_aco_20261030/scripts/write_curation.py
Note text avoids apostrophes, quotes, ampersands and angle brackets: libSBML writes them as entities, which breaks the
SBML roundtrip test.
"""
import json
from pathlib import Path

REPORT = Path("model/reports/r24_homocitrate_aco_20261030")
REVIEW = Path("model/reports/r24_homocitrate_review_20261029")
OUT = Path("model/curation/gpr/r24_homocitrate_aconitase.json")
before = json.loads((REPORT / "before.json").read_text())
LYS4, ACO1, ACO2 = "YALI1E03387g", "YALI1D11984g", "YALI1E17935g"
TARGET = "R24"
AFTER = f"{ACO1} or {ACO2}"
AUTH = ("2026-10-30 user: Option C, then Approve plan (reports/r24_homocitrate_aco_20261030: R24 to YALI1D11984g or "
        "YALI1E17935g, flag --r24-homocitrate-aconitase, built after --aconitase-aco2-split). No push.")
TIER = (f"{ACO1}: homology plus S. cerevisiae and fungal biochemistry and genetics. {ACO2}: homology plus "
        "S. cerevisiae Aco2p biochemistry only; inactive Aco2-clade homologues exist in Aspergillus. No Y. lipolytica "
        "enzymology, knockout or localisation for either aconitase on homocitrate.")
EVIDENCE = (
    f"{LYS4} (YALI0E02728g, Swiss-Prot Q6C791 LYS4, by similarity) is a homoaconitase; its only reviewed activity is "
    "homoisocitrate to cis-homoaconitate (RHEA:15485), the R395 step, and it stays on R395. Yeast and A. fumigatus "
    "homoaconitases show no homocitrate activity; homocitrate is converted by aconitase plus homoaconitase (Fazius "
    "2012, PMC3556520, Table 1 and Results). In S. cerevisiae each aco1 or aco2 single mutant grows without lysine, "
    "while the aco1 aco2 double mutant is completely lysine auxotrophic (same paper), so either aconitase alone "
    f"supports the step; R24 therefore reads {ACO1} (ACO1-type, Q6C9P6) or {ACO2} (ACO2-type, Q6C5V2). The ACO2-type "
    "member stands on transfer from S. cerevisiae: Aspergillus Aco2-clade AcoB proteins were inactive, and the "
    "Aco2-clade lysine at Aco2p 610 is shared by those inactive proteins. Report reports/r24_homocitrate_aco_20261030; "
    "review reports/r24_homocitrate_review_20261029.")

r = before["reactions"][TARGET]
notes = {
    "r24_homocitrate_aco_curation_id": "R24-HOMOCITRATE-ACO-20261030",
    "r24_homocitrate_aco_status": "user_authorized_gpr_correction",
    "r24_homocitrate_aco_authorization": AUTH,
    "r24_homocitrate_aco_previous_gpr": f"{LYS4}; the LYS4-like homoaconitase, kept on R395",
    "r24_homocitrate_aco_evidence": EVIDENCE,
    "r24_homocitrate_aco_evidence_tier": TIER,
}
reactions = {TARGET: {
    "reaction_name": r["name"], "before_gpr": r["gpr"], "after_gpr": AFTER, "bounds": r["bounds"],
    "stoichiometry": r["stoichiometry"],
    "species": {m: {k: v[k] for k in ("formula", "charge", "compartment")} for m, v in r["species"].items()},
    "annotation": r["annotation"], "before_notes": r["notes"], "notes": notes,
}}
after = {LYS4: ["R395"], ACO1: sorted([*before["genes"][ACO1]["reactions"], TARGET]), ACO2: [TARGET]}
spec = {
    "schema_version": 1,
    "assignment_id": "R24-HOMOCITRATE-ACO-20261030",
    "status": "user_authorized_gpr_correction",
    "authorization": AUTH,
    "reactions": reactions,
    "genes": {g: {"annotation": before["genes"][g]["annotation"], "reactions_before": before["genes"][g]["reactions"],
                  "reactions_after": after[g]} for g in (LYS4, ACO1, ACO2)},
    "reviewed_model": {"path": "model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_"
                               "thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml",
                       "sha256": before["baseline_model"]["sha256"]},
    "not_in_scope": "R24 name, cross-references and PROTEIN_CLASS note; m76 cross-references; the 3 H+ of R24 and "
                    "R395 (review options N and H)",
    "task": "model/reports/r24_homocitrate_aco_20261030/TASK.md",
    "report": "model/reports/r24_homocitrate_aco_20261030/REPORT.md",
    "reviewer": "Claude Code session 2026-10-29/30 (evidence and implementation); decision by the user; independent "
                "audit in the report",
    "evidence_tier": TIER,
    "evidence": [
        {"url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3556520/",
         "claim": "Fazius 2012 full text: Lys4p and LysF show no activity with homocitrate (Table 1); aconitase plus "
                  "homoaconitase convert homocitrate to homoisocitrate (yeast part not shown); aco1 and aco2 single "
                  "mutants grow without lysine, the double mutant is completely lysine auxotrophic; Aspergillus AcoB "
                  "proteins inactive",
         "tier": "S. cerevisiae and fungal biochemistry and genetics (full text)",
         "raw": f"{REVIEW}/evidence/PMC3556520_fulltext.txt"},
        {"url": "https://rest.uniprot.org/uniprotkb/Q6C791.json",
         "claim": f"{LYS4}: Swiss-Prot Homoaconitase, mitochondrial, LYS4, EC 4.2.1.36, catalytic activity "
                  "RHEA:15485 only, by similarity", "tier": "homology (reviewed entry)",
         "raw": f"{REVIEW}/evidence/uniprot_Q6C791.json"},
        {"url": "https://rest.uniprot.org/uniprotkb/Q6C9P6.json",
         "claim": f"{ACO1}: ACO1-type aconitase, PANTHER PTHR43160:SF3", "tier": "homology",
         "raw": f"{REVIEW}/evidence/uniprot_Q6C9P6.json"},
        {"url": "https://rest.uniprot.org/uniprotkb/Q6C5V2.json",
         "claim": f"{ACO2}: ACO2-type aconitase, PANTHER PTHR43160:SF2 homocitrate dehydratase", "tier": "homology",
         "raw": f"{REVIEW}/evidence/uniprot_Q6C5V2.json"},
        {"url": "https://rest.uniprot.org/uniprotkb/P39533.fasta",
         "claim": f"residue aligned to Aco2p K610: R in Aco1p, {ACO1} and Af AcoA; K in {ACO2} and in the inactive "
                  "Af and An AcoB (not diagnostic of activity)", "tier": "homology",
         "raw": f"{REVIEW}/logs/aco2_k610_position.log"},
    ],
}
assert not any(c in "".join(notes.values()) for c in "'\"&<>"), "note text must not need XML entities"
OUT.write_text(json.dumps(spec, indent=1) + "\n")
print(OUT, sorted(reactions), {g: v["reactions_after"] for g, v in spec["genes"].items()})
