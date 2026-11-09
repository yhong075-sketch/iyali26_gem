"""Write model/curation/gpr/vatpase_a1_gpr.json from before.json and the decision (A1 + E-1).

Usage (from the repo root): python3 model/reports/vatpase_gpr_a1_20261103/scripts/write_curation.py
"""
import json
from pathlib import Path

REPORT = Path("model/reports/vatpase_gpr_a1_20261103")
OUT = Path("model/curation/gpr/vatpase_a1_gpr.json")
before = json.loads((REPORT / "before.json").read_text())
CORE = ["YALI1A11258g", "YALI1E14125g", "YALI1D00581g", "YALI1E10492g", "YALI1A09766g", "YALI1F20965g",
        "YALI0E16192g", "YALI1F31854g", "YALI1F21690g", "YALI1E37063g", "YALI1E32332g", "YALI1F13017g"]
NEW_GENE = "YALI1F25663g"
A_SUBUNITS = ["YALI1E12482g", "YALI1F38820g"]
AFTER = " and ".join(CORE + [NEW_GENE]) + f" and ({' or '.join(A_SUBUNITS)})"
SUBUNITS = {
    "YALI1A09766g": "V1 A", "YALI1E32332g": "V1 B", "YALI1A11258g": "V1 C", "YALI1D00581g": "V1 D",
    "YALI1F20965g": "V1 E", "YALI0E16192g": "V1 F", "YALI1E14125g": "V1 G", "YALI1E10492g": "V1 H",
    "YALI1F31854g": "V0 c", "YALI1E37063g": "V0 c′", "YALI1F13017g": "V0 c″", "YALI1F21690g": "V0 d",
    "YALI1E12482g": "V0 a (a1)", "YALI1F38820g": "V0 a (a2, Vph1-like)", NEW_GENE: "V0 e",
}
# Note text avoids ' " & < >: the SBML writer escapes them, so a written model would not round-trip.
AUTH = ("2026-11-03 user: Yes, with E-1 + Scope-1, answering Run the V-ATPase round with this plan? after accepting "
        "the recommended line (Std-§5, Fgate-lift, A1, Label-1, Scope-1, ID-keep, Verify-1) and the E-2 guide check: "
        "R794 and R795 get the 12 V1/V0 core subunits and V0 e YALI1F25663g AND (a1 YALI1E12482g or a2 YALI1F38820g); "
        "flag --vatpase-a1-gpr on baseline 2bc263f2. No commit or push.")
EVIDENCE = (
    "13 of the 15 rule genes are reciprocal best hits of S. cerevisiae V-ATPase subunits (BLAST 2026-11-01; H weakest, "
    "26.8%). The two exceptions: a1 YALI1E12482g, whose best S. cerevisiae hit is VPH1 (876 vs 811 bits over STV1), so "
    "its Stv1-type assignment is not reciprocal; and V1 F, which has no W29 protein (CLIB122 Q6C5Q2 best hit VMA7, "
    "forward only). The former second branch (a2 and c\u2033) had no V1 head and no c/c\u2032 ring; c\u2033 null "
    "loses all activity (Hirata 1997, PMID 9030535), so c\u2033 joins the core. Both S. cerevisiae a isoforms form "
    "active pumps (Manolson 1994, PMID 7514599; Kawasaki-Nishi 2001, PMID 11592965; Perzov 2002, PMID 11948198), "
    "which meets the AGENTS.md section 5 OR test when carried to Y. lipolytica by homology; same-compartment "
    "substitution is partial and one-way (Deschamps 2020, PMC7002768), so the 2026-09-21 same-compartment test is "
    "deliberately not met until Y. lipolytica localisation data exist. V0 e YALI1F25663g: reciprocal best hit of Vma9 "
    "(61.4%); the e null fails V1-V0 assembly (Compton 2006, PMID 16569636); against: once assembled, e is not "
    "needed for coupled pumping in vitro (Bueler 2015, PMID 25546637); its Cas12a guides did not cut (CS -1.60, "
    "-0.36) while four Cas9 guides cut and deplete. V1 F YALI0E16192g: the W29 pseudogene call is taken to be a "
    "missed 2-nt first exon plus a 56-nt intron, giving an intact ORF when spliced (inferred; no RNA evidence "
    "checked) (Fgate-lift).")
TIER = ("orthology (reciprocal BLAST) and S. cerevisiae mutant and biochemical literature; no Y. lipolytica "
        "localisation, knockout or complex data; homology-based candidate rule, not a verified native rule")
notes = {
    "vatpase_a1_gpr_curation_id": "VATPASE-A1-GPR-20261103",
    "vatpase_a1_gpr_status": "user_authorized_gpr_correction",
    "vatpase_a1_gpr_authorization": AUTH,
    "vatpase_a1_gpr_previous_gpr": before["reactions"]["R794"]["gpr"],
    "vatpase_a1_gpr_standard": ("meets AGENTS.md section 5 by homology (each S. cerevisiae a isoform forms an active "
                                "pump on its own; isoform assignment of a1 and a2 has moderate support); does not "
                                "meet the 2026-09-21 same-compartment substitution recommendation"),
    "vatpase_a1_gpr_evidence": EVIDENCE,
    "vatpase_a1_gpr_evidence_tier": TIER,
}
reactions = {}
for rid in ("R794", "R795"):
    r = before["reactions"][rid]
    assert r["gpr"] == before["reactions"]["R794"]["gpr"]
    reactions[rid] = {"reaction_name": r["name"], "bounds": r["bounds"], "stoichiometry": r["stoichiometry"],
                      "species": r["species"], "annotation": r["annotation"], "before_gpr": r["gpr"],
                      "after_gpr": AFTER, "before_notes": r["notes"], "notes": notes}
genes = {g: {"subunit": SUBUNITS[g], "annotation": v["annotation"], "reactions_before": v["reactions"],
             "reactions_after": v["reactions"]} for g, v in before["genes"].items()}
assert set(genes) == set(CORE + A_SUBUNITS) and not before["new_gene_present"]
spec = {
    "schema_version": 1,
    "assignment_id": "VATPASE-A1-GPR-20261103",
    "status": "user_authorized_gpr_correction",
    "authorization": AUTH,
    "baseline_sha256": before["baseline_model"]["sha256"],
    "reactions": reactions,
    "genes": genes,
    "new_gene": {
        "id": NEW_GENE,
        "subunit": SUBUNITS[NEW_GENE],
        "annotation": {"sbo": "SBO:0000243", "ncbigene": "90949792", "refseq": "XP_065950526.1"},
        "reactions_after": ["R794", "R795"],
        "identity": ("NCBI Gene 90949792, YALI1_F25663g (old tag YALI2_F00374g), NC_090775.1:2566174-2566753 minus, "
                     "XP_065950526.1, 75 aa, RefSeq 'partial'; no CLIB122 ID, UniProt or KEGG entry verified; no "
                     "verified name. Source: reference/ncbi/clib89_feature_table.txt"),
        "orthology": "reciprocal best hit of S. cerevisiae VMA9 Q3E7B6: 61.4% identity over 70 aa, 5e-26 (2026-11-01)",
        "screen": ("Cas9 FS -5.23 essential (6 guides; 4 cut with CS 5.3-7.2); Cas12a FS +0.45 non-essential from 2 "
                   "guides that did not cut (CS -1.60, -0.36); consensus essential"),
    },
    "notes": notes,
}
OUT.write_text(json.dumps(spec, indent=1, ensure_ascii=False) + "\n")
print(OUT, AFTER)
