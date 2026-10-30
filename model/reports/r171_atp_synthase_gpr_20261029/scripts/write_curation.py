"""Write model/curation/gpr/r171_atp_synthase_gpr.json from before.json and the decision (Option D).

Usage (from the repo root): python3 model/reports/r171_atp_synthase_gpr_20261029/scripts/write_curation.py
"""
import json
from pathlib import Path

REPORT = Path("model/reports/r171_atp_synthase_gpr_20261029")
OUT = Path("model/curation/gpr/r171_atp_synthase_gpr.json")
before = json.loads((REPORT / "before.json").read_text())
CORE = ["YALI1M00260g", "YALI1M00267g", "YALI1M00318g", "YALI1F04433g", "YALI1F04110g", "YALI1D27826g",
        "YALI1D15734g", "YALI1D14639g", "YALI1B08885g", "YALI1F07384g", "YALI1B05364g", "YALI1F26967g"]
UNLINKED = {
    "YALI1B28025g": {"subunit": "g, ATP20-like", "identity": "YALI0B21527g, UniProt B5FVB8, 138 aa",
                     "previous_role": "single-gene OR branch",
                     "reason": "a single accessory subunit has no evidence of catalysing ATP synthesis alone; "
                               "dimer-specific (Hahn 2016), not found in the purified Y. lipolytica enzyme (Liu 2015)"},
    "YALI1B15733g": {"subunit": "k, ATP19-like", "identity": "YALI0B11913g, UniProt B5FVB3, 70 aa",
                     "previous_role": "single-gene OR branch",
                     "reason": "as g; not found in the purified Y. lipolytica enzyme (Liu 2015)"},
    "YALI1E38131g": {"subunit": "e, TIM11/ATP21-like", "identity": "YALI0E32164g, UniProt B5FVG3, 90 aa",
                     "previous_role": "AND member",
                     "reason": "dimer-specific; the Y. lipolytica monomer lacking e, g, k hydrolyses ATP at about "
                               "the rate of the dimer (Hahn 2016); S. cerevisiae e, g, k not required for an "
                               "enzymatically active ATP synthase (Arnold 1998)"},
    "YALI0D17490g": {"subunit": "j, ATP18-like", "identity": "CLIB122 id YALI0D17490g (S2 map: YALI1_D21476g), "
                                                             "UniProt Q6C8S0, 62 aa",
                     "previous_role": "AND member",
                     "reason": "S. cerevisiae null viable with about half the specific ATPase activity and a low "
                               "ADP/O ratio (Vaillier 1999, abstract): loss impairs, does not abolish, the activity; "
                               "one S. cerevisiae abstract, no Y. lipolytica data, not screened"},
}
AUTH = ("2026-10-29 user: Option D (Recommended) for R171, answering the question which R171 option to build, "
        "with the plan of TASK.md, on baseline d8d1bac7: 12-subunit AND, single-gene OR branches g and k removed, subunits e and j "
        "removed from the rule, the four genes kept in the model unlinked; flag --r171-atp-synthase-gpr. No commit or "
        "push.")
after_gpr = " and ".join(CORE)
r = before["reaction"]
notes = {
    "r171_atp_synthase_gpr_curation_id": "R171-ATP-SYNTHASE-GPR-20261029",
    "r171_atp_synthase_gpr_status": "user_authorized_gpr_correction",
    "r171_atp_synthase_gpr_authorization": AUTH,
    "r171_atp_synthase_gpr_previous_gpr": r["gpr"],
    "r171_atp_synthase_gpr_unlinked_genes": "; ".join(f"{g} (subunit {v['subunit']})" for g, v in UNLINKED.items())
                                            + "; kept in the model with no reaction",
    "r171_atp_synthase_gpr_evidence": (
        "AND members are the subunits whose loss is taken to close the reaction: ATP6/8/9 (mtDNA) and the nine "
        "nuclear core subunits (alpha, beta, gamma, delta, OSCP, b, d, h, f). Subunits e, g, k are dimer-specific "
        "and the Y. lipolytica monomer lacking them hydrolyses ATP at about the rate of the dimer (ATP hydrolysis, not "
        "synthesis, was measured; Hahn 2016, PMID 27373333); g and k were "
        "not found in the purified Y. lipolytica enzyme (Liu 2015, PMID 25759169); S. cerevisiae e, g, k are not "
        "required for an enzymatically active ATP synthase (Arnold 1998, PMID 9857174); the S. cerevisiae subunit-j "
        "null keeps about half the ATPase activity (Vaillier 1999, PMID 9867878). OR-rule review 2026-10-29, REPORT "
        "section 7.1."),
    "r171_atp_synthase_gpr_evidence_tier": (
        "literature on the Y. lipolytica enzyme (structure, purified-enzyme composition, hydrolysis rather than "
        "synthesis) and S. cerevisiae mutants; no Y. lipolytica knockout of e, g, k or j; keep-or-drop of j is a "
        "judgement"),
}
spec = {
    "schema_version": 1,
    "assignment_id": "R171-ATP-SYNTHASE-GPR-20261029",
    "status": "user_authorized_gpr_correction",
    "authorization": AUTH,
    "reaction": {"id": "R171", "reaction_name": r["name"], "before_gpr": r["gpr"], "after_gpr": after_gpr,
                 "bounds": r["bounds"], "stoichiometry": r["stoichiometry"], "species": r["species"],
                 "annotation": r["annotation"], "before_notes": r["notes"], "notes": notes},
    "genes": {g: {"annotation": v["annotation"], "reactions_before": v["reactions"],
                  "reactions_after": [] if g in UNLINKED else v["reactions"]}
              for g, v in before["genes"].items()},
    "unlinked_genes": UNLINKED,
    "reviewed_model": before["baseline_model"],
    "review": "model/reports/or_rule_review_r171_aars_20261029/REPORT.md",
    "report": "model/reports/r171_atp_synthase_gpr_20261029/REPORT.md",
    "reviewer": "Claude Code session 2026-10-29 (review and implementation); option chosen by the user; independent "
                "audit in the report",
    "evidence_tier": notes["r171_atp_synthase_gpr_evidence_tier"],
    "evidence": [
        {"url": "https://europepmc.org/article/MED/27373333",
         "claim": "Y. lipolytica ATP synthase dimer (cryo-EM): e, g, k dimer-specific; DDM monomer lacking e, g, k "
                  "hydrolyses ATP at about the rate of the dimer", "tier": "Y. lipolytica structure and enzymology"},
        {"url": "https://europepmc.org/article/MED/25759169",
         "claim": "purified Y. lipolytica ATP synthase: subunits g, k and l not found; e listed (absent from the DDM "
                  "preparation)", "tier": "Y. lipolytica purified-enzyme composition"},
        {"url": "https://europepmc.org/article/MED/9857174",
         "claim": "S. cerevisiae: e, g, k not required for enzymatically active ATP synthase; e and g needed for "
                  "the dimer", "tier": "S. cerevisiae mutants (abstract)"},
        {"url": "https://europepmc.org/article/MED/9867878",
         "claim": "S. cerevisiae atp18 null: viable, low growth on non-fermentable medium, about 2-fold less "
                  "specific ATPase activity, low ADP/O", "tier": "S. cerevisiae mutant (abstract)"},
        {"url": "https://www.uniprot.org/uniprotkb/B5FVB8/entry", "claim": "YALI0B21527g: ATP synthase subunit g",
         "tier": "database"},
        {"url": "https://www.uniprot.org/uniprotkb/B5FVB3/entry", "claim": "YALI0B11913g: ATP synthase subunit K",
         "tier": "database"},
        {"url": "https://www.uniprot.org/uniprotkb/B5FVG3/entry", "claim": "YALI0E32164g: ATP synthase subunit e",
         "tier": "database"},
        {"url": "https://www.uniprot.org/uniprotkb/Q6C8S0/entry", "claim": "YALI0D17490g: ATP synthase subunit J",
         "tier": "database"},
    ],
}
OUT.write_text(json.dumps(spec, indent=1) + "\n")
print(OUT, after_gpr)
