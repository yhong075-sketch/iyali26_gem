"""Write model/curation/gpr/r673_prolyl_split.json from before.json and the decision (note wording revised 2026-10-29,
reports/r171_r673_note_revision_20261029).

Usage (from the repo root): python3 model/reports/r673_prolyl_split_20261029/scripts/write_curation.py
Note text avoids apostrophes and quotes: libSBML writes them as entities, which breaks the SBML roundtrip test.
"""
import json
from pathlib import Path

REPORT = Path("model/reports/r673_prolyl_split_20261029")
PRED = "model/reports/r673_prolyl_localization_20261029"
OUT = Path("model/curation/gpr/r673_prolyl_split.json")
before = json.loads((REPORT / "before.json").read_text())
CYTO, MITO = "YALI1E05929g", "YALI1D19110g"
AUTH = ("2026-10-29 user: Approve as written (TASK.md section 4 of reports/r673_prolyl_localization_20261029: R673 to "
        "YALI1E05929g, YALI1D19110g kept unlinked, flag --r673-prolyl-split on baseline ca43ae0c). No commit or push.")
TIER = ("localisation predictions (TargetP 2.0, DeepLoc 2.1) on the W29 sequences plus genome-annotation homology "
        "names (mixed; the ancestral-locus tags come from the same annotation); no "
        "Y. lipolytica experimental localisation; the W29 and CLIB122 gene models of YALI1E05929g differ in start site")
r = before["reaction"]
notes = {
    "r673_prolyl_split_curation_id": "R673-PROLYL-SPLIT-20261029",
    "r673_prolyl_split_status": "user_authorized_gpr_correction",
    "r673_prolyl_split_authorization": AUTH,
    "r673_prolyl_split_previous_gpr": f"{r['gpr']}; {MITO} is the AIM10-like paralog, predicted mitochondrial, and is now unlinked",
    "r673_prolyl_split_evidence": (
        f"{CYTO} (W29 XP_503566.2, 528 aa; annotation similar to S. cerevisiae YHR020W, a putative prolyl-tRNA "
        "synthetase with no UniProt subcellular-location annotation, ancestral locus Anc_1.352): TargetP 2.0 Other 0.99995, DeepLoc 2.1 "
        "Cytoplasm 0.752 with Mitochondrion 0.222 below its threshold 0.622. "
        f"{MITO} (XP_502869.1, 445 aa; annotation similar to C. albicans P78600, a prolyl-tRNA synthetase annotated "
        "cytoplasmic, and to S. cerevisiae YER087W AIM10, a probable mitochondrial prolyl-tRNA synthetase, ancestral "
        "locus Anc_7.368): TargetP mTP 0.625 with cleavage at 13-14, DeepLoc Mitochondrion 0.844. Homology names are "
        "mixed; the assignment rests on the two predictors; the ancestral-locus tags come from the same homology "
        "annotation. Prediction report "
        "reports/r673_prolyl_localization_20261029; OR-rule review 2026-10-29 section 7.2."),
    "r673_prolyl_split_evidence_tier": TIER,
}
spec = {
    "schema_version": 1,
    "assignment_id": "R673-PROLYL-SPLIT-20261029",
    "status": "user_authorized_gpr_correction",
    "authorization": AUTH,
    "reaction": {"id": "R673", "reaction_name": r["name"], "before_gpr": r["gpr"], "after_gpr": CYTO,
                 "bounds": r["bounds"], "stoichiometry": r["stoichiometry"], "species": r["species"],
                 "annotation": r["annotation"], "before_notes": r["notes"], "notes": notes},
    "genes": {g: {"annotation": v["annotation"], "reactions_before": v["reactions"],
                  "reactions_after": [] if g == MITO else v["reactions"]} for g, v in before["genes"].items()},
    "unlinked_gene": {"id": MITO, "after": "kept in the model with no reaction (as YALI1B07957g after the aaRS split)",
                      "identity": "YALI0D15620g, W29 XP_502869.1, UniProt A0A1D8NEQ1 / Q6C8Z3, 445 aa, AIM10-like "
                                  "prolyl-tRNA synthetase; no mitochondrial prolyl-tRNA reaction in the model"},
    "reviewed_model": {"path": "model/candidates/E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_"
                               "gpr_lipoylation_acylacp_cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_"
                               "nuclear_diffusion_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_20261029.xml",
                       "sha256": before["baseline_model"]["sha256"]},
    "review": "model/reports/or_rule_review_r171_aars_20261029/REPORT.md",
    "predictions": f"{PRED}/REPORT.md",
    "report": "model/reports/r673_prolyl_split_20261029/REPORT.md",
    "reviewer": "Claude Code session 2026-10-29 (predictions and implementation); decision by the user; independent "
                "audit in the report",
    "evidence_tier": TIER,
    "evidence": [
        {"url": "https://services.healthtech.dtu.dk/services/TargetP-2.0/",
         "claim": f"job 6AE3BC76002509BB05EB4B1C: {CYTO} Other 0.999947; {MITO} mTP 0.625112, cleavage 13-14",
         "tier": "prediction", "raw": f"{PRED}/raw/targetp/"},
        {"url": "https://services.healthtech.dtu.dk/services/DeepLoc-2.1/",
         "claim": f"job 6AE3BCA7002509D23808033B: {CYTO} Cytoplasm 0.7524, Mitochondrion 0.2215 (threshold 0.6220); "
                  f"{MITO} Mitochondrion 0.8443", "tier": "prediction", "raw": f"{PRED}/raw/deeploc/"},
        {"url": "https://www.uniprot.org/uniprotkb/P38708/entry",
         "claim": "S. cerevisiae YHR020W: Putative proline--tRNA ligase YHR020W, no subcellular location annotated "
                  "(best S. cerevisiae hit of YALI1E05929g)", "tier": "ortholog reference"},
        {"url": "https://www.uniprot.org/uniprotkb/P39965/entry",
         "claim": "S. cerevisiae AIM10 / YER087W: Probable proline--tRNA ligase, mitochondrial; Mitochondrion "
                  "(ECO:0000269, PMID 16823961); synteny partner of YALI1D19110g", "tier": "ortholog reference"},
        {"url": "https://www.uniprot.org/uniprotkb/P78600/entry",
         "claim": "C. albicans Proline--tRNA ligase, cytoplasmic; Cytoplasm (ECO:0000305); named in the annotation of "
                  "YALI1D19110g", "tier": "ortholog reference"},
    ],
}
OUT.write_text(json.dumps(spec, indent=1) + "\n")
text = OUT.read_text()
assert not any(c in "".join(notes.values()) for c in "'\"&<>"), "note text must not need XML entities"
print(OUT, CYTO)
