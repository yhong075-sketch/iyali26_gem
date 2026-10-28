"""Write model/curation/thiamine/tpp_requirement.json from before.json plus the reviewed evidence text.

The before-values (stoichiometry, annotation, notes) are copied from before.json so the guard compares
exact strings; the evidence and rationale are written here by hand.
Usage (from the worktree root): python make_curation_json.py BEFORE.json OUT.json
"""
import json
import sys
from xml.sax.saxutils import unescape

before_path, out = sys.argv[1:3]
before = json.load(open(before_path))
assert before["sha256"] == "cf0825b86c3533d8a12f697a994eba3d740671a68a367143ee32cbfdc8686b7c"
biomass, r1342 = before["targets"]["biomass_C"], before["targets"]["R1342"]
COEFFICIENT = -1.2e-06
# The whole cytosolic TPP neighbourhood: its reactions besides the targets and their dead-end products.
SPECIES = ["m975[C_cy]", "m976[C_cy]", "m1116[C_mi]"]
REQUIRED = ["R752", "R753", "R933"]

spec = {
    "schema_version": 1,
    "curation_id": "TPP-REQUIREMENT-20261022",
    "status": "user_authorized_biomass_term_and_closure",
    "authorization": ("2026-10-22 user: Approved: implement the TPP plan as written: --tpp-requirement on top of "
                      "cf0825b8, TPP in biomass_C at the yeast-GEM value 1.2e-6 mmol/gDW, R1342 closed (not "
                      "deleted), ThiL BLAST check, KO comparison included. No commit or push."),
    "baseline_model_sha256": before["sha256"],
    "report": "model/reports/tpp_requirement_20261022/",
    # reactions_before: producers and consumers other than the biomass_C term this change adds.
    "species": {mid: {key: before["species"][mid][key] for key in ("name", "formula", "charge", "compartment")}
                | {"reactions_before": before["species"][mid]["reactions"]} for mid in SPECIES},
    "requires_reactions": {rid: {"gpr": before["context"][rid]["gpr"], "bounds": before["context"][rid]["bounds"],
                                 "stoichiometry": before["context"][rid]["stoichiometry"]} for rid in REQUIRED},
    "targets": {
        "biomass_C": {
            "name": biomass["name"], "gpr": biomass["gpr"], "bounds_before": biomass["bounds"],
            "bounds_after": biomass["bounds"], "annotation": biomass["annotation"],
            # Unescaped, as the notes are during a build; the guard unescapes before comparing.
            "notes_before": {key: unescape(value, {"&apos;": "'", "&quot;": '"'})
                             for key, value in biomass["notes"].items()},
            "before_stoichiometry": biomass["stoichiometry"],
            "after_stoichiometry": {**biomass["stoichiometry"], "m975[C_cy]": COEFFICIENT},
            "evidence": ("yeast-GEM 9.1.1 cofactor pseudoreaction r_4598 consumes s_1475, named TDP (thiamine "
                         "diphosphate: KEGG C00068, BiGG thmpp, C12H16N4O7P2S), at 1.2e-6 mmol/gDW; "
                         "https://github.com/SysBioChalmers/yeast-GEM/blob/main/model/yeast-GEM.yml (raw file "
                         "SHA-256 2ddfaaa73ead41ab243127d1b34fa56db2c86052f6eabcd12bb4eea93947dc17, read 2026-10-22)"),
            "evidence_tier": ("coefficient verified in yeast-GEM 9.1.1 on 2026-10-22; it is an S. cerevisiae value, "
                              "not measured in Y. lipolytica, and is provisional. That TPP is a required cofactor "
                              "(pyruvate dehydrogenase, 2-oxoglutarate dehydrogenase, transketolase, pyruvate "
                              "decarboxylase) is established/uncited. m975 has no formula or cross-reference; its "
                              "identity as thiamine diphosphate is inferred from its reactions R752 (KEGG R00619), "
                              "R753 (KEGG R00616) and R933 (thiamine diphosphate transport)."),
            "rationale": ("No reaction drained cytosolic TPP (m975) to a sink: its consumers R753 (to thiamine "
                          "triphosphate, m976) and R933 (into mitochondria, m1116) end in dead ends. So TPP "
                          "synthesis and THI80 (R752) carried no flux and a THI80 knockout had no effect. Diluting "
                          "TPP into new biomass makes its synthesis required, as yeast-GEM does for TDP."),
            "previous_state": "biomass_C had no TPP term",
        },
        "R1342": {
            "name": r1342["name"], "gpr": r1342["gpr"], "bounds_before": r1342["bounds"],
            "bounds_after": [0.0, 0.0], "annotation": r1342["annotation"], "notes_before": r1342["notes"],
            "before_stoichiometry": r1342["stoichiometry"], "after_stoichiometry": r1342["stoichiometry"],
            "evidence": ("KEGG REST link/yli/K00946 and link/sce/K00946 (thiL, thiamine-monophosphate kinase, EC "
                         "2.7.4.16) return no genes; link/yli/K00949 (thiamine pyrophosphokinase) returns "
                         "yli:2911892 = YALI2_E01340g, the old locus tag of W29 YALI1E25443g in GCF_001761485.1 "
                         "(read 2026-10-22). blastp of all 27 Swiss-Prot ThiL entries (release 2026_03) against "
                         "the W29 proteome finds no hit at E of at most 1e-5 (best E = 0.006 over 27% of the query); the "
                         "THI80 control P35202 finds YALI1E25443g at E = 1.8e-24 "
                         "(model/reports/tpp_requirement_20261022/thil_blast_w29.csv)."),
            "evidence_tier": ("verified on 2026-10-22: no KEGG K00946 assignment in yli (the CLIB122 genome) or sce, "
                              "and no ThiL BLAST hit in W29. The THI80 control is a different protein family; it "
                              "shows that the database and search work, not that a divergent ThiL would be found. "
                              "Absence of a sequence homologue does not prove that no enzyme phosphorylates TMP. "
                              "The yeast route, TMP to thiamine to TPP through THI80, is established/uncited."),
            "rationale": ("R1342 is gene-less and, once TPP is required, was the only route to TPP that bypasses "
                          "THI80 (R752). It is closed, not deleted, so its equation and cross-references stay "
                          "for review."),
            "previous_state": "bounds (0, 1000)",
        },
    },
}
# libSBML drops a note holding < or > and escapes ' and " on the next read, so new note texts avoid them.
for row in spec["targets"].values():
    for key in ("evidence", "evidence_tier", "rationale", "previous_state"):
        assert not set(row[key]) & set("<>&'\""), (key, row[key])
assert not set(spec["authorization"]) & set("<>&'\"")
with open(out, "w") as handle:
    json.dump(spec, handle, indent=1)
    handle.write("\n")
print("wrote", out, "biomass terms", len(spec["targets"]["biomass_C"]["after_stoichiometry"]))
