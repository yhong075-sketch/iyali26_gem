"""Write model/curation/biomass/trehalose_biomass.json from the reviewed baseline (read-only on the model).

The accepted before-state of biomass_C (every coefficient except trehalose, bounds, note keys) and the trehalose
species context are read from the baseline file, so the guard pins exactly what was reviewed.

Usage (from the repo root): python model/reports/trehalose_biomass_20261030/scripts/write_curation.py BASELINE.xml OUT.json
"""
import hashlib
import json
import sys
from pathlib import Path

from cobra.io import read_sbml_model

BASELINE_SHA = "e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897"
BASELINE_PATH = ("model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split"
                 "_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml")
BIOMASS, TERM, COEF = "biomass_C", "m294[C_cy]", -0.00686
AUTHORIZATION = ("2026-10-30 user: Approve C on e1e7c868 (reports/trehalose_biomass_20261030 TASK.md section 9: "
                 "trehalose m294[C_cy] removed from biomass_C, flag --trehalose-biomass-removal; option D, the R143 rule, "
                 "deferred to its own round). No commit or push.")
PREVIOUS = ("m294[C_cy] trehalose at 0.00686 mmol/gDW (6.86 nmol/mg DW, 2.35 mg/gDW), inherited from the source model. "
            "The same number is hard-coded without a source in the eciYali5-GEM build script (xCARBOHYDRATE); it is the "
            "iYali4 xBIOMASS trehalose share (0.0032) scaled about 2.14 times with glucan, mannan and chitin. "
            "No Yarrowia lipolytica measurement behind it was found.")
EVIDENCE = ("Flores et al. 2011, PLoS One, PMID 21931609, PMC3171402 (full text): Y. lipolytica wild type contains below "
            "1 nmol trehalose per mg dry weight in minimal glucose medium (Table 4; per Methods, YNB without amino acids "
            "and ammonium sulfate plus 0.1 percent glutamate, generally 2 percent glucose, 30 C), and below 1 when grown "
            "in glucose up to stationary phase or in glycerol. The tps1 disruptant (YlTPS1 = YALI0E14685g = "
            "YALI1E17616g, UniProt O74932, cDNA AJ011032) grows on glucose at least as fast as the wild type (doubling "
            "time 139 vs 149 min); no trehalose is detected in tps1 with nth1 also disrupted or after heat shock (the "
            "content of the single tps1 mutant at 30 C is not reported). Trehalose is therefore not an obligate biomass "
            "component under this condition. Report reports/trehalose_biomass_20261030.")
TIER = ("direct Y. lipolytica measurement (an upper bound: reported as below 1 nmol/mg DW; the paper states no "
        "detection limit) plus a Y. lipolytica mutant phenotype; strain CJM 645 on YNB-glutamate glucose, not PO1f on "
        "the Task A medium; heat shock (not modelled) raises trehalose to 7-20 nmol/mg DW")
SCOPE = ("Only the biomass_C trehalose term is removed; biomass_C is not re-normalised. The trehalose terms of the "
         "non-objective xBIOMASS, newBiom and R1372, the closed R1387 and R1710, and the R143 rule are unchanged.")

model_path, out = Path(sys.argv[1]), Path(sys.argv[2])
if hashlib.sha256(model_path.read_bytes()).hexdigest() != BASELINE_SHA:
    raise SystemExit("baseline SHA-256 differs")
model = read_sbml_model(model_path)  # Path object: cobra rejects long str paths on macOS
biomass = model.reactions.get_by_id(BIOMASS)
term = model.metabolites.get_by_id(TERM)
stoichiometry = {m.id: c for m, c in sorted(biomass.metabolites.items(), key=lambda x: x[0].id)}
if stoichiometry.pop(TERM) != COEF:
    raise SystemExit("trehalose coefficient differs")

notes = {
    "trehalose_biomass_curation_id": "TREHALOSE-BIOMASS-20261030",
    "trehalose_biomass_status": "user_authorized_biomass_term_removal",
    "trehalose_biomass_authorization": AUTHORIZATION,
    "trehalose_biomass_previous_term": PREVIOUS,
    "trehalose_biomass_evidence": EVIDENCE,
    "trehalose_biomass_evidence_tier": TIER,
    "trehalose_biomass_scope": SCOPE,
}
if any(c in v for v in notes.values() for c in "<>&"):
    raise SystemExit("note text must stay SBML-safe")

spec = {
    "schema_version": 1,
    "curation_id": "TREHALOSE-BIOMASS-20261030",
    "status": "user_authorized_biomass_term_removal",
    "authorization": AUTHORIZATION,
    "reviewed_model": {"path": BASELINE_PATH, "sha256": BASELINE_SHA},
    "task": "model/reports/trehalose_biomass_20261030/TASK.md",
    "report": "model/reports/trehalose_biomass_20261030/REPORT.md",
    "reviewer": "Claude Code session 2026-10-30 (evidence and implementation); decision by the user; independent audit in the report",
    "requires_objective": BIOMASS,
    "biomass_term": {
        "biomass_reaction": BIOMASS,
        "metabolite_id": TERM,
        "name": term.name, "formula": term.formula, "charge": term.charge, "compartment": term.compartment,
        "coefficient_before": COEF,
        "reactions_before": sorted(r.id for r in term.reactions),
        "reactions_after": sorted(r.id for r in term.reactions if r.id != BIOMASS),
        "notes": notes,
    },
    "accepted_before_state": {
        "id": "e1e7c868_biomass_C",
        "biomass_bounds": list(biomass.bounds),
        "biomass_stoichiometry_without_term": stoichiometry,
        "biomass_note_keys": sorted(biomass.notes),
    },
    "evidence_tier": TIER,
    "evidence": [
        {"url": "https://europepmc.org/article/MED/21931609",
         "full_text": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3171402/fullTextXML",
         "full_text_sha256": "c5529d322a24927ab4179ffd2164a66dd91bc4a9f6b57539795896a31f7f539f",
         "claim": "WT trehalose below 1 nmol/mg dry weight on minimal glucose medium (Table 4) and when grown in glucose "
                  "up to stationary phase or in glycerol; the tps1 disruptant grows on glucose at least as fast as the "
                  "WT; no trehalose is detected in nth1 tps1 or heat-shocked tps1",
         "tier": "direct Y. lipolytica measurement and mutant phenotype (full text read 2026-10-30)"},
        {"url": "https://rest.uniprot.org/uniprotkb/O74932",
         "claim": "YALI0E14685g = TPS1, EMBL AJ011032 (the Flores 2011 cDNA), cites PMID 21931609",
         "tier": "identity cross-check"},
        {"url": "https://github.com/bioengtaltech/eciYali5-GEM",
         "claim": "code/debiaggi2025/generate_iYali5GEM.m line 56 hard-codes xCARBOHYDRATE trehalose -0.00686 with no "
                  "source; iYali4 xBIOMASS trehalose 0.0032 (local copies in the research workspace)",
         "tier": "provenance of the removed value"},
    ],
    "out_of_scope": SCOPE,
}
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(spec, indent=1) + "\n")
print(out, len(stoichiometry), "other biomass_C entries,", len(spec["accepted_before_state"]["biomass_note_keys"]),
      "note keys; species reactions", spec["biomass_term"]["reactions_before"])
