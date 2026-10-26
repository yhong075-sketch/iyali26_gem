"""Write model/curation/lipoylation/lipoylation_module.json from the approved design (run from the repo root)."""
import json
import cobra

M = cobra.io.read_sbml_model("model/candidates/E5_coq9_resp_pkg_r132_combined_20261020.xml")
g = M.metabolites.get_by_id
A = 1e-5
SRC = {"schonauer2009": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2749097/",
       "hermes2013": "https://pubmed.ncbi.nlm.nih.gov/23960015/",
       "rhea17665": "https://www.rhea-db.org/rhea/17665",
       "rhea16585": "https://www.rhea-db.org/rhea/16585",
       "rhea16413": "https://www.rhea-db.org/rhea/16413"}
new_mets = {
    "m2062[C_mi]": {"name": "apo-glycine cleavage H-protein (L-lysyl-[GcvH])_HR", "formula": "HR", "charge": 0, "compartment": "C_mi"},
    "m2063[C_mi]": {"name": "octanoyl-H-protein (N6-octanoyl-L-lysyl-[GcvH])_C8H15OR", "formula": "C8H15OR", "charge": 0, "compartment": "C_mi"},
    "m2064[C_mi]": {"name": "apo-E2 lipoyl domain (L-lysyl-[Lat1/Kgd2])_HR", "formula": "HR", "charge": 0, "compartment": "C_mi"},
    "m2065[C_mi]": {"name": "lipoyl-E2 (N6-lipoyl-L-lysyl-[Lat1/Kgd2])_C8H13OS2R", "formula": "C8H13OS2R", "charge": 0, "compartment": "C_mi"},
    "m2066[C_mi]": {"name": "L-cysteine_C3H7NO2S", "formula": "C3H7NO2S", "charge": 0, "compartment": "C_mi"},
    "m2067[C_mi]": {"name": "5'-deoxyadenosine_C10H13N5O3", "formula": "C10H13N5O3", "charge": 0, "compartment": "C_mi"},
}
new_rxns = {
    "R2286": {"name": "octanoyl-[ACP]:H-protein N-octanoyltransferase, mitochondrial (LIP2-type)",
              "stoichiometry": {"m1421[C_mi]": -1.0, "m2062[C_mi]": -1.0, "m2063[C_mi]": 1.0, "m179[C_mi]": 1.0},
              "bounds": [0.0, 1000.0], "gpr": "YALI1F34729g", "ec": "2.3.1.181", "rhea": "17665",
              "evidence": "Lip2 is required for lipoylation of Lat1, Kgd2 and Gcv3 in S. cerevisiae (Schonauer et al. 2009). "
                          "YALI1F34729g: 30% identity / 84% query coverage to Lip2 (Q06005). Reciprocity is family-level only: the reverse "
                          "top hits are LipB-family octanoyltransferases (S. pombe O36017 first, a bacterial LipB second, K. lactis "
                          "O13476 third); S. cerevisiae Lip2 is not among the top five. Not on the consensus essential list "
                          "(Cas9 essential q 0.0023, Cas12a nonessential q 0.68). Octanoyl is placed on the H-protein "
                          "first by analogy to the Rhea/B. subtilis route; the yeast step order is not established. Returns ACP1 "
                          "in the mitochondrion.",
              "evidence_tier": "B (weak: 30% identity)", "sources": [SRC["schonauer2009"], SRC["rhea17665"]]},
    "R2287": {"name": "lipoyl synthase, mitochondrial (LIP5, lumped sulfur and electron donors)",
              "stoichiometry": {"m2063[C_mi]": -1.0, "m60[C_mi]": -2.0, "m2066[C_mi]": -2.0, "m176[C_mi]": -4.0,
                                "m588[C_mi]": 1.0, "m808[C_mi]": 2.0, "m2067[C_mi]": 2.0, "m756[C_mi]": 2.0, "m178[C_mi]": 4.0},
              "bounds": [0.0, 1000.0], "gpr": "YALI1E13472g", "ec": "2.8.1.8", "rhea": "16585",
              "evidence": "Lip5 is required for lipoylation in S. cerevisiae (Schonauer et al. 2009). YALI1E13472g: 70% identity / "
                          "82% coverage to Lip5 (P32875); reverse top hit Swiss-Prot Q6C6C6 LIPA_YARLI (100%; PE3, HAMAP-inferred, the "
                          "same protein in strain CLIB122, so not independent identity evidence). Rhea 16585 writes the dihydrolipoyl-"
                          "protein product, 2 SAM to 2 Met + 2 5'-deoxyadenosine, with an auxiliary [4Fe-4S] cluster as sulfur source "
                          "(released as 2 hydrogen sulfide + 4 Fe3+) and ferredoxin. Lumped here: the 2 inserted sulfur atoms as "
                          "L-cysteine to L-alanine (cysteine desulfurase) and the electrons as NADPH. The chemistry needs 2 NADPH "
                          "(4 electrons) with ChEBI charges; it is written as 4 NADPH because the model's NADP(+) (m178, "
                          "C21H29N7O17P3) carries charge 0, so each NADPH to NADP(+) releases one H here. Written charge-balanced.",
              "evidence_tier": "B", "sources": [SRC["schonauer2009"], SRC["rhea16585"]]},
    "R2288": {"name": "LIP3-dependent lipoylation of E2 lipoyl domains, mitochondrial (lumped; mechanism unresolved)",
              "stoichiometry": {"m586[C_mi]": -1.0, "m2064[C_mi]": -1.0, "m2062[C_mi]": 1.0, "m2065[C_mi]": 1.0},
              "bounds": [0.0, 1000.0], "gpr": "YALI1C04991g and YALI1E00574g", "ec": "", "rhea": "",
              "evidence": "Lip3 and lipoylated Gcv3 are required for Lat1 and Kgd2 lipoylation in S. cerevisiae (Schonauer et al. "
                          "2009). Mechanism unresolved: Hermes and Cronan 2013 found Lip3 octanoyl-CoA:protein transferase activity "
                          "but no amidotransfer and no ligase activity in vitro, and proposed that Lip3 octanoylates PDH and OGDH "
                          "from octanoyl-CoA downstream of Gcv3 lipoylation; UniProt P47051 is PE3, a putative lipoate-protein ligase. "
                          "This reaction lumps E2 lipoylation as transfer of the lipoyl moiety made on the H-protein; per E2 lipoyl "
                          "the net cost (one mtFAS octanoyl and two sulfur insertions) is the same under both proposed mechanisms. "
                          "No EC or Rhea identifier is asserted. YALI1C04991g: 37% identity / 70% coverage to Lip3/AIM22 (P47051); "
                          "reverse top hits LPLA S. pombe and S. cerevisiae. YALI1E00574g (GCV3 ortholog, 39%/98%, reverse top hit "
                          "GCSH_YEAST) is AND-ed because lipoylated Gcv3 is required. One pooled E2 lipoyl-domain species stands "
                          "for Lat1 and Kgd2, and it is not connected to the PDH or OGDH reactions.",
              "evidence_tier": "C (mechanism unresolved)", "sources": [SRC["schonauer2009"], SRC["hermes2013"]]},
    "R2289": {"name": "L-cysteine transport, cytosol to mitochondrion (ASSUMPTION)",
              "stoichiometry": {"m443[C_cy]": -1.0, "m2066[C_mi]": 1.0}, "bounds": [0.0, 1000.0], "gpr": "", "ec": "", "rhea": "",
              "evidence": "ASSUMPTION: supplies the lipoyl synthase sulfur donor; no mitochondrial cysteine transporter is "
                          "assigned and none was reviewed. Unidirectional.",
              "evidence_tier": "assumption", "sources": []},
    "R2290": {"name": "5'-deoxyadenosine sink, mitochondrion (ASSUMPTION)",
              "stoichiometry": {"m2067[C_mi]": -1.0}, "bounds": [0.0, 1000.0], "gpr": "", "ec": "", "rhea": "",
              "evidence": "ASSUMPTION: 5'-deoxyadenosine from lipoyl synthase has no salvage route in the model; this labelled "
                          "sink removes it. Its flux equals 2 x the lipoyl synthase flux.",
              "evidence_tier": "assumption", "sources": []},
}
new_genes = {
    "YALI1F34729g": {"refseq": "XP_505943.3", "function": "LIP2-type octanoyltransferase candidate (LipB family; family-level reciprocity only)", "tier": "B (weak)"},
    "YALI1E13472g": {"refseq": "XP_503786.1", "function": "lipoyl synthase, LIP5 ortholog (reverse top hit Q6C6C6 LIPA_YARLI, PE3)", "tier": "B"},
    "YALI1C04991g": {"refseq": "XP_501412.3", "function": "LIP3 ortholog (lipoate-protein ligase family; yeast Lip3 shows octanoyl-CoA:protein transferase activity in vitro)", "tier": "C"},
}
edited = {}
for rid in ["R362", "R364"]:
    r = M.reactions.get_by_id(rid)
    before = {m.id: float(c) for m, c in r.metabolites.items()}
    after = {k: v for k, v in before.items() if k != "m28[C_mi]"}
    edited[rid] = {"before": {"stoichiometry": before, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule},
                   "after": {"stoichiometry": after, "bounds": list(r.bounds), "gpr": r.gene_reaction_rule},
                   "reason": "H+ removed: balanced once lipoyl-H-protein m586 has the formula C8H13OS2R (charge 0)."}
used = {m for r in new_rxns.values() for m in r["stoichiometry"]} | {m for e in edited.values() for m in e["before"]["stoichiometry"]}
existing = sorted((used | {"m586[C_mi]", "m179[C_mi]"}) - set(new_mets))
species = {mid: {"formula": g(mid).formula, "charge": g(mid).charge, "compartment": g(mid).compartment} for mid in existing}
residuals = {
    "R1394": {"charge": -1.0, "H": -1.0, "R": -1.0},
    "R1395": {"charge": -1.0, "H": -3.0, "C": -2.0, "O": -1.0},
    "R1396": {"charge": -1.0, "H": -14.0, "C": -8.0, "O": -2.0, "S": -1.0},
    "R2218": {"charge": -3.0, "H": -20.0, "C": -9.0, "O": -1.0, "P": 1.0, "S": -1.0, "R": -1.0, "A": 1.0},
    "R2221": {"charge": -3.0, "H": -5.0, "C": -2.0},
    "R2224": {"charge": -3.0, "H": -5.0, "C": -2.0},
    "R2227": {"charge": -4.0, "H": -3.0, "C": 138.0, "P": -1.0, "A": -1.0},
    "R2230": {"charge": -3.0, "H": -5.0, "C": -2.0},
    "R2233": {"charge": -4.0, "H": -3.0, "C": 158.0, "P": -1.0, "A": -1.0},
}
import json as _json
_before = _json.load(open("model/reports/lipoylation_module_20261020/before.json"))["reactions"]
residuals_before = {rid: _before[rid]["mass_balance"] for rid in
                    ["R78", "R539", "R1394", "R1395", "R1396", "R2218", "R2221", "R2224", "R2227", "R2230", "R2233"]}
spec = {
    "schema_version": 1, "change_id": "LIPOYLATION-MODULE-20261020", "status": "user_authorized_candidate",
    "reviewer": "human_user",
    "authorization": "2026-10-20: user approved the opt-in lipoylation module (D1a fix m586, D2a cysteine donor + NADPH + "
                     "labelled 5'-deoxyadenosine sink, D3 alpha 1e-5, D4 GCV3 AND in LA3 only) and the ACP1 m179 formula "
                     "fill (HSR, charge 0); R2216-R2241 left unchanged.",
    "report": "model/reports/lipoylation_module_20261020/REPORT.md",
    "biomass_reaction": "biomass_C", "alpha_mmol_per_gDW": A,
    "alpha_source": "user-selected provisional value; no Yarrowia lipolytica measurement found; sensitivity at 1e-6, 1e-5, 1e-4",
    "biomass_terms": {"m586[C_mi]": -A, "m2062[C_mi]": A, "m2065[C_mi]": -A, "m2064[C_mi]": A},
    "formula_fills": {
        "m179[C_mi]": {"before": {"formula": None, "charge": 0}, "after": {"formula": "HSR", "charge": 0},
                       "reason": "holo-ACP thiol, consistent with octanoyl-ACP C8H15OSR; balances R78 and R539"},
        "m586[C_mi]": {"before": {"formula": "C12H18N2O4S2R2(C2H2NOR)n", "charge": 0, "name": "lipoylprotein_C12H18N2O4S2R2(C2H2NOR)n"},
                       "after": {"formula": "C8H13OS2R", "charge": 0, "name": "lipoylprotein_C8H13OS2R"},
                       "reason": "lipoyl-lysyl-H-protein, consistent with dihydrolipoyl m588 C8H15OS2R and aminomethyl m587 C9H18NOS2R"}},
    "formula_fill_users": {"m179[C_mi]": ["R78", "R539", "R1394", "R1395", "R1396", "R2218", "R2221", "R2224", "R2227", "R2230", "R2233"],
                           "m586[C_mi]": ["R362", "R364"]},
    "balanced_after": ["R78", "R539", "R362", "R364", "R2286", "R2287", "R2288", "R2289"],
    "preexisting_residuals_after_fill": residuals,
    "acp1_user_residuals_before_fill": residuals_before,
    "edited_reactions": edited, "new_metabolites": new_mets, "new_reactions": new_rxns, "new_genes": new_genes,
    "context_genes": ["YALI1E00574g"], "species": species,
}
json.dump(spec, open("model/curation/lipoylation/lipoylation_module.json", "w"), indent=1)
print("written; existing species recorded:", len(species), existing)
