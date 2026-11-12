"""Write the R2121 closure curation file from the reviewed e584151f state (static read; the decision text is below)."""
import hashlib, json, sys, pathlib, cobra
model_p, out_p = sys.argv[1:3]
sha = hashlib.sha256(open(model_p, 'rb').read()).hexdigest()
assert sha == 'e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2', sha
m = cobra.io.read_sbml_model(pathlib.Path(model_p))
def state(rid):
    x = m.reactions.get_by_id(rid)
    return ({k.id: v for k, v in x.metabolites.items()},
            {k.id: {'formula': k.formula, 'charge': k.charge, 'compartment': k.compartment} for k in x.metabolites})
r, k = m.reactions.R2121, m.reactions.R2119
s21, sp21 = state('R2121'); s19, sp19 = state('R2119')
spec = {
 "schema_version": 1,
 "curation_id": "R2121-MALONYL-COA-CLOSE-20261112",
 "status": "user_authorized",
 "authorization": "2026-11-12 user, after the Stage 1 results (model/reports/r2121_malonyl_coa_20261110/REPORT.md): \"Approve Stage 2\" = close R2121 to (0, 0), keeping the reaction and its GPR, with a note giving the KEGG evidence, through a guarded --r2121-malonyl-coa-close step. Stage 1 was approved on 2026-11-10 (\"A2 malonyl-CoA first\", \"Approve Stage 1\").",
 "report": "model/reports/r2121_malonyl_coa_20261110/",
 "reviewer": "Claude (reaction-activation-options session, 2026-11-10/12); independent audit in the report",
 "baseline_model_sha256": sha,
 "baseline_model_note": "Provenance only; the guard pins the reaction states below, not the file SHA.",
 "reaction_id": "R2121",
 "reaction_name": r.name,
 "gpr": r.gene_reaction_rule,
 "before_bounds": [-1000, 1000],
 "after_bounds": [0, 0],
 "stoichiometry": s21,
 "species": sp21,
 "annotation": r.annotation,
 "before_notes": dict(r.notes),
 "balance_note": "As inherited: C and H balance; charge -1 on the substrate side (malonate semialdehyde m1976[C_cy] is -1) and 0 on the product side. Not changed: closing the reaction does not touch its chemistry.",
 "retained_reaction": {
  "reaction_id": "R2119",
  "reaction_name": k.name,
  "gpr": k.gene_reaction_rule,
  "bounds": list(k.bounds),
  "annotation": k.annotation,
  "stoichiometry": s19,
  "species": sp19,
  "role": "KEGG R00705 (EC 1.2.1.18 / 1.2.1.27), the reaction of the K00140 ortholog group that YALI1C02550g belongs to; kept unchanged"
 },
 "evidence": [
  {"claim": "R2121 is KEGG R00740, EC 1.2.1.75, malonyl-CoA reductase (malonate semialdehyde-forming) of the 3-hydroxypropanoate CO2-fixation cycles of green non-sulfur bacteria and thermoacidophilic archaea",
   "source_url": "https://www.kegg.jp/entry/R00740; https://www.kegg.jp/entry/1.2.1.75", "evidence_tier": "verified (KEGG REST read 2026-11-10)"},
  {"claim": "YALI1C02550g is yli:2909116, ortholog K00140 (EC 1.2.1.18 / 1.2.1.27), whose reaction R00705 (decarboxylating, CoA-acetylating) is R2119 here",
   "source_url": "https://www.kegg.jp/entry/yli:2909116; https://www.kegg.jp/entry/R00705", "evidence_tier": "verified (KEGG REST read 2026-11-10; model annotation ec-code 1.2.1.27)"},
  {"claim": "KEGG assigns no K15017 (malonyl-CoA reductase) gene to Y. lipolytica (yli) or S. cerevisiae (sce)",
   "source_url": "https://www.kegg.jp/entry/K15017", "evidence_tier": "verified, annotation-level absence (KEGG REST link query 2026-11-10); not an enzyme assay"},
  {"claim": "Yeast-GEM carries R00705 (r_4573, irreversible) and not R00740; eciYali5-GEM carries neither",
   "source_url": "$IYALI26_RESEARCH_ROOT/reference/external_models/yeast_gem/Yeast-GEM.xml (SHA-256 9fd2c572...); chalmers_yarrowia/eciYali5-GEM/model/eciYali5-GEM.xml (628df440...)", "evidence_tier": "verified (local copies scanned by KEGG id and name, 2026-11-10)"},
  {"claim": "Acetyl-CoA carboxylase ACC1 is essential in S. cerevisiae",
   "source_url": "https://pubmed.ncbi.nlm.nih.gov/8098706/", "evidence_tier": "verified (abstract)"},
  {"claim": "R2121 is inherited unchanged (reversible, YALI1C02550g) from iYli21 and the immutable source",
   "source_url": "model/source/iyli21.xml; model/source/iyali26.xml", "evidence_tier": "verified (read 2026-11-10)"},
  {"claim": "Stage 1 in memory on e584151f: WT -0.139 %; R88 required (FVA minimum 0.0077); no-uptake ATP maximum 0; sweep 15 % 141/185/43/694 to 143/183/43/694",
   "source_url": "model/reports/r2121_malonyl_coa_20261110/REPORT.md", "evidence_tier": "observed (computation, 2026-11-10)"}
 ],
 "notes": {
  "r2121_malonyl_coa_close_curation": "R2121-MALONYL-COA-CLOSE-20261112",
  "r2121_previous_bounds": "[-1000, 1000]",
  "r2121_malonyl_coa_close_evidence": "KEGG R00740 / EC 1.2.1.75 is the malonyl-CoA reductase (malonate semialdehyde-forming) of the 3-hydroxypropanoate CO2-fixation cycles of green non-sulfur bacteria and thermoacidophilic archaea. YALI1C02550g is KEGG yli:2909116, K00140 (EC 1.2.1.18 / 1.2.1.27); its reaction R00705 is R2119 in this model, kept. KEGG lists no K15017 gene for Y. lipolytica or S. cerevisiae; Yeast-GEM carries R00705 and not R00740, eciYali5-GEM neither (checked 2026-11-10).",
  "r2121_malonyl_coa_close_effect": "Stage 1 on e584151f (Task A medium + PO1f, Gurobi 1e-9): WT growth -0.139 %; cytosolic malonyl-CoA moves to the acetyl-CoA carboxylase R88; a pFBA TAG lipolysis, beta-oxidation and elongation cycle on glucose stops.",
  "r2121_malonyl_coa_close_limits": "Inherited from iYli21 with no recorded reason. Absence from KEGG and from the yeast models is annotation-level evidence, not an enzyme assay; no Y. lipolytica enzyme data. Chemistry unchanged (charge -1 on the substrate side as inherited). The inherited ec-code annotation 1.2.1.27 is the EC of the K00140 group of YALI1C02550g, not of R00740 (EC 1.2.1.75); annotations are not changed here.",
  "r2121_malonyl_coa_close_sources": "https://www.kegg.jp/entry/R00740; https://www.kegg.jp/entry/1.2.1.75; https://www.kegg.jp/entry/yli:2909116; https://www.kegg.jp/entry/R00705; https://www.kegg.jp/entry/K15017",
  "r2121_malonyl_coa_close_authorization": "User authorized closing R2121 to (0, 0), reaction and GPR kept, on 2026-11-12 after the Stage 1 report (model/reports/r2121_malonyl_coa_20261110/)."
 }
}
for v in spec['notes'].values():
    assert not any(c in v for c in '<>&\'"'), v  # SBML escapes these on write
json.dump(spec, open(out_p, 'w'), indent=2, ensure_ascii=False); open(out_p, 'a').write('\n')
print('written', out_p, 'before_notes', list(spec['before_notes']))
