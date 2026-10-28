import cobra, warnings, sys, re
warnings.simplefilter("ignore")
Y = cobra.io.read_sbml_model(sys.argv[1])
pat = re.compile(r"dolich|prenyl", re.I)
rx = [r for r in Y.reactions if any(pat.search(m.name) for m in r.metabolites)]
for r in sorted(rx, key=lambda r: r.id):
    print(r.id, "|", r.name, "|", r.bounds, "|", r.gene_reaction_rule, "|", r.build_reaction_string(use_metabolite_names=True)[:200])
print("--- biomass/lipid pseudo terms with dolichol ---")
for r in Y.reactions:
    if "pseudoreaction" in r.name.lower() or r.id in ("r_4041","r_4598","r_2108","r_4063"):
        hits = [(m.name, c) for m, c in r.metabolites.items() if pat.search(m.name)]
        if hits: print(r.id, r.name, hits)
print("--- genes ---")
for g in ("YBR002C","YMR101C","YDL193W","YMR013C","YGR036C","YIL049W"):
    if g in Y.genes: print(g, sorted(r.id for r in Y.genes.get_by_id(g).reactions))
    else: print(g, "not in yeast-GEM")
