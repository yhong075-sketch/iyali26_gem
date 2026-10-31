import pathlib, cobra, re, itertools
m = cobra.io.read_sbml_model(pathlib.Path('../m059.xml'))
def show(rid):
    r=m.reactions.get_by_id(rid)
    print(rid, r.name, r.bounds, r.build_reaction_string(), '| GPR:', r.gene_reaction_rule)
# 1) Boolean GPR on R794/R795
for rid in ['R794','R795']:
    r=m.reactions.get_by_id(rid); genes=[g.id for g in r.genes]
    single=[g for g in genes if not r.gpr.eval({g})]
    print(rid,'n genes',len(genes),'single KOs falsifying:',single)
# 2) all iron species
print('--- iron species')
fe=[x for x in m.metabolites if (x.formula and re.fullmatch(r'Fe\d*',x.formula)) or 'iron' in (x.name or '').lower() or x.name.lower().startswith('fe')]
for x in fe:
    print(x.id, x.name, x.formula, x.charge, x.compartment, [ (r.id, r.get_coefficient(x.id), r.bounds) for r in x.reactions])
