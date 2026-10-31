import pathlib, cobra, re
m = cobra.io.read_sbml_model(pathlib.Path('../m059.xml'))
def show(rid):
    r=m.reactions.get_by_id(rid)
    print(rid,'|',r.name,'|',r.bounds,'|',r.build_reaction_string(),'| GPR:',r.gene_reaction_rule)
for rid in ['R2193','R1189','R863','R727','R302','R384']: show(rid)
# any metabolite whose formula contains element Fe
print('--- Fe-containing formulas')
for x in m.metabolites:
    f=x.formula or ''
    if re.search(r'Fe(?![a-z])',f): print(x.id,x.name,f,x.compartment,[r.id for r in x.reactions])
# biomass components with Fe or heme
b=m.reactions.get_by_id('biomass_C')
print('--- biomass Fe/heme/siroheme/cyt comps')
for x,c in b.metabolites.items():
    if re.search(r'Fe(?![a-z])',x.formula or '') or re.search('heme|iron|siro|cytochrome',x.name,re.I): print(x.id,x.name,c)
# R2193 Boolean with core coupling
r=m.reactions.get_by_id('R2193'); print('R2193 genes',[g.id for g in r.genes])
for g in r.genes: print(g.id, 'KO closes R2193:', not r.gpr.eval({g.id}))
