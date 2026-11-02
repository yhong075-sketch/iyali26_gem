import pathlib, cobra
m = cobra.io.read_sbml_model(pathlib.Path('../m059.xml'))
def show(rid):
    r=m.reactions.get_by_id(rid)
    print(' ',rid,'|',r.name,'|',r.bounds,'|',r.build_reaction_string(use_metabolite_names=False),'| GPR:',r.gene_reaction_rule,'| bal:',r.check_mass_balance())
for rid in ['R2047','R2046','R1255','R1359','R2048','R1240']: show(rid)
for mid in ['m676[C_er]']:
    x=m.metabolites.get_by_id(mid); print(mid,x.name)
    for r in x.reactions: show(r.id)
for rid in ['R2047','R2046','R1255']:
    for x in m.reactions.get_by_id(rid).metabolites:
        if x.compartment!='C_go': print(' other-side species',x.id,x.name,[ (r.id,r.bounds) for r in x.reactions][:12])
print('mass balance R2048:',m.reactions.R2048.check_mass_balance())
for mid in ['m1327[C_go]','m1886[C_go]','m1314[C_go]','m293[C_go]']:
    x=m.metabolites.get_by_id(mid); print(mid,x.formula,x.charge)
