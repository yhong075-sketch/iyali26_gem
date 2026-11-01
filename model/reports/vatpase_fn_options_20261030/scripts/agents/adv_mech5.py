import pathlib, cobra
m = cobra.io.read_sbml_model(pathlib.Path('../m059.xml'))
def show(rid):
    r=m.reactions.get_by_id(rid)
    print(' ',rid,'|',r.name,'|',r.bounds,'|',r.build_reaction_string(use_metabolite_names=True),'| GPR:',r.gene_reaction_rule)
for mid in ['m836[C_go]','m1635[C_cy]']:
    print(mid); [show(r.id) for r in m.metabolites.get_by_id(mid).reactions]
show('R439'); show('R383')
for mid in ['m390[C_er]','m674[C_er]']:
    x=m.metabolites.get_by_id(mid); print(mid,x.name,[ (r.id,r.bounds,r.get_coefficient(mid)) for r in x.reactions])
b=m.reactions.biomass_C
print('biomass bounds',b.bounds,'objective',[r.id for r in m.reactions if r.objective_coefficient])
print('xMAINTENANCE',m.reactions.xMAINTENANCE.bounds if 'xMAINTENANCE' in m.reactions else 'absent')
# what fraction of iron cap does WT need (arithmetic)
wt=0.3116608216336117
need=3.04e-5*wt+1e-6*wt
print('WT Fe need (m948 direct + heme a) =',need,' cap 1.61e-5; ratio',need/1.61e-5)
for tol in (1e-9,1e-7,1e-6):
    g=tol/3.04e-5; print('tol',tol,'max leak growth',g,'KO/WT',g/wt)
