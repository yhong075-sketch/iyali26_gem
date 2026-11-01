import pathlib, re, cobra, hashlib, sys
p=pathlib.Path(sys.argv[1])
print(hashlib.sha256(p.read_bytes()).hexdigest())
m=cobra.io.read_sbml_model(p)
r=m.reactions.get_by_id('R2193'); print('R2193',r.name,r.bounds,r.build_reaction_string(True),'|',r.gene_reaction_rule)
fe=[x for x in m.metabolites if x.formula and re.search(r'Fe(?![a-z])',x.formula)]
print('Fe-containing species',len(fe))
for x in fe:
    print(' ',x.id,x.name,x.compartment,x.formula,[(rr.id,rr.bounds) for rr in x.reactions][:8])
for rr in m.exchanges:
    if any(x in fe for x in rr.metabolites): print('EX',rr.id,rr.name,rr.bounds)
