import pathlib, cobra, sys
m=cobra.io.read_sbml_model(pathlib.Path(sys.argv[1]))
hits=[x for x in m.metabolites if 'iron' in (x.name or '').lower() or 'ferr' in (x.name or '').lower() or 'sideroph' in (x.name or '').lower() or 'hemin' in (x.name or '').lower() or x.id=='m948']
for x in hits:
    print(x.id,repr(x.name),x.compartment,repr(x.formula),x.charge)
    for rr in x.reactions: print('    ',rr.id,rr.bounds,rr.build_reaction_string(True)[:160],'|',rr.gene_reaction_rule[:80])
