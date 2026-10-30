import pathlib, cobra
m = cobra.io.read_sbml_model(pathlib.Path("../m059.xml"))
b = m.reactions.biomass_C
print("biomass_C:", b.build_reaction_string())
comps = {x.compartment for x in b.metabolites}; print("biomass comps", comps)
# expand one level into pseudo-reactions producing biomass_C reactants
for x in b.metabolites:
    if b.metabolites[x] < 0:
        prods = [r for r in x.reactions if r.metabolites[x] > 0 or r.reversibility]
        for r in prods:
            cc = {y.compartment for y in r.metabolites}
            if cc & {"C_va","C_go"}: print("  biomass precursor touches va/go:", x.id, r.id, r.name, cc)
# all biomass-like reactions: any va/go metabolites?
for rid in ["biomass_C","xBIOMASS","R1372","R1387","R1710","R1373"]:
    r=m.reactions.get_by_id(rid); cc={y.compartment for y in r.metabolites}
    print(rid, r.bounds, cc, "obj" if r.objective_coefficient else "")
# vacuolar solutes of H+-coupled transporters: other reactions
H=m.metabolites.get_by_id("m1007[C_va]")
for r in sorted(H.reactions,key=lambda r:r.id):
    if r.id in("R795","R1163"): continue
    for x in r.metabolites:
        if x.compartment=="C_va" and x.id!="m1007[C_va]":
            others=[(o.id,o.bounds) for o in x.reactions if o.id!=r.id]
            print(f"{r.id} {r.bounds} solute {x.id} {x.name} -> other rxns {others}")
# Golgi: all reactions crossing the Golgi boundary
print("\nGolgi boundary-crossing reactions:")
for r in m.reactions:
    cc={y.compartment for y in r.metabolites}
    if "C_go" in cc and len(cc)>1:
        print(" ", r.id, r.name, r.bounds, r.build_reaction_string(), "| GPR:", r.gene_reaction_rule or "-")
# endosome H+ cross-membrane reactions
print("\nEndosome H+ cross-membrane:")
for r in m.metabolites.get_by_id("m24[C_en]").reactions:
    cc={y.compartment for y in r.metabolites}
    if len(cc)>1: print(" ", r.id, r.name, r.bounds, r.build_reaction_string())
# any other ATP-driven H+ translocators anywhere
print("\nATP hydrolysis + H+ translocation (any membrane):")
Hs={x.id for x in m.metabolites if x.formula=="H" and x.charge==1}
for r in m.reactions:
    atp=[x for x in r.metabolites if x.name.startswith("ATP") and r.metabolites[x]<0]
    hs=[x for x in r.metabolites if x.id in Hs]
    if atp and len({x.compartment for x in hs})>1:
        print(" ", r.id, r.name, r.bounds, r.build_reaction_string(), "|", r.gene_reaction_rule[:80])
