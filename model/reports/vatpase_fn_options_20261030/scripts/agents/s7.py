import pathlib, cobra
m = cobra.io.read_sbml_model(pathlib.Path("../m059.xml"))
mip = [x for x in m.metabolites if "phosphomannosylinositol" in x.name.lower() or x.name.lower().startswith("m(ip)2c")]
for x in mip:
    print("##",x.id,x.name[:60],x.compartment)
    for r in x.reactions:
        print("    ",r.id,r.name[:40],r.bounds,"coef",r.metabolites[x],"|",r.gene_reaction_rule[:40] or "-")
# biomass lipid components names
b=m.reactions.biomass_C
for x in b.metabolites:
    if b.metabolites[x]<0 and not x.id.startswith("trna"):
        print("BIO",x.id,x.name[:50],b.metabolites[x])
