import pathlib, cobra
m = cobra.io.read_sbml_model(pathlib.Path("../m059.xml"))
def show(mid):
    x=m.metabolites.get_by_id(mid); print("##",mid,x.name)
    for r in x.reactions: print("   ",r.id,r.name,r.bounds,r.build_reaction_string(),"|",r.gene_reaction_rule[:60] or "-")
for mid in ["m836[C_go]","m837[C_go]","m838[C_go]","m1314[C_go]","m1327[C_go]","m1886[C_go]","m289[C_go]","m291[C_go]","m292[C_go]","m290[C_go]","m847[C_go]","m1889[C_go]"]:
    show(mid)
show("m676[C_er]")
