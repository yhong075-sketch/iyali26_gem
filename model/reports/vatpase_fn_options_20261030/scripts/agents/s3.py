import pathlib, cobra
m = cobra.io.read_sbml_model(pathlib.Path("../m059.xml"))
H = [x for x in m.metabolites if x.formula == "H" and x.charge == 1]
for h in H: print(h.id, h.name, h.compartment, len(h.reactions))
bio = [r for r in m.reactions if "biomass" in r.name.lower() or "growth" in r.name.lower()]
print("objective:", [r.id for r in m.reactions if r.objective_coefficient])
for r in bio: print("bio-ish", r.id, r.name)
