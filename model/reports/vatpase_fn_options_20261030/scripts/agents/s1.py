import pathlib, cobra, csv, json
m = cobra.io.read_sbml_model(pathlib.Path("../m059.xml"))
print(m.compartments)
for rid in ["R794","R795","R1159"]:
    r = m.reactions.get_by_id(rid)
    print("==",rid, r.name)
    print(" ", r.build_reaction_string(use_metabolite_names=False))
    print(" ", r.build_reaction_string(use_metabolite_names=True))
    print("  comps", {mm.id: mm.compartment for mm in r.metabolites})
    print("  bounds", r.bounds)
    print("  GPR", r.gene_reaction_rule)
    print("  notes", r.notes)
    print("  annot", r.annotation)
