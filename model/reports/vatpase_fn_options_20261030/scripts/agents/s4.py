import pathlib, cobra
m = cobra.io.read_sbml_model(pathlib.Path("../m059.xml"))
for hid in ["m1007[C_va]","m293[C_go]","m19[C_vm]","m871[C_gm]"]:
    h = m.metabolites.get_by_id(hid)
    print("\n######", hid, len(h.reactions))
    for r in sorted(h.reactions, key=lambda r: r.id):
        comps = sorted({x.compartment for x in r.metabolites})
        print(f"{r.id} | {r.name} | lb={r.lower_bound} ub={r.upper_bound} | H coef={r.metabolites[h]} | comps={comps}")
        print("   ", r.build_reaction_string(use_metabolite_names=False))
        print("   ", r.build_reaction_string(use_metabolite_names=True))
        print("    GPR:", r.gene_reaction_rule or "(none)")
