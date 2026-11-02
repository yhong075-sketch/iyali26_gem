import cobra, sys, warnings
from pathlib import Path
warnings.simplefilter("ignore")
m = cobra.io.read_sbml_model(Path(sys.argv[1]))
def show(r):
    print(f"== {r.id} | {r.name} | bounds {r.bounds} | comp {r.compartments}")
    print("   ", r.build_reaction_string(use_metabolite_names=False))
    print("   ", r.build_reaction_string(use_metabolite_names=True))
    print("   GPR:", r.gene_reaction_rule, "| balance:", r.check_mass_balance())
    print("   ann:", {k:v for k,v in r.annotation.items()})
    print("   notes:", r.notes)
keys = ["homocitr","homoacon","homoisocitr","oxoadipate","2-oxoadipate","aminoadipate","methylcitr","methylisocitr","aconit"]
mets=set()
for met in m.metabolites:
    if any(k in (met.name or "").lower() for k in keys):
        mets.add(met)
for met in sorted(mets, key=lambda x:x.id):
    print(f"MET {met.id} {met.name} {met.formula} charge={met.charge} comp={met.compartment} ann={dict(met.annotation)}")
rxns=set()
for met in mets: rxns |= set(met.reactions)
for r in sorted(rxns, key=lambda x:x.id): show(r)
for g in ["YALI1E03387g","YALI1D11984g","YALI1E17935g","YALI1F03803g"]:
    gg=m.genes.get_by_id(g); print("GENE",g, gg.name, dict(gg.annotation), sorted(r.id for r in gg.reactions))
