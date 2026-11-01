import pathlib, cobra, hashlib
p=pathlib.Path("/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad/vat/m059.xml")
print(hashlib.sha256(p.read_bytes()).hexdigest())
m=cobra.io.read_sbml_model(p)
print("genes",len(m.genes),"rxns",len(m.reactions))
for rid in ["R1189","R2193","R863","R727"]:
    r=m.reactions.get_by_id(rid); print(rid,r.name,r.bounds,r.build_reaction_string(),"| GPR:",r.gene_reaction_rule)
fe=[x for x in m.metabolites if x.formula in ("Fe",) or "iron" in (x.name or "").lower() or "Fe" == (x.formula or "")]
for x in fe:
    print("MET",x.id,x.name,x.compartment,x.formula,x.charge,[ (r.id,r.bounds, r.metabolites[x]) for r in x.reactions])
b=m.reactions.get_by_id("biomass_C")
for x,c in b.metabolites.items():
    if x.id in ("m948",) or (x.formula and "Fe" in x.formula): print("BIOMASS",x.id,x.name,c)
# all metabolites with Fe in formula: producers
for x in m.metabolites:
    if x.formula and "Fe" in x.formula and x not in fe:
        print("FEFORM",x.id,x.name,x.compartment,x.formula,[r.id for r in x.reactions][:12])
