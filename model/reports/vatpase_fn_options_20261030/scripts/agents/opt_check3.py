import pathlib, cobra
p=pathlib.Path("/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad/vat/m059.xml")
m=cobra.io.read_sbml_model(p)
for rid in ["R302","R384"]:
    r=m.reactions.get_by_id(rid); print(rid,r.name,r.bounds,r.build_reaction_string(use_metabolite_names=True),r.gene_reaction_rule,r.check_mass_balance())
fe=[x for x in m.metabolites if x.elements.get("Fe") and len(x.elements)==1]
print("pure Fe species:",[(x.id,x.name,x.compartment) for x in fe])
