import pathlib, cobra
p=pathlib.Path("/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad/vat/m059.xml")
m=cobra.io.read_sbml_model(p)
r=m.reactions.R302
for met in r.metabolites:
    if met.name.startswith("iron"):
        print(met.id,met.name,met.compartment,met.formula,met.charge)
        for rr in met.reactions: print("  ",rr.id,rr.name,rr.bounds,rr.build_reaction_string(),rr.gene_reaction_rule)
        # trace other compartments with same name
        for x in m.metabolites:
            if x.name==met.name and x is not met:
                print(" other",x.id,x.compartment,[(q.id,q.bounds,q.build_reaction_string(),q.gene_reaction_rule) for q in x.reactions])
