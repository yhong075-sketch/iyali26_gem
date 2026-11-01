import pathlib, cobra
p=pathlib.Path("/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad/vat/m059.xml")
m=cobra.io.read_sbml_model(p)
for rid in ["R2048","R1240","R1189","R727"]:
    r=m.reactions.get_by_id(rid); print(rid,r.name,r.bounds,r.build_reaction_string(),r.gene_reaction_rule)
for met in m.reactions.R2048.metabolites:
    print(met.id,met.name[:30],met.compartment,[(x.id,x.bounds) for x in met.reactions])
