import pathlib, cobra, csv
p=pathlib.Path("/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad/vat/m059.xml")
m=cobra.io.read_sbml_model(p)
cons=open("/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/state/essentiality/repository/consensus_essential_genes.csv").read()
r=m.reactions.R2193
print(r.name, r.notes if r.notes else "")
for g in r.genes:
    print(g.id, repr(g.name), "CONSENSUS" if g.id in cons else "-", [x.id for x in g.reactions])
    for line in cons.splitlines():
        if g.id in line: print("   ",line[:250])
print([(x.id,x.name,x.compartment) for x in r.metabolites])
