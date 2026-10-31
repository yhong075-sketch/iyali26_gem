import pathlib, csv, cobra
p=pathlib.Path('../m059.xml').resolve()
m=cobra.io.read_sbml_model(p)
for r in ['R794','R795','R1159','R1163','R1363','R871','R876','R882','R883','R_NTP1','R2018','R2021']:
    if r in m.reactions:
        x=m.reactions.get_by_id(r)
        print(r, x.bounds, x.build_reaction_string(use_metabolite_names=False), '| GPR:', x.gene_reaction_rule[:400])
    else: print(r,'absent')
bio=m.reactions.get_by_id('biomass_C')
print('biomass compartments', sorted({mm.compartment for mm in bio.metabolites}))
genes=['YALI1A11258g','YALI1E14125g','YALI1D00581g','YALI1E10492g','YALI1A09766g','YALI1F20965g','YALI0E16192g','YALI1F31854g','YALI1F21690g','YALI1E37063g','YALI1E32332g','YALI1E12482g','YALI1F38820g','YALI1F13017g']
cons='/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/state/essentiality/repository/consensus_essential_genes.csv'
txt=open(cons).read()
print('consensus header', txt.splitlines()[0][:200])
for g in genes:
    rx=[r.id for r in m.genes.get_by_id(g).reactions] if g in m.genes else 'not in model'
    print(g, 'in_consensus' if g in txt else 'NOT_in_consensus', rx)
# boolean: which single KO makes R794/R795 GPR false
from cobra.core.gene import GPR
for rid in ['R794','R795']:
    gpr=m.reactions.get_by_id(rid).gpr
    allg={g.id for g in m.reactions.get_by_id(rid).genes}
    print(rid,'single KOs making GPR false:',[g for g in sorted(allg) if not gpr.eval(knockouts={g})])
