import pathlib, cobra, re, csv
from cobra.core.gene import GPR
m = cobra.io.read_sbml_model(pathlib.Path('../m059.xml'))
def show(rid):
    r=m.reactions.get_by_id(rid)
    print(' ',rid,'|',r.name,'|',r.bounds,'|',r.build_reaction_string(),'| GPR:',r.gene_reaction_rule)
# medium exchanges: Fe content
print('== Task A medium exchange metabolites with Fe or heme/siderophore names')
med='/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/essential-gene-prediction-e76be3/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv'
for row in csv.DictReader(open(med)):
    r=m.reactions.get_by_id(row['exchange'])
    for x in r.metabolites:
        if re.search(r'Fe(?![a-z])',x.formula or '') or re.search('iron|heme|hemin|sideroph|ferri|ferro',x.name,re.I):
            print(' ',r.id,x.id,x.name,x.formula,row['uptake'])
# any exchange-like reaction (single metabolite, C_ex) touching Fe/heme anywhere
print('== all boundary reactions on Fe/heme species')
for r in m.boundary:
    for x in r.metabolites:
        if re.search(r'Fe(?![a-z])',x.formula or '') or re.search('iron|heme|hemin|sideroph',x.name,re.I): print(' ',r.id,r.bounds,x.id,x.name)
# Golgi H+
print('== m293 reactions'); [show(r.id) for r in m.metabolites.get_by_id('m293[C_go]').reactions]
# Golgi boundary reactions carrying any H+ species
hids={x.id for x in m.metabolites if x.formula=='H' and x.charge==1}
print('H+ species:',sorted(hids))
gb=[r for r in m.reactions if any(x.compartment=='C_go' for x in r.metabolites) and any(x.compartment!='C_go' for x in r.metabolites)]
print('Golgi-boundary reactions:',len(gb),'; carrying H+ :',[r.id for r in gb if any(x.id in hids for x in r.metabolites)])
# R2048 trace
print('== R2048 and its species')
show('R2048')
for x in m.reactions.get_by_id('R2048').metabolites:
    print(' species',x.id,x.name,x.compartment,'->',[(r.id,r.bounds,r.get_coefficient(x.id)) for r in x.reactions])
for rid in ['R1240','R1159','R1163','R586','R_NDP7g','R141','R_NTP1','R795','R794']: show(rid)
# biomass compartments and sphingolipid presence
b=m.reactions.get_by_id('biomass_C')
print('biomass compartments',sorted({x.compartment for x in b.metabolites}))
print('biomass sphingo/ceramide/IPC comps',[x.name for x in b.metabolites if re.search('ceramide|IPC|inositol.?phospho.?ceramide|sphing',x.name,re.I)])
# vacuolar cysteine
print('vacuolar cysteine species',[ (x.id,x.name) for x in m.metabolites if x.compartment=='C_va' and re.search('cyste',x.name,re.I)])
for rid in ['R2020','R2028','R2030','R2038','R2040','R882']: show(rid)
# Boolean O1 rule
core=['YALI1A09766g','YALI1A11258g','YALI1D00581g','YALI1E10492g','YALI1F20965g','YALI0E16192g','YALI1E14125g','YALI1F31854g','YALI1E37063g','YALI1F13017g','YALI1F21690g','YALI1E32332g']
o1=GPR.from_string('('+' and '.join(core)+') and (YALI1E12482g or YALI1F38820g)')
genes14=core+['YALI1E12482g','YALI1F38820g']
print('O1 rule single-KO closes:',[g for g in genes14 if not o1.eval({g})])
print('O1 rule single-KO keeps open:',[g for g in genes14 if o1.eval({g})])
# O2a coupled R2193 rule
r2193=m.reactions.get_by_id('R2193').gene_reaction_rule
cpl=GPR.from_string('('+r2193+') and ('+' and '.join(core)+') and (YALI1E12482g or YALI1F38820g)')
print('O2a coupled R2193 single-KO closes:',[g for g in genes14+['YALI1D09277g','YALI1D09308g','YALI1D08684g','YALI1D08564g'] if not cpl.eval({g})])
cons={row[0] for row in csv.reader(open('/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/state/essentiality/repository/consensus_essential_genes.csv'))}
print('closers on consensus:',[g for g in genes14 if not cpl.eval({g}) and g in cons]); print('closers NOT on consensus:',[g for g in genes14 if not cpl.eval({g}) and g not in cons])
# other reactions carrying the 14 genes
for g in genes14:
    print(' gene',g,'reactions',[r.id for r in m.genes.get_by_id(g).reactions])
