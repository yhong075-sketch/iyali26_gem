import pathlib, re, cobra, hashlib
p=pathlib.Path("/private/tmp/claude-501/-Users-david-Desktop-Lab-Ian-wheeldon-code-iyali26-gem--claude-worktrees-essential-gene-prediction-e76be3/0a4ec4c9-d926-4768-b75c-c14c6fe4ecec/scratchpad/vat/m059.xml")
print(hashlib.sha256(p.read_bytes()).hexdigest())
m=cobra.io.read_sbml_model(p)
bio=m.reactions.get_by_id("biomass_C") if "biomass_C" in m.reactions else None
# iron / copper / calcium / zinc / cysteine species by formula, any compartment
def show(pred,label):
    print("==",label)
    for met in m.metabolites:
        if pred(met):
            rx=[(r.id,r.bounds) for r in met.reactions]
            print(met.id,met.name,met.compartment,met.formula,len(rx),rx[:8])
show(lambda x: x.formula in ("Fe","Cu","Ca","Zn") , "Fe/Cu/Ca/Zn ions")
show(lambda x: x.compartment in ("C_va",) and "ystein" in (x.name or ""), "vacuolar cysteine")
# biomass components with Fe
if bio:
    print("== biomass Fe/S/heme-containing comps")
    for met,c in bio.metabolites.items():
        if met.formula and ("Fe" in met.formula):
            print(met.id,met.name,met.compartment,met.formula,c)
# genes named FET3/FTR1 like
for g in m.genes:
    if re.search(r"FET|FTR|ferroxidase|iron permease|multicopper",g.name or "",re.I):
        print("gene",g.id,g.name,[r.id for r in g.reactions])
print("R_NTP1",m.reactions.R_NTP1.bounds,m.reactions.R_NTP1.reaction, m.reactions.R_NTP1.gene_reaction_rule)
for rid in ["xMAINTENANCE"]:
    if rid in m.reactions: r=m.reactions.get_by_id(rid); print(rid,r.bounds,r.reaction)
