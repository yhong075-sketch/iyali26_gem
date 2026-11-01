import pathlib, csv, re, cobra, openpyxl
R="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/state/essentiality/repository/"
genes = {  # locus-map S. cerevisiae ortholog annotation -> subunit guess
 "YALI1A09766g":"V1 A (VMA1/TFP1)","YALI1E32332g":"V1 B (VMA2)","YALI1A11258g":"V1 C (VMA5)",
 "YALI1D00581g":"V1 D (VMA8)","YALI1F20965g":"V1 E (VMA4)","YALI1E19360g":"V1 F (VMA7), W29 Pseudo mRNA",
 "YALI0E16192g":"V1 F (VMA7), old CLIB122 ID","YALI1E14125g":"V1 G (VMA10)","YALI1E10492g":"V1 H (VMA13)",
 "YALI1E12482g":"V0 a (STV1/VPH1-like)","YALI1F38820g":"V0 a (VPH1-like)","YALI1F31854g":"V0 c (VMA3/CUP5)",
 "YALI1E37063g":"V0 c' (VMA11/TFP3)","YALI1F13017g":"V0 c'' (VMA16/PPA1)","YALI1F21690g":"V0 d (VMA6)",
 "YALI1F25663g":"V0 e (VMA9)","YALI1B30163g":"assembly VMA21","YALI1A01298g":"assembly VMA12/VPH2",
 "YALI1F37528g":"assembly PKR1","YALI1C03969g":"assembly VOA1 (YGR106C)","YALI1C23873g":"RAVE RAV1"}
m = cobra.io.read_sbml_model(pathlib.Path("m059.xml"))
mg = {g.id for g in m.genes}
for r in ("R794","R795"):
    rx=m.reactions.get_by_id(r); print(r, rx.name, rx.reaction, rx.bounds); print("  GPR:", rx.gene_reaction_rule)
cons = {row["gene_id"]:row for row in csv.DictReader(open(R+"consensus_essential_genes.csv"))}
fit = {}
for row in csv.DictReader(open(R+"cas9_cas12a_fitness.csv")):
    fit.setdefault(row["gene_id"],[]).append(row)
# full-model scan for any other V-ATPase-annotated gene
for g in m.genes:
    txt=(g.name or "")+" "+str(g.annotation)+" "+str(getattr(g,"notes",""))
    if re.search(r"V-ATPase|vacuolar.{0,15}ATP|VMA\d|VPH|H\+-ATPase", txt, re.I) and g.id not in genes:
        print("EXTRA model gene:", g.id, txt[:200])
print()
for g,s in genes.items():
    inm = g in mg
    rx = sorted(r.id for r in m.genes.get_by_id(g).reactions) if inm else []
    c = g in cons
    calls = "; ".join(f"{f['assay']} FS={f['fitness_score']} q={f['q_value']} call={f['experimental_call']} (row {f['source_row']})" for f in fit.get(g,[])) or "no Cas9/Cas12a row"
    print(f"{g}\t{s}\tin_model={inm} {rx}\tconsensus={c}\t{calls}")
