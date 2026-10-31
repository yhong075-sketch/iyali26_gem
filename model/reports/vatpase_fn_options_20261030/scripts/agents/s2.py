import pathlib, cobra, csv, json
from cobra.core.gene import GPR
m = cobra.io.read_sbml_model(pathlib.Path("../m059.xml"))
ess = {}
with open("/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research/state/essentiality/repository/consensus_essential_genes.csv") as f:
    for row in csv.DictReader(f): ess[row["gene_id"]] = row
fn = {x["gene"]: x for x in json.load(open('../../next/current_fn_fp.json'))["fn"]}
fp = {x["gene"] for x in json.load(open('../../next/current_fn_fp.json'))["fp"]}
r794 = m.reactions.R794; r795 = m.reactions.R795
print("GPR identical:", r794.gene_reaction_rule == r795.gene_reaction_rule)
br1 = "YALI1A11258g YALI1E14125g YALI1D00581g YALI1E10492g YALI1A09766g YALI1F20965g YALI0E16192g YALI1F31854g YALI1F21690g YALI1E37063g YALI1E32332g YALI1E12482g".split()
br2 = "YALI1F38820g YALI1F13017g".split()
allg = sorted({g.id for g in r794.genes} | {g.id for g in r795.genes})
gids = {g.id for g in m.genes}
for gid in allg:
    g = m.genes.get_by_id(gid)
    off = {gid}
    k794 = not r794.gpr.eval(off); k795 = not r795.gpr.eval(off)
    others = sorted(r.id for r in g.reactions)
    br = "1 (12-gene AND)" if gid in br1 else ("2 (2-gene AND)" if gid in br2 else "?")
    print(f"{gid}\tbranch={br}\tess={'Y' if gid in ess else 'N'}\tKOclose794={k794}\tKOclose795={k795}\trxns={others}\tFN={fn.get(gid,{}).get('cause')}/{fn.get(gid,{}).get('sub')} ko={fn.get(gid,{}).get('ko')}\tFP={gid in fp}")
    print("   model name:", g.name, "| annot:", dict(g.annotation), "| notes:", g.notes)
    if gid in ess: print("   consensus fn:", ess[gid]["function"][:200], "|", ess[gid]["confidence"])
# both branches knocked
print("eval with branch2 both off:", r794.gpr.eval(set(br2)))
# pairs
import itertools
cl=[(a,b) for a in br1 for b in br2 if not r794.gpr.eval({a,b})]
print("closing pairs count:", len(cl), "of", len(br1)*len(br2))
