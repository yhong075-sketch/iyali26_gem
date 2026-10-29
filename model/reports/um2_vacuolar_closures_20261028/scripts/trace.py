"""Per-compartment producers/consumers/transporters of the UM2 compounds (read-only). Species matched by
the compound part of the name (before '_') and formula-free stem, then listed per compartment; reactions are
classified by stoichiometry (transport = the compound appears in >1 compartment in the reaction)."""
import sys, warnings, re
from collections import defaultdict
from pathlib import Path
import cobra
warnings.simplefilter("ignore")
m = cobra.io.read_sbml_model(Path(sys.argv[1]))
STEMS = ["H+", "H2O", "L-asparagine", "L-aspartate", "L-glutamate", "L-glutamine", "L-isoleucine", "L-leucine",
         "L-tyrosine", "putrescine", "spermidine", "spermine", "L-glycine", "L-alanine", "L-proline", "Gly-Asp",
         "Gly_Glu", "Ala-Gly", "Gly-Pro"]
def stem(x): return (x.name or "").split("_")[0] if not (x.name or "").startswith("Gly_Glu") else "Gly_Glu"
by = defaultdict(list)
for x in m.metabolites:
    by[stem(x)].append(x)
def side(r, x):
    c = r.metabolites[x]
    if r.lower_bound >= 0: return "produces" if c > 0 else "consumes"
    if r.upper_bound <= 0: return "produces" if c < 0 else "consumes"
    return "both"
for s in STEMS:
    mets = sorted(by.get(s, []), key=lambda x: x.compartment)
    print(f"\n### {s}: species {[(x.id) for x in mets]}")
    for x in mets:
        if s in ("H+", "H2O") and x.compartment not in ("C_va", "C_vm"):
            # only summarise the big pools
            rs = x.reactions
            nopen = sum(1 for r in rs if r.bounds != (0, 0))
            print(f"  {x.id}: {len(rs)} reactions ({nopen} open)")
            continue
        rows = []
        for r in sorted(x.reactions, key=lambda r: r.id):
            comps = {y.compartment for y in r.metabolites if stem(y) == s}
            kind = "transport" if len(comps) > 1 else ("exchange" if r.boundary else "reaction")
            rows.append(f"    {r.id:8s} {'CLOSED' if r.bounds == (0, 0) else 'open  '} [{r.lower_bound:g},{r.upper_bound:g}] {side(r, x):8s} {kind:9s} {r.name[:45]:45s} GPR={r.gene_reaction_rule[:60]}")
        print(f"  {x.id} ({x.compartment}): {len(rows)} reactions, open producers "
              f"{sum(1 for r in x.reactions if r.bounds != (0, 0) and side(r, x) in ('produces', 'both'))}, open consumers "
              f"{sum(1 for r in x.reactions if r.bounds != (0, 0) and side(r, x) in ('consumes', 'both'))}")
        if x.compartment in ("C_va", "C_vm", "C_ex") or len(rows) <= 40:
            print("\n".join(rows) if x.compartment in ("C_va", "C_ex") else "\n".join(l for l in rows if "transport" in l))
