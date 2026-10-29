"""Vacuolar transports in external models (read-only): match by compound name + compartment, print bounds/GPR."""
import sys, re, warnings
from pathlib import Path
import cobra
warnings.simplefilter("ignore")
KEY = re.compile(r"^(H\+|H2O|water|L-asparagine|L-aspartate|L-glutamate|L-glutamine|L-isoleucine|L-leucine|L-tyrosine|putrescine|spermidine|spermine|glycine|L-glycine|L-alanine|L-proline|ATP|.*Gly.*)", re.I)
for p in sys.argv[1:]:
    m = cobra.io.read_sbml_model(Path(p))
    vac = [c for c, n in m.compartments.items() if n and "vacuol" in n.lower() and "membrane" not in n.lower()]
    print(f"\n=== {Path(p).name}: {len(m.reactions)} rxns; vacuole compartments {vac} ({[m.compartments[c] for c in vac]})")
    for r in m.reactions:
        comps = {x.compartment for x in r.metabolites}
        if not (comps & set(vac)) or len(comps) < 2:
            continue
        vmets = [x for x in r.metabolites if x.compartment in vac]
        if not any(KEY.match(x.name or "") for x in vmets):
            continue
        print(f"{r.id} [{r.lower_bound},{r.upper_bound}] {r.name} | GPR {r.gene_reaction_rule[:120]}")
        print("     ", r.build_reaction_string(use_metabolite_names=True)[:220], "|", sorted(comps))
