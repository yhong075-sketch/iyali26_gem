"""Find homocitrate -> homoaconitate -> homoisocitrate reactions in external models by metabolite name (read-only)."""
import cobra, sys, warnings, re
from pathlib import Path
warnings.simplefilter("ignore")
pat = re.compile(r"homocitrat|homoaconit|homoisocitr|but-1-ene-1,2,4|hydroxybutane-1,2,4|methylcitrate dehydratase", re.I)
for p in sys.argv[1:]:
    m = cobra.io.read_sbml_model(Path(p))
    print("#####", p, len(m.reactions))
    hits = {r for r in m.reactions if pat.search(r.name or "") or any(pat.search(x.name or "") for x in r.metabolites)}
    for r in sorted(hits, key=lambda r: r.id):
        if any(k in r.name.lower() for k in ("transport","exchange")): continue
        print(f"{r.id} | {r.name} | {r.bounds} | {r.gene_reaction_rule}")
        print("    ", r.build_reaction_string(use_metabolite_names=True))
