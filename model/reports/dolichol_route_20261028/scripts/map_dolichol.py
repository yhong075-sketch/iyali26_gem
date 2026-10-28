"""Read-only map of the dolichol-route fragments in a model (no solves)."""
import sys, re, json, warnings
import cobra
warnings.simplefilter("ignore")
path = sys.argv[1]
M = cobra.io.read_sbml_model(path)
def rxn_row(r):
    gpr = r.gene_reaction_rule
    notes = {k: v for k, v in r.notes.items()}
    ann = {k: v for k, v in r.annotation.items() if k != 'sbo'}
    return f"{r.id:8s} {r.name[:55]:55s} {str(r.bounds):16s} GPR={gpr or '-'}\n           {r.reaction}\n           notes={notes} ann={ann}"
ids = [f"R{i}" for i in range(217, 236)] + [f"R{i}" for i in range(2177, 2183)] + ["R287"]
print("=== named reactions ===")
for rid in ids:
    if rid in M.reactions:
        print(rxn_row(M.reactions.get_by_id(rid)))
print("\n=== metabolites matching dolich|prenyl|farnesyl|isopentenyl|dimethylallyl|geranyl ===")
pat = re.compile(r"dolich|prenyl|farnesyl|isopentenyl|dimethylallyl|geranyl", re.I)
hits = [m for m in M.metabolites if pat.search(m.name) or pat.search(m.id)]
for m in sorted(hits, key=lambda m: (m.name, m.id)):
    ann = {k: v for k, v in m.annotation.items() if k != 'sbo'}
    print(f"{m.id:14s} {m.name[:60]:60s} formula={m.formula} charge={m.charge} rxns={sorted(r.id for r in m.reactions)} ann={ann}")
print("\n=== reactions touching those metabolites (not already listed) ===")
seen = set(ids)
for m in hits:
    for r in m.reactions:
        if r.id not in seen:
            seen.add(r.id)
            print(rxn_row(r))
print("\n=== genes on the fragments ===")
genes = set()
for rid in ids:
    if rid in M.reactions:
        genes |= {g.id for g in M.reactions.get_by_id(rid).genes}
for r in seen:
    genes |= {g.id for g in M.reactions.get_by_id(r).genes}
for g in sorted(genes):
    G = M.genes.get_by_id(g)
    print(g, G.name, sorted(r.id for r in G.reactions))
print("\nbiomass_C dolichol/prenyl terms:", [(m.id, c) for m, c in M.reactions.get_by_id("biomass_C").metabolites.items() if pat.search(m.name)])
