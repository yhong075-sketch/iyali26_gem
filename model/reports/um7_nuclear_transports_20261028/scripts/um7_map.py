"""Read-only map of the UM7 reactions and their metabolites, per compartment (no solves)."""
import sys, re, warnings
import cobra
warnings.simplefilter("ignore")
M = cobra.io.read_sbml_model(sys.argv[1])
UM7 = ["R92", "R396", "R567", "R568", "R569", "R579", "R580", "R965", "R1000", "R1001", "R1017", "R1055",
       "R1082", "R1212", "R1214", "R1267", "R1270"]


def row(r):
    ann = {k: v for k, v in r.annotation.items() if k != "sbo"}
    return (f"{r.id:7s} {r.name[:50]:50s} {str(r.bounds):18s} subsys={r.subsystem!r}\n"
            f"        {r.build_reaction_string(use_metabolite_names=False)}\n"
            f"        {r.build_reaction_string(use_metabolite_names=True)}\n"
            f"        GPR={r.gene_reaction_rule or '-'} ann={ann}\n        notes={dict(r.notes)}")


print("=== UM7 reactions ===")
for rid in UM7:
    print(row(M.reactions.get_by_id(rid)))
genes = sorted({g.id for rid in UM7 for g in M.reactions.get_by_id(rid).genes})
print("\n=== genes on UM7 ===")
for g in genes:
    G = M.genes.get_by_id(g)
    print(g, repr(G.name), dict(G.annotation), sorted(r.id for r in G.reactions))

# Base ids of the metabolites to trace (looked up by name+formula, then all compartments of the same base id)
print("\n=== species of interest: every compartment copy (by base id) ===")
mets = {m for rid in UM7 for m in M.reactions.get_by_id(rid).metabolites}
bases = sorted({m.id.split("[")[0] for m in mets}, key=lambda s: int(s[1:]))
for b in bases:
    for m in sorted(M.metabolites.query(lambda x: x.id.split("[")[0] == b), key=lambda m: m.compartment):
        print(f"{m.id:14s} {m.name[:45]:45s} {m.formula} q={m.charge} n_rxn={len(m.reactions)}")
