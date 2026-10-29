"""UM11 inventory: reactions, GPRs, and every reaction touching the involved species in nu/cy (read-only)."""
import sys, warnings, cobra
warnings.simplefilter("ignore")
m = cobra.io.read_sbml_model(sys.argv[1])
UM11 = ["R707", "R709", "R711", "R1914", "R782", "R1033", "R1076", "R1077", "R1093", "R1136", "R1353"]
def eq(r):
    return r.build_reaction_string(use_metabolite_names=False)
def eqn(r):
    return r.build_reaction_string(use_metabolite_names=True)
mets = set()
for rid in UM11:
    r = m.reactions.get_by_id(rid)
    print(f"{rid} | {r.name} | {r.bounds} | GPR: {r.gene_reaction_rule} | subsys: {r.subsystem} | comps {sorted(r.compartments)}")
    print("    ", eq(r)); print("    ", eqn(r))
    print("     notes:", {k: v for k, v in r.notes.items()} if r.notes else {}, "| annot:", dict(r.annotation))
    mets |= {x.id for x in r.metabolites}
print()
print("=== species in UM11 and their compartment copies ===")
names = {}
for mid in sorted(mets):
    x = m.metabolites.get_by_id(mid)
    names.setdefault(x.name, []).append(mid)
skip = {"H+", "H2O", "ATP", "ADP"}
for name in sorted(names):
    copies = [x for x in m.metabolites if x.name == name]
    for x in copies:
        print(f"--- {x.id} {x.name} formula={x.formula} charge={x.charge} comp={x.compartment} n_rxns={len(x.reactions)}")
        if x.name in ("H+", "H2O") : continue
        if len(x.reactions) > 40:
            print("      (many reactions, skipped)"); continue
        for r in sorted(x.reactions, key=lambda r: r.id):
            c = r.metabolites[x]
            print(f"      {r.id:8s} coef {c:+g} b={r.bounds} gpr=[{r.gene_reaction_rule}] {eq(r)}")
