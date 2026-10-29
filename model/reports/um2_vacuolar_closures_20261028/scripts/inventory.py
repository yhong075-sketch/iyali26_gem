"""UM2 Stage 1 inventory (read-only): every reaction touching C_va / C_vm species, plus the 37 UM2 reactions."""
import sys, warnings, re, json
from pathlib import Path
import cobra
warnings.simplefilter("ignore")
m = cobra.io.read_sbml_model(Path(sys.argv[1]))
UM2 = "R795 R869 R870 R871 R876 R878 R879 R882 R883 R887 R888 R897 R898 R918 R925 R927 R1163 R1363 R2018 R2019 R2020 R2021 R2026 R2027 R2028 R2029 R2030 R2031 R2032 R2033 R2034 R2035 R2036 R2037 R2038 R2039 R2040".split()
def eq(r):
    return r.build_reaction_string(use_metabolite_names=False) + "   ||   " + r.build_reaction_string(use_metabolite_names=True)
def notes(r):
    return {k: (v if len(str(v)) < 400 else str(v)[:400] + "...") for k, v in (r.notes or {}).items()}
print("compartments:", m.compartments)
va = [x for x in m.metabolites if x.compartment in ("C_va", "C_vm")]
print("vacuolar species:", len(va))
rx = sorted({r for x in va for r in x.reactions}, key=lambda r: int(re.sub(r"\D", "", r.id) or 0))
print("reactions touching C_va/C_vm:", len(rx))
for r in rx:
    tag = "UM2" if r.id in UM2 else "   "
    print(f"{tag} {r.id} [{r.lower_bound},{r.upper_bound}] {r.name} | {r.subsystem} | GPR: {r.gene_reaction_rule}")
    print("      ", eq(r))
    n = notes(r)
    if n: print("       notes:", json.dumps(n)[:600])
print("\nUM2 not touching vacuole:")
for rid in UM2:
    r = m.reactions.get_by_id(rid)
    if r not in rx:
        print(f"{r.id} [{r.lower_bound},{r.upper_bound}] {r.name} | {r.subsystem} | {r.gene_reaction_rule}\n      ", eq(r))
        n = notes(r)
        if n: print("       notes:", json.dumps(n)[:600])
print("\nvacuolar species list:")
for x in sorted(va, key=lambda x: x.id):
    print(x.id, x.name, x.formula, x.charge, len(x.reactions), sorted(r.id for r in x.reactions))
