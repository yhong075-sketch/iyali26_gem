"""Per-compartment producers / consumers / transporters of the UM7 metabolites (read-only, no solves).

Copies of one compound are matched by formula plus a shared KEGG or ChEBI identifier (not by name).
A reaction is listed as producer / consumer by the sign of the coefficient and its bounds
(reversible reactions are listed as both); a reaction whose species span two compartments is a transporter.
"""
import sys, warnings
import cobra
warnings.simplefilter("ignore")
M = cobra.io.read_sbml_model(sys.argv[1])
# seed species (one copy each); their copies in every compartment are found below
SEEDS = {"homocitrate": "m628[C_nu]", "diphosphate": "m210[C_nu]", "ADP-ribose": "m818[C_nu]", "acetate": "m205[C_nu]",
         "NAD": "m817[C_nu]", "NADP": "m41[C_cy]", "nicotinamide": "m819[C_nu]", "acetyl-CoA": "m208[C_nu]",
         "NMN": "m832[C_nu]", "deamido-NAD": "m821[C_nu]"}
COMPS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["C_nu", "C_cy", "C_mi"]


def keys(m):
    out = set()
    for k in ("kegg.compound", "chebi", "metanetx.chemical"):
        v = m.annotation.get(k)
        if v:
            out |= {f"{k}:{x}" for x in ([v] if isinstance(v, str) else v)}
    return out


def copies(seed):
    s = M.metabolites.get_by_id(seed)
    ks = keys(s)
    hits = [m for m in M.metabolites if m.formula == s.formula and (keys(m) & ks or m.id == seed)]
    return sorted(hits, key=lambda m: m.compartment), ks


def role(r, m):
    c = r.metabolites[m]
    rev = r.lower_bound < 0 and r.upper_bound > 0
    if r.bounds == (0, 0):
        return "CLOSED"
    if rev:
        return "both"
    fwd = r.upper_bound > 0
    return "produces" if (c > 0) == fwd else "consumes"


for label, seed in SEEDS.items():
    if seed not in M.metabolites:
        print(f"\n##### {label}: seed {seed} not in model"); continue
    cps, ks = copies(seed)
    s = M.metabolites.get_by_id(seed)
    print(f"\n##### {label}  (formula {s.formula}; match keys {sorted(ks)[:4]})")
    print("   copies:", ", ".join(f"{m.id}" for m in cps))
    for m in cps:
        if m.compartment not in COMPS:
            continue
        print(f"  -- {m.id} {m.name} q={m.charge}")
        for r in sorted(m.reactions, key=lambda r: r.id):
            comps = sorted({x.compartment for x in r.metabolites})
            kind = "TRANSPORT" if len(comps) > 1 else "reaction"
            print(f"     {r.id:9s} {role(r, m):8s} {kind:9s} {str(r.bounds):16s} GPR={r.gene_reaction_rule or '-':40s} "
                  f"{r.name[:40]}\n                {r.build_reaction_string()}")
