"""Read-only static checks behind the V-ATPase FN options report. No LP is solved."""
import csv, hashlib, json, pathlib, sys
import cobra
from cobra.core.gene import GPR

model_path = pathlib.Path(sys.argv[1])
consensus_path = pathlib.Path(sys.argv[2])
out = {"model_sha256": hashlib.sha256(model_path.read_bytes()).hexdigest(),
       "consensus_sha256": hashlib.sha256(consensus_path.read_bytes()).hexdigest(),
       "cobra_version": cobra.__version__}
m = cobra.io.read_sbml_model(model_path)
out["n_genes"], out["n_reactions"] = len(m.genes), len(m.reactions)
with consensus_path.open() as fh:
    rows = list(csv.DictReader(fh))
key = [k for k in rows[0] if "gene" in k.lower() or "locus" in k.lower()][0]
consensus = {r[key].strip() for r in rows}
out["consensus_key"], out["consensus_n"] = key, len(consensus)

def rxn(rid):
    r = m.reactions.get_by_id(rid)
    return {"id": rid, "name": r.name, "bounds": list(r.bounds),
            "equation": r.build_reaction_string(), "gpr": r.gene_reaction_rule,
            "imbalance": {k: float(v) for k, v in r.check_mass_balance().items()}}

out["R794"], out["R795"] = rxn("R794"), rxn("R795")
assert out["R794"]["gpr"] == out["R795"]["gpr"]
genes = sorted(g.id for g in m.reactions.R794.genes)
out["vatpase_genes"] = {g: {"consensus_essential": g in consensus} for g in genes}

a1, a2 = "YALI1E12482g", "YALI1F38820g"
core = [g for g in genes if g not in (a1, a2)]
corrected = "(" + " and ".join(core) + f") and ({a1} or {a2})"
rules = {"current": GPR.from_string(out["R794"]["gpr"]), "corrected_core_and_a1_or_a2": GPR.from_string(corrected)}
out["corrected_rule"] = corrected
for name, gpr in rules.items():
    out[f"single_ko_closes_{name}"] = {g: (not gpr.eval({g})) for g in genes}

def neighbours(mid):
    x = m.metabolites.get_by_id(mid)
    return [{"id": r.id, "coef": float(r.get_coefficient(mid)), "bounds": list(r.bounds),
             "gpr": r.gene_reaction_rule} for r in x.reactions]

for mid in ["m293[C_go]", "m1007[C_va]", "m676[C_er]", "m674[C_er]", "m1314[C_go]", "m948[C_cy]"]:
    out[f"neighbours_{mid}"] = neighbours(mid)

bio = m.reactions.biomass_C
out["biomass_species_in_go_or_va"] = [x.id for x in bio.metabolites if x.compartment in ("C_go", "C_va")]
out["biomass_m948_coef"] = float(bio.get_coefficient("m948[C_cy]"))
out["boundary_fe_species"] = [{"rxn": r.id, "bounds": list(r.bounds), "species": [x.id for x in r.metabolites]}
                              for r in m.boundary if any(x.formula and "Fe" in x.formula for x in r.metabolites)]
for rid in ["R2193", "R1189", "R2048", "R1240", "R439", "R1159", "R1163", "R_NTP1", "R_NDP7g", "R586"]:
    out[rid] = rxn(rid)
pathlib.Path(sys.argv[3]).write_text(json.dumps(out, indent=1))
