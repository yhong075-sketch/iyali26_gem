"""In-memory UM11 option variants (read-only; never writes a model file).

Labels are joined with '+'. Test reaction ids start with T_ so they cannot collide with model ids.
  as_built   no change (= option B: nuclear copies stay blocked, documented)
  trx        T_TRX1tn  reduced thioredoxin  m632[C_cy] <-> m931[C_nu]  (gene-less pore diffusion)
  trxss      T_TRXSStn thioredoxin disulphide m633[C_cy] <-> m933[C_nu] (gene-less pore diffusion)
  h2o1914    R1914 water written in the nucleus (m32[C_cy] -> m626[C_nu]); compartment fix only
  ntrr       nuclear thioredoxin reductase: NADPH/NADP nucleus species + diffusions + R755 copy (same GPR)
  close_nu   the 4 nuclear RNR copies and R782 closed (0, 0) (documentation form of option B)
"""
from cobra import Metabolite, Reaction

NUCLEAR_RNR = ["R707", "R709", "R711", "R1914"]


def add(model, rid, name, mets, bounds=(-1000, 1000), gpr=""):
    r = Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in mets.items()})
    if gpr:
        r.gene_reaction_rule = gpr
    return r


def nuclear_copy(model, cy_id, new_id):
    src = model.metabolites.get_by_id(cy_id)
    if new_id not in model.metabolites:
        model.add_metabolites([Metabolite(new_id, formula=src.formula, charge=src.charge, name=src.name + " (test)",
                                          compartment="C_nu")])


def apply(model, label):
    for p in label.split("+"):
        if p == "as_built":
            continue
        elif p == "trx":
            add(model, "T_TRX1tn", "reduced thioredoxin nucleus<->cytosol (pore diffusion)",
                {"m632[C_cy]": -1, "m931[C_nu]": 1})
        elif p == "trxss":
            add(model, "T_TRXSStn", "thioredoxin disulphide nucleus<->cytosol (pore diffusion)",
                {"m633[C_cy]": -1, "m933[C_nu]": 1})
        elif p == "h2o1914":
            r = model.reactions.R1914
            r.add_metabolites({model.metabolites.get_by_id("m32[C_cy]"): -1,
                               model.metabolites.get_by_id("m626[C_nu]"): 1})
        elif p == "ntrr":
            nuclear_copy(model, "m40[C_cy]", "m_T_nadph[C_nu]")
            nuclear_copy(model, "m41[C_cy]", "m_T_nadp[C_nu]")
            add(model, "T_NADPHtn", "NADPH nucleus<->cytosol (test)", {"m40[C_cy]": -1, "m_T_nadph[C_nu]": 1})
            add(model, "T_NADPtn", "NADP nucleus<->cytosol (test)", {"m41[C_cy]": -1, "m_T_nadp[C_nu]": 1})
            add(model, "T_R755n", "thioredoxin reductase (NADPH), nuclear copy (test)",
                {"m627[C_nu]": -1, "m_T_nadph[C_nu]": -1, "m933[C_nu]": -1, "m_T_nadp[C_nu]": 1, "m931[C_nu]": 1},
                (0, 1000), model.reactions.R755.gene_reaction_rule)
        elif p == "close_nu":
            for rid in NUCLEAR_RNR + ["R782"]:
                model.reactions.get_by_id(rid).bounds = (0, 0)
        else:
            raise ValueError(p)
