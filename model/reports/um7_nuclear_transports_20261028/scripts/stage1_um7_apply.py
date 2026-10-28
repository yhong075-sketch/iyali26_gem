"""apply() of stage1_um7.py, importable without running the LP script (copied verbatim)."""
from cobra import Metabolite, Reaction
# test reactions (gene-less unless stated); ids prefixed T_ so they cannot collide with model ids
LINKS = {
    "T_ADPRtn": ("ADP-ribose nucleus<->cytosol (pore diffusion)", {"m818[C_nu]": -1, "m1725[C_cy]": 1}),
    "T_PPItn": ("diphosphate nucleus<->cytosol (pore diffusion)", {"m210[C_nu]": -1, "m203[C_cy]": 1}),
    "T_NADtn": ("NAD nucleus<->cytosol (pore diffusion)", {"m817[C_nu]": -1, "m122[C_cy]": 1}),
}


def add(model, rid, name, mets, bounds=(-1000, 1000)):
    r = Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([r])
    r.add_metabolites({model.metabolites.get_by_id(k): v for k, v in mets.items()})
    return r


def homocitrate_cy(model):
    if "m_T_hcit[C_cy]" not in model.metabolites:
        model.add_metabolites([Metabolite("m_T_hcit[C_cy]", formula="C7H10O7", charge=0, name="homocitrate (test)",
                                          compartment="C_cy")])
    add(model, "T_HCITtn", "homocitrate nucleus<->cytosol (pore diffusion)", {"m628[C_nu]": -1, "m_T_hcit[C_cy]": 1})
    add(model, "T_HCITtm", "homocitrate cytosol<->mitochondrion (carrier unknown)", {"m_T_hcit[C_cy]": -1, "m75[C_mi]": 1})


def apply(model, label):
    parts = label.split("+")
    for p in parts:
        if p == "as_built":
            continue
        elif p in ("adpr", "ppi", "nad"):
            rid = {"adpr": "T_ADPRtn", "ppi": "T_PPItn", "nad": "T_NADtn"}[p]
            add(model, rid, *LINKS[rid])
        elif p == "close396":
            model.reactions.R396.bounds = (0, 0)
        elif p == "hcit":
            homocitrate_cy(model)
        elif p == "open92":
            model.reactions.R92.bounds = (0, 1000)
        elif p == "censusH11":  # census targeted test H11: homocitrate nu<->mi + sinks on PPi[nu], ADP-ribose[nu]
            add(model, "T_HCIT_nm", "homocitrate nu<->mi (census H11)", {"m628[C_nu]": -1, "m75[C_mi]": 1})
            add(model, "T_SK_ppi_nu", "sink PPi[nu]", {"m210[C_nu]": -1})
            add(model, "T_SK_adpr_nu", "sink ADP-ribose[nu]", {"m818[C_nu]": -1})
        elif p == "collapse":  # R567 rewritten onto the cytosolic species; nuclear R567/R569/R579 closed
            add(model, "T_R567c", "sirtuin NAD hydrolysis, cytosolic copy (test)",
                {"m32[C_cy]": -1, "m122[C_cy]": -1, "m10[C_cy]": 2, "m1725[C_cy]": 1, "m829[C_cy]": 1}, (0, 1000))
            for rid in ("R567", "R569", "R579"):
                model.reactions.get_by_id(rid).bounds = (0, 0)
        else:
            raise ValueError(p)


