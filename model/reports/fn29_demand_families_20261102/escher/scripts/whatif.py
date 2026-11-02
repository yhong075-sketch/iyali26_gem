"""In-memory what-if variants of the five demand families, taken from the scouts' own builders (nothing is written to
any model file; none of these variants was implemented).

ARCHIVE = $IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/fn29 (the gap session's scratchpad,
archived 2026-11-02). Two builders are imported from the archive unchanged; three scout scripts run their LPs on import
and hard-code scratchpad paths, so the variant lines are copied verbatim below with their source and SHA-256.

- sph_V2   sphingolipid/05_variants.py apply(m, 'V2_ipc_golgi_route', 0.000417)         (316c72dc..., copied)
- gpi_b1   gpi/s2_variants.py apply(model, 'b1') with COEF from Yin 2007                (884c3f8d..., copied)
- cl_a     cardiolipin/variants.py apply(model, {'kind': 'biomass', 'coef': c})          (96a664f8..., imported)
- ng_C2    nglycan/variants2.py apply_C2(m, 0.01)                                        (adabb9eb..., imported)
- ng_B     nglycan/variants2.py apply_B(m)                                               (adabb9eb..., imported)
- mt_rec   mito_translation/run_20261102/variants2.py variant_a2(m, 0.01863252) + met_fix(m) (1a287be4..., copied),
           plus R1218 reversible as in mito_translation/verify_lp/verify2.py variant_a2(r1218=True) (962df2d6...)
"""
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path

import cobra
from cobra import Reaction

ARCHIVE = Path(os.environ["IYALI26_RESEARCH_ROOT"]) / "artifacts/tasks/fn29_demand_families_20261102/fn29"
SOURCES = {
    "sphingolipid/05_variants.py": "316c72dcb72644ec380a4120466ca04d2fbf2b3d376d5a047803bbaee8697b4a",
    "gpi/s2_variants.py": "884c3f8d340c3adb1b3e7ed329cd1ce5243e6760adeb7f4d1f6a18188f38ee9c",
    "cardiolipin/variants.py": "96a664f8c35b567bd4658b2e82d5e99d73f47b882973748206f51730fe8b1469",
    "nglycan/variants2.py": "adabb9ebd342f4b07245f3d52148509ad328680b1421866be39e1ce2e9ca3d8c",
    "mito_translation/run_20261102/variants2.py": "1a287be459f56526d2c378c87472882ecc763ac2f53dcc44e1d76d9ab59f7996",
    "mito_translation/run_20261102/mt_composition.json": "28021f8373a716fccb6894688d7ba2b45379164eb504c548a16852ecd3f4a2b9",
    "mito_translation/verify_lp/verify2.py": "962df2d699d653a5b8eac0041159b2b7b1e777b1813b4fcef0efb575a333d96f",
}
for rel, want in SOURCES.items():
    got = hashlib.sha256((ARCHIVE / rel).read_bytes()).hexdigest()
    if got != want:
        raise SystemExit(f"STOP: {rel} SHA-256 {got} differs from {want}")


def _import(rel, name):
    spec = importlib.util.spec_from_file_location(name, ARCHIVE / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CL = _import("cardiolipin/variants.py", "fn29_cardiolipin_variants")
NG = _import("nglycan/variants2.py", "fn29_nglycan_variants2")

# ---- sphingolipid V2 (05_variants.py lines 11-37, V2 branch, C_SRC = 0.000417)
C_SRC = 0.000417
IPC_LINK = {'m1707[C_go]': 'm1265[C_go]', 'm1708[C_go]': 'm1270[C_go]', 'm1709[C_go]': 'm1266[C_go]', 'm1710[C_go]': 'm1271[C_go]',
            'm1711[C_go]': 'm1263[C_go]', 'm1712[C_go]': 'm1268[C_go]', 'm1713[C_go]': 'm1264[C_go]', 'm1714[C_go]': 'm1269[C_go]',
            'm1715[C_go]': 'm1267[C_go]', 'm1716[C_go]': 'm1272[C_go]'}


def add_rxn(m, rid, mets, lb, ub, name):
    r = Reaction(rid); r.name = name; m.add_reactions([r])
    r.add_metabolites({m.metabolites.get_by_id(k): v for k, v in mets.items()}); r.bounds = (lb, ub)
    return r


def sph_V2(m):
    bm = m.reactions.get_by_id('biomass_C')
    for i in range(1825, 1835):
        m.reactions.get_by_id('R%d' % i).bounds = (0, 1000)
    for a, b in IPC_LINK.items():
        add_rxn(m, 'L_%s_%s' % (a.split('[')[0], b.split('[')[0]), {a: -1, b: 1}, 0, 1000, 'duplicate IPC species link')
    add_rxn(m, 'L_DAG_gm_em', {'m1704[C_gm]': -1, 'm1649[C_em]': 1}, 0, 1000, 'DAG Golgi membrane -> endosomal membrane')
    bm.add_metabolites({m.metabolites.get_by_id('m1271[C_go]'): -C_SRC})
    return {"biomass_term": {"m1271[C_go]": -C_SRC}}


# ---- GPI b1 (s2_variants.py lines 10-13 and 34-41; exit R_DM_GPI left open)
N_CWP_PER_CELL = 1.9e6
PG_PER_CELL = 18.6
NA = 6.02214076e23
GPI_COEF = N_CWP_PER_CELL / NA / (PG_PER_CELL * 1e-12) * 1e3   # mmol GPI-anchored protein per gDW


def gpi_b1(m):
    m.reactions.get_by_id("biomass_C").add_metabolites({m.metabolites.get_by_id("m603[C_er]"): -GPI_COEF})
    return {"biomass_term": {"m603[C_er]": -GPI_COEF}}


# ---- cardiolipin variant a (11_variants.log spec {"kind": "biomass", "coef": 0.0007161586749409294}
#      = 0.08 x 12.98 mg membrane lipid/gDW / 1449.958 g/mol, Sekova 2019)
CL_COEF = 0.08 * 12.98 / 1449.958
assert abs(CL_COEF / 0.0007161586749409294 - 1) < 1e-12, CL_COEF


def cl_a(m):
    CL.apply(m, {"kind": "biomass", "coef": CL_COEF})
    return {"biomass_term": {CL.CL: -CL_COEF}}


# ---- N-glycan C2 (A2 + B, c_n = 0.01 placeholder) and B alone
def ng_C2(m):
    NG.apply_C2(m, 0.01)
    return {"biomass_term": {"m_NGLYC[C_er]": -0.01}}


def ng_B(m):
    NG.apply_B(m)
    return {"biomass_term": {}}


# ---- mitochondrial translation, recommended = a2 + R1251 reversible + R1218 reversible
#      (run_20261102/variants2.py lines 3-23 and 43-53; R1218 as verify_lp/verify2.py line 183)
FRAC = json.load(open(ARCHIVE / "mito_translation/run_20261102/mt_composition.json"))["fraction"]
THREE = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','Q':'Gln','E':'Glu','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'}
MT_COEF = 0.01863252  # = 3.726504 x 0.005 (placeholder, no source)


def trna_pairs(m):
    """Find charged / free mitochondrial tRNA species by name (not by copied id)."""
    pairs = {}
    for x in m.metabolites:
        if x.compartment != 'C_mi' or 'trna' not in x.name.lower(): continue
        mm = re.match(r'^([A-Za-z]{3,4})-tRNA\(([A-Za-z]{3,4})\)', x.name)
        if mm and mm.group(1) != 'fMet': pairs.setdefault(mm.group(1), {})['charged'] = x.id
        mm2 = re.match(r'^tRNA\(([A-Za-z]{3})\), mitchondrial', x.name)
        if mm2: pairs.setdefault(mm2.group(1), {})['free'] = x.id
        if x.name.startswith('fMet-tRNA'): pairs.setdefault('fMet', {})['charged'] = x.id
    return pairs


def add_species(m):
    p = cobra.Metabolite('mtprot[C_mi]', name='mitochondrially synthesised protein (pseudo, residue-equivalent)', compartment='C_mi'); m.add_metabolites([p]); return p


def met_fix(m):
    """Mitochondrial methionine import: R1251 (geneless, mi -> cy only) opened to reversible. A direction decision, labelled."""
    m.reactions.R1251.bounds = (-1000, 1000)


def variant_a2(m, coef):
    """(a2) 14 charged mito tRNAs only, fractions renormalised to the 14 represented residues; no fMet, no free amino acids."""
    pairs = trna_pairs(m); p = add_species(m); st = {}
    have = [a for a in FRAC if THREE[a] in pairs and 'charged' in pairs[THREE[a]]]
    tot = sum(FRAC[a] for a in have)
    for a in have:
        t = THREE[a]; f = FRAC[a] / tot
        st[m.metabolites.get_by_id(pairs[t]['charged'])] = -f; st[m.metabolites.get_by_id(pairs[t]['free'])] = f
    st[p] = 1.0
    rx = cobra.Reaction('MITO_TRANSLATION_14', name='mitochondrial translation pseudo-reaction (14 tRNAs, renormalised)', lower_bound=0, upper_bound=1000)
    rx.add_metabolites(st); m.add_reactions([rx]); m.reactions.biomass_C.add_metabolites({p: -coef}); return rx


def mt_rec(m):
    variant_a2(m, MT_COEF)
    met_fix(m)
    m.reactions.R1218.bounds = (-1000, 1000)
    return {"biomass_term": {"mtprot[C_mi]": -MT_COEF}}


VARIANTS = {"sph_V2": sph_V2, "gpi_b1": gpi_b1, "cl_a": cl_a, "ng_C2": ng_C2, "ng_B": ng_B, "mt_rec": mt_rec}
