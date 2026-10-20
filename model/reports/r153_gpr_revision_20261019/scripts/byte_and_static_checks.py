"""Read-only byte-level and static checks for the R153 GPR revision.

Usage: python -I byte_and_static_checks.py BASELINE.xml CONTROL_SEED0.xml CANDIDATE.xml
Prints JSON: which reaction blocks differ between the files (and whether the differing lines are
identical once sorted), whether all text outside reaction blocks is identical, and for the candidate
the R153 GPR, element/charge balance, gene count and genes without reactions.
"""
import hashlib
import json
import re
import sys

REACTION = re.compile(r'<reaction [^>]*id="([^"]+)".*?</reaction>', re.S)


def blocks(text):
    return {m.group(1): m.group(0) for m in REACTION.finditer(text)}


def compare(a, b):
    ba, bb = blocks(a), blocks(b)
    differing = [rid for rid in ba if ba[rid] != bb.get(rid)]
    return {"same_reaction_ids_and_order": list(ba) == list(bb),
            "differing_reactions": {rid: {"same_lines_when_sorted": sorted(ba[rid].splitlines())
                                          == sorted(bb.get(rid, "").splitlines())} for rid in differing},
            "non_reaction_text_identical": REACTION.sub("", a) == REACTION.sub("", b)}


def formula(f):
    return {e: int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", f or "")}


def static(text):
    attr = lambda tag, a: (re.search(r'\b%s="([^"]*)"' % a, tag) or [None, None])[1]
    species = {attr(t, "id"): (attr(t, "fbc:chemicalFormula"), attr(t, "fbc:charge"))
               for t in re.findall(r"<species\b[^>]*>", text)}
    r153 = blocks(text)["R_R153"]
    sides = {tag: re.findall(r'species="([^"]+)" stoichiometry="([^"]+)"',
                             re.search(r"<%s>(.*?)</%s>" % (tag, tag), r153, re.S).group(1))
             for tag in ("listOfReactants", "listOfProducts")}
    elements, charge = {}, 0.0
    for tag, sign in (("listOfReactants", -1), ("listOfProducts", 1)):
        for sid, coef in sides[tag]:
            f, z = species[sid]
            for e, n in formula(f).items():
                elements[e] = elements.get(e, 0) + sign * float(coef) * n
            charge += sign * float(coef) * float(z or 0)
    genes = set(re.findall(r'<fbc:geneProduct [^>]*fbc:id="(G_[^"]+)"', text))
    used = set(re.findall(r'geneProductRef fbc:geneProduct="(G_[^"]+)"', text))
    return {"R153_has_gene_association": "geneProductAssociation" in r153,
            "R153_element_imbalance": {e: v for e, v in elements.items() if abs(v) > 1e-9},
            "R153_charge_imbalance": charge,
            "gene_products": len(genes), "genes_without_reactions": sorted(genes - used),
            "G_YALI1D17462g_present": "G_YALI1D17462g" in genes,
            "G_YALI1D17462g_referenced": "G_YALI1D17462g" in used}


def main(baseline, control, candidate):
    texts = {p: open(p, "rb").read() for p in (baseline, control, candidate)}
    a, b, c = (texts[p].decode() for p in (baseline, control, candidate))
    print(json.dumps({
        "sha256": {p.split("/")[-1]: hashlib.sha256(t).hexdigest() for p, t in texts.items()},
        "control_seed0_vs_baseline": compare(a, b),
        "candidate_vs_control_seed0": compare(b, c),
        "baseline_static": static(a),
        "candidate_static": static(c),
    }, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])
