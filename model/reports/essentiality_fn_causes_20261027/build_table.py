"""Merge the screen, classification and control outputs into the per-gene FN table (no optimization).

Inputs ($IYALI26_RESEARCH_ROOT/artifacts/tasks/essentiality_fn_causes_20261027/fn_classification/):
fn_classification.json, trace_controls.json, trace_extra.json. Output: fn_causes_a1237821.csv next to this file.

Cause categories follow the task brief:
  a  OR rule / isozyme bypass: the gene KO leaves a reaction open through another gene, and closing the gene's
     reactions outright (reaction-level KO) is lethal (< 15 % of WT, or no solution).
  b  alternative or duplicate reaction: the gene's reactions are closed (or closing them all is not lethal) and
     another route carries the flux; growth may drop.
  c  product supplied by the medium: the KO turns essential (< 15 %) once the 11 amino-acid + adenine uptakes are
     closed (trace_controls part 1).
  d  product not demanded: every reaction the gene closes is blocked in WT (no flux possible) or carries no flux
     at the WT optimum (FVA at 100 %), and the reaction-level KO is not lethal.
  e  gene in the model with no reaction.
  f  conditional essentiality flag (added to c): the call depends on an assumed uptake cap.
  g  other.
The rule order is: e, then c (control), then a, then b/d. Sub-labels are mechanical; the REPORT explains them.
"""
import csv, json, os, re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
task = Path(os.environ["IYALI26_RESEARCH_ROOT"]) / "artifacts/tasks" / HERE.name / "fn_classification"
C = json.loads((task / "fn_classification.json").read_text())
T = json.loads((task / "trace_controls.json").read_text())
X = json.loads((task / "trace_extra.json").read_text())

VERIFIED_BYPASS = {  # bypass-closed controls run in trace_extra.py (gene KO + named reaction closed)
    "YALI1B01136g": "R2024 (propanoate:CoA ligase, AMP-forming, reversible) run backwards regenerates ATP from AMP + PPi; with R2024 forward-only the KO grows 0 (WT unchanged)",
    "YALI1F18105g": "R2024 (propanoate:CoA ligase, AMP-forming, reversible) run backwards consumes PPi; with R2024 forward-only the KO grows 0 (WT unchanged)",
    "YALI1F38116g": "leucine catabolism R1968-R1973 makes mitochondrial HMG-CoA, exported by gene-less R1172; with R1172 closed the KO grows 0 (WT 0.99695 of the original)",
    "YALI1C08262g": "pentose phosphate pathway + methylglyoxal route (gene-less R1252); with R1252 closed the KO still grows 0.693 of WT, so R1252 is part of the bypass only",
    "YALI1B11532g": "acetoacetyl-CoA from leucine catabolism R1968-R1973 via mitochondrial HMG-CoA, exported by gene-less R1172; with R1172 closed the KO grows 0",
    "YALI1F00821g": "duplicate reaction R2081 (same equation as R334) keeps glutamine synthesis open through YALI1D16151g; reaction-level KO of R334+R2081 grows 0 (also reports/r2081_closed_f00821_ko_20260911)",
}
GROUP_NOTES = {
    "R171": "R171 ATP synthase rule is (14-subunit complex) OR YALI1B28025g (ATP20) OR YALI1B15733g (ATP19): each single accessory subunit keeps R171 open",
    "R759": "cytosolic R759 (and mitochondrial R760) both carry THS1-like YALI1F14882g OR MST1-like YALI1A09121g; R760 is blocked in WT",
    "R723": "cytosolic R723 carries SES1-like YALI1F03908g OR DIA4-like YALI1B07957g (DIA4 is the mitochondrial SerRS in S. cerevisiae, recalled)",
    "R673": "cytosolic R673 carries cytosolic-ProRS-like YALI1D19110g OR YHR020W-like YALI1E05929g",
    "R2193": "R2193 iron import is (YALI1D09277g and YALI1D09308g) or (YALI1D08684g and YALI1D08564g) (STATE 2026-10-22 open problem 7)",
    "R2062": "complex I closed; cytosolic NADH is oxidised by NDH2-like R570 (YALI1F32476g) and matrix NADH by shuttles (reports/complexI_shuttle_assumption_screen_20261014: closing the R533 or R831 outlet gives KO/WT 0.069, recalled)",
}
COFACTOR = re.compile(r"FOL[0-9]|FMN1|FAD1|QNS1|BNA[0-9]|BPL1|riboflavin|pantothen|coenzyme A|CAB[0-9]|pyridox|biotin|NAD synth|nicotin|folate|ECM31|RIB[0-9]|thiamin|SNZ[0-9]|Q03148|Sno-type", re.I)
MITO_TRANSL = re.compile(r"MSK1|MSM1|MSD1|MST1|MSE1|MSF1|MSR1|MSY1|MSW1|ISM1|SLM5|MSS116|mitochondrial .*tRNA synthetase|DIA4|PET112", re.I)
SPHINGO = re.compile(r"LCB1|LCB2|TSC10|TSC13|SUR2|LAG1|LAC1|AUR1|ceramide|sphing", re.I)
GPI = re.compile(r"GPI[0-9]+|GAA1|GAB1|PIG-|MCD4|SPT14|glycosylphosphatidylinositol", re.I)
GLYCO = re.compile(r"ALG[0-9]+|SEC59|DPM1|PMT[0-9]|glycosyltransferase|mannosyltransferase", re.I)
VATP = re.compile(r"VMA[0-9]+|STV1|VPH1|TFP1|PPA1 H\+-ATPase|P32842|V-ATPase|vacuolar.*ATP", re.I)


def sublabel(fn):
    for lab, rx in (("cofactor_synthesis", COFACTOR), ("mito_translation", MITO_TRANSL), ("sphingolipid", SPHINGO),
                    ("gpi_anchor", GPI), ("protein_glycosylation", GLYCO), ("v_atpase", VATP)):
        if rx.search(fn):
            return lab
    return "other"


def short(fn):
    m = re.search(r"Saccharomyces cerevisiae (?:Y[A-Z]{2}\d+[WC](?:-[A-Z])? )?([A-Z0-9]{3,8}) ?\(?(Y[A-Z]{2}\d+[WC])?", fn)
    sym = m.group(1) if m else None
    return (f"homolog of S. cerevisiae {sym}; " if sym and not sym.startswith("Y") else "") + fn.strip()[:120]


rows, counts = [], Counter()
for g in C["genes"]:
    gid = g["gene_id"]
    p1 = T["part1"][gid]
    aa_ess = p1["aa_closed"]["essential_15"]
    rk = g["reaction_ko_ratio"]
    rk_lethal = g["reaction_ko_status"] != "optimal" or (rk is not None and rk < 0.15)
    info = g["reaction_info"]
    closed = g["closed_by_gene_ko"]
    partners = sorted({x for p in g["or_partners"].values() for x in p["other_genes"]})
    rxset = set(g["reactions"])
    notes = [v for k, v in GROUP_NOTES.items() if k in rxset]
    flags = []
    if not g["reactions"]:
        cat, sub = "e", "gene_without_reaction"
    elif aa_ess:
        cat, sub = "c", "amino_acid_uptake_cap"
        flags.append("f: call depends on the concentration-scaled amino-acid caps (benchmark.md: assumptions, not uptake rates)")
    elif gid in VERIFIED_BYPASS and gid != "YALI1F00821g":
        cat, sub = "b", "bypass_route_verified"
        notes.append(VERIFIED_BYPASS[gid])
    elif gid == "YALI1F00821g":
        cat, sub = "b", "duplicate_reaction_verified"
        notes.append(VERIFIED_BYPASS[gid])
    elif rk_lethal and partners:
        cat, sub = "a", "or_partner_carries_flux"
    elif "R2062" in rxset:
        cat, sub = "b", "alternative_nadh_oxidation"
    else:
        used = any(info[r]["used_at_optimum"] for r in g["reactions"])
        if rk is not None and rk < 0.999:
            cat, sub = "b", "alternative_route_growth_cost"
        elif used:
            cat, sub = "b", "alternative_route_no_growth_cost"
        elif all(info[r]["blocked"] for r in (closed or g["reactions"])):
            cat, sub = "d", "blocked_in_wt:" + sublabel(g["list_function"])
        else:
            cat, sub = "d", "unused_at_optimum:" + sublabel(g["list_function"])
    loops = [r for r in g["reactions"] if max(abs(info[r]["fva_min"]), abs(info[r]["fva_max"])) > 400]
    if loops:
        flags.append("loop: " + ",".join(loops) + " has |FVA| > 400 at the WT optimum (internal cycle)")
    if gid in VERIFIED_BYPASS and VERIFIED_BYPASS[gid] not in notes:
        notes.append(VERIFIED_BYPASS[gid])
    dead = {k: v for k, v in p1["blocked_reaction_dead_ends"].items() if v}
    counts[(cat, sub.split(":")[0])] += 1
    rows.append(dict(
        gene_id=gid, verified_name="none verified", function_short=short(g["list_function"]),
        evidence_level="list annotation by homology; model GPR assignment (not experimental validation)",
        ko_wt_ratio=f"{g['ko_wt_ratio']:.4f}", cause=cat, sub_cause=sub, flags="; ".join(flags),
        reactions=";".join(g["reactions"]), closed_by_gene_ko=";".join(closed), or_partner_genes=";".join(partners),
        reaction_ko_ratio="" if rk is None else f"{rk:.4f}", reaction_ko_status=g["reaction_ko_status"],
        closed_wt_fva=";".join(f"{r}:[{info[r]['fva_min']:.3g},{info[r]['fva_max']:.3g}]{'(blocked)' if info[r]['blocked'] else ''}" for r in closed),
        aa_closed_ratio="" if p1["aa_closed"]["ratio"] is None else f"{p1['aa_closed']['ratio']:.4f}",
        vit_closed_ratio="" if p1["vit_closed"]["ratio"] is None else f"{p1['vit_closed']['ratio']:.4f}",
        dead_end_metabolites=";".join(f"{k}:{'|'.join(v)}" for k, v in dead.items())[:400],
        mechanism_notes=" | ".join(notes)))
fields = list(rows[0])
with (HERE / "fn_causes_a1237821.csv").open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(sorted(rows, key=lambda r: (r["cause"], r["sub_cause"], r["gene_id"])))
print(len(rows))
for k, v in sorted(counts.items()):
    print(k, v)
print(Counter(r["cause"] for r in rows))
print("d sublabels", Counter(r["sub_cause"] for r in rows if r["cause"] == "d"))
print("a by reaction group", Counter(next((k for k in GROUP_NOTES if k in r["reactions"].split(";")), "other") for r in rows if r["cause"] == "a"))
