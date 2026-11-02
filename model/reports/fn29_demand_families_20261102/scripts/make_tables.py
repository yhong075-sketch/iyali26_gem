"""Per-gene and per-family CSV tables for the 29 demand-family false negatives (no LP; reads archived JSON only).

Usage: python make_tables.py OUTDIR
Reads $IYALI26_RESEARCH_ROOT/artifacts/tasks/fn29_demand_families_20261102/fn29/ (fails closed if unset or missing).
Every number is parsed from the JSON the scouts and the LP verifiers wrote on e1e7c868; the text columns (family,
why no effect, verdict) restate the classifications of REPORT.md sections 2-3 and contain no numbers.
The knockout values come from the LP lens (independent reproduction) wherever it ran the variant, otherwise from the
scout (column ko_source says which file).
"""
import csv
import json
import os
import sys
from pathlib import Path

root = os.environ.get("IYALI26_RESEARCH_ROOT")
if not root:
    sys.exit("STOP: IYALI26_RESEARCH_ROOT is unset")
A = Path(root) / "artifacts/tasks/fn29_demand_families_20261102/fn29"
if not A.is_dir():
    sys.exit(f"STOP: archive {A} not found")
out = Path(sys.argv[1])
load = lambda rel: json.load(open(A / rel))
CUT = 0.15

SPH = load("sphingolipid/verify_lp/verify2.json")
SPH_SCOUT = load("sphingolipid/05_variants.json")
SPH_V1 = next(r for r in SPH_SCOUT["runs"] if r["variant"] == "V1_ceramide_er_C26" and abs(r["coef"] - 0.000417) < 1e-12)
SPH_SWEEP = load("sphingolipid/08_compare.json")["summary"]
GPI = load("gpi/verify_lp/repro_gpi_lp.json")
GPI_SCOUT = load("gpi/s2_variants.json")["variants"]
CL = load("cardiolipin/verify_lp/verify_cardiolipin.json")
NG = load("nglycan/verify_lp/verify_c2.json")
MT = load("mito_translation/verify_lp/verify2.json")

ALL_SPH = "no sphingolipid term in biomass_C (legacy R1371 is (0,0)), so the ER chain from R722 carries no flux"
ELONG = "no C20-C26 acyl-CoA is demanded because nothing demands a sphingolipid"
NG_TEXT = {
    "YALI1B21618g": "ER chain R2108-R2105 ends at m1969 with no consumer (UM18); no OST and no N-glycan demand",
    "YALI1E40471g": "R173 sits on cytosolic GlcNAc2-PP-Dol m344 (no producer) making m345 (no consumer); no N-glycan demand",
    "YALI1F28833g": "mannan is met by the PMT lump R289 (which creates carbon), so Golgi GDP-mannose is never needed; "
                    "gene-less R1137 duplicates R850",
}
GENES = [
    # gene, family, why no effect as built, verdict (REPORT.md sections 2-3)
    ("YALI1C04265g", "sphingolipid", ALL_SPH, "solvable after decision"),
    ("YALI1F20637g", "sphingolipid", ALL_SPH, "solvable after decision"),
    ("YALI1B22883g", "sphingolipid", ALL_SPH, "solvable after decision"),
    ("YALI1A05116g", "sphingolipid", ELONG, "solvable after decision"),
    ("YALI1B26350g", "sphingolipid", ELONG, "solvable after decision"),
    ("YALI1F15781g", "sphingolipid", ELONG, "solvable after decision"),
    ("YALI1D08000g", "sphingolipid", "GPR mis-assignment: annotated delta(4)-desaturase on the cytosolic IPC synthase set",
     "not by curation (GPR correction and a GlcCer branch needed)"),
    ("YALI1E19647g", "sphingolipid", "Golgi AUR1 set R1825-R1834 inherited (0,0) and dead-ended (duplicate IPC species, "
     "Golgi-membrane DAG with no exit)", "solvable after decision"),
    ("YALI1E20626g", "sphingolipid", ALL_SPH, "solvable after decision"),
    ("YALI1F17082g", "sphingolipid", ALL_SPH, "solvable after decision"),
] + [(g, "GPI anchor", "chain feasible, but GPI m603[C_er] is drained only by the unforced R_DM_GPI; no GPI term in "
      "biomass_C", "solvable after decision")
     for g in ["YALI1C32810g", "YALI1E25932g", "YALI1A21970g", "YALI1F27343g", "YALI1F15403g", "YALI1D05270g"]] + [
    (g, "cardiolipin", "no cardiolipin or PG term in biomass_C; the linear route R1791-R197-R1836-R1722-R1721 is unused",
     "solvable after decision") for g in ["YALI1E17162g", "YALI1F31283g"]] + [
    ("YALI1B21618g", "N-glycan", NG_TEXT["YALI1B21618g"], "not by curation now (blocked on composition evidence)"),
    ("YALI1E40471g", "N-glycan", NG_TEXT["YALI1E40471g"], "not by curation now (blocked on composition evidence)"),
    ("YALI1F28833g", "N-glycan", NG_TEXT["YALI1F28833g"], "solvable after decision (variant B)"),
] + [(g, "mitochondrial translation", "charged mitochondrial tRNA has no consumer: no translation reaction and no "
      "mitochondrial-protein term", "not by curation now (blocked on coefficient source and transport evidence)")
     for g in ["YALI1C25782g", "YALI1D04566g", "YALI1D30570g", "YALI1E05833g", "YALI1E21230g", "YALI1E34408g",
               "YALI1F33819g", "YALI1F03527g"]]
assert len(GENES) == 29 and len({g for g, *_ in GENES}) == 29


def ko_row(fam, g):
    """(closed reactions, ko_wt baseline, best variant, ko_wt best, source, alt variant, ko_wt alt, alt source)."""
    if fam == "sphingolipid":
        b, v = SPH["baseline"]["genes"][g], SPH["v2"]["genes"][g]
        return (b["closed"], b["ko_wt"], "V2", v["ko_wt"], "sphingolipid/verify_lp/verify2.json",
                "V1 (4.17e-4)", SPH_V1["genes"][g]["ko_wt"], "sphingolipid/05_variants.json (scout)")
    if fam == "GPI anchor":
        b, v = GPI["baseline"]["ko"][g], GPI["b2"]["ko"][g]
        return (b["knocked"], b["ko_wt"], "b2 (= b1)", v["ko_wt"], "gpi/verify_lp/repro_gpi_lp.json",
                "b1", GPI_SCOUT["b1"]["knockouts"][g]["ko_wt"], "gpi/s2_variants.json (scout)")
    if fam == "cardiolipin":
        b, v = CL["baseline"]["ko"][g], CL["variant_a"]["ko"][g]
        return (b["closed_reactions"], b["ko_wt"], "a", v["ko_wt"], "cardiolipin/verify_lp/verify_cardiolipin.json", "", "", "")
    if fam == "N-glycan":
        b, v, alt = NG["baseline"]["ko"][g], NG["C2"]["ko"][g], NG["B"]["ko"][g]
        return (b["closed"], b["ko_wt"], "C2 (= A2 + B)", v["ko_wt"], "nglycan/verify_lp/verify_c2.json",
                "B", alt["ko_wt"], "nglycan/verify_lp/verify_c2.json")
    b, v = MT["ko_baseline"][g], MT["recommended"]["ko"][g]
    return (b["closed"], b["ko_wt"], "a2 + R1251 rev + R1218 rev", v["ko_wt"], "mito_translation/verify_lp/verify2.json", "", "", "")


flip = lambda x: "" if x in ("", None) else ("yes" if x < CUT else "no")
fmt = lambda x: "" if x in ("", None) else repr(float(x))
out.mkdir(parents=True, exist_ok=True)
rows = []
for g, fam, why, verdict in GENES:
    closed, kb, var, kv, src, alt, ka, asrc = ko_row(fam, g)
    rows.append({"gene": g, "family": fam, "reactions_closed_by_ko": " ".join(sorted(closed, key=lambda r: (len(r), r))), "why_no_effect_as_built": why,
                 "best_variant": var, "ko_wt_baseline": fmt(kb), "ko_wt_best_variant": fmt(kv), "flips_at_15pct": flip(kv),
                 "ko_source": src, "alt_variant": alt, "ko_wt_alt_variant": fmt(ka), "alt_flips_at_15pct": flip(ka),
                 "alt_source": asrc, "verdict": verdict})
with open(out / "genes_29.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)


def tally(c):
    return "/".join(str(c[k]) for k in ("TP", "FN", "FP", "TN")) + (f" ({c['unresolved']} unresolved)" if "unresolved" in c else "")


def changes(pairs):
    pairs = list(pairs)
    tp =sorted(g for g, a in pairs if a == "TP")
    fp = sorted(g for g, a in pairs if a == "FP")
    return tp, fp


sph_ch = changes((c["gene"], c["after"]) for c in SPH["class_changes"] if c["cut"] == CUT)
gpi_ch = changes((c["gene"], c["after"]) for c in GPI["b2"]["class_changes"] if c["cutoff"] == CUT)
cl_ch = changes((c["gene"], c["after_class"]["0.15"]) for c in CL["class_changes"])
ng_ch = changes((c["gene"], c["becomes_at_15pct"]) for c in NG["class_changes_own_numbers"])
mt_ch = changes((g, c["classes"]["0.15"][1]) for g, c in MT["class_changes"].items())
fam_genes = lambda fam: {g for g, f, *_ in GENES if f == fam}
FAMS = [
    ("sphingolipid", "V2", SPH["baseline"]["wt_growth"], SPH["v2"]["wt_growth"], SPH["v2"]["genes"],
     tally(SPH_SWEEP["before"]["0.15"]), tally(SPH["v2"]["sweep_counts"]["0.15"]), "rule A (infeasible unresolved)",
     sph_ch, SPH["lp_count"], "sphingolipid/verify_lp/verify2.json; 08_compare.json (before)",
     "solvable after decision (9); YALI1D08000g not by curation"),
    ("GPI anchor", "b2 (= b1)", GPI["baseline"]["wt"], GPI["b2"]["wt"], GPI["b2"]["ko"],
     tally(GPI["b2"]["sweep_classes"]["0.15"]["before"]), tally(GPI["b2"]["sweep_classes"]["0.15"]["after"]),
     "mt-excl (16 mtDNA genes excluded, infeasible unresolved)", gpi_ch, GPI["lp_count"],
     "gpi/verify_lp/repro_gpi_lp.json", "solvable after decision"),
    ("cardiolipin", "a", CL["baseline"]["wt"]["growth"], CL["variant_a"]["wt"]["growth"], CL["variant_a"]["ko"],
     tally(CL["before_hybrid"]["totals"]["0.15"]), tally(CL["after_sweep"]["totals"]["0.15"]),
     "mt-excl (16 mtDNA genes excluded, infeasible unresolved)", cl_ch, CL["lp_count"],
     "cardiolipin/verify_lp/verify_cardiolipin.json", "solvable after decision"),
    ("N-glycan", "C2 (= A2 + B)", NG["baseline"]["wt_growth"], NG["C2"]["wt_growth"], NG["C2"]["ko"],
     tally(NG["baseline_confusion_mixed_provenance"]["0.15"]), tally(NG["C2_sweep"]["confusion"]["0.15"]),
     "infeasible as growth 0, mtDNA + COX9 in scope", ng_ch, NG["lp_count"], "nglycan/verify_lp/verify_c2.json",
     "split: VRG4 solvable after decision (B); ALG7, ALG1 blocked on composition evidence"),
    ("mitochondrial translation", "a2 + R1251 rev + R1218 rev", MT["wt_baseline"]["growth"], MT["recommended"]["wt"]["growth"],
     MT["recommended"]["ko"], tally(MT["sweep_counts"]["0.15"]["scout_baseline"]),
     tally(MT["sweep_counts"]["0.15"]["after_recommended"]), "rule A (infeasible unresolved)", mt_ch, MT["lp_total"],
     "mito_translation/verify_lp/verify2.json", "not by curation now (blocked on coefficient source and transport evidence)"),
]
frows = []
for fam, var, wt0, wt1, ko, before, after, rule, (tp, fp), lps, src, verdict in FAMS:
    members = fam_genes(fam)
    flipped = sorted(g for g in members if ko[g]["ko_wt"] is not None and ko[g]["ko_wt"] < CUT)
    frows.append({"family": fam, "n_genes": len(members), "best_variant": var, "wt_baseline": repr(wt0), "wt_variant": repr(wt1),
                  "wt_change_pct": f"{100 * (wt1 - wt0) / wt0:.4f}", "family_fn_to_tp_15pct": len(flipped),
                  "new_fp_15pct": " ".join(fp), "side_effect_tp_15pct": " ".join(sorted(set(tp) - members)),
                  "sweep_15pct_before_tp_fn_fp_tn": before, "sweep_15pct_after_tp_fn_fp_tn": after, "sweep_rule": rule,
                  "lp_lens_lps_in_json": lps, "source": src, "verdict": verdict})
with open(out / "families_5.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(frows[0]))
    w.writeheader()
    w.writerows(frows)
for r in frows:
    print(r["family"], r["family_fn_to_tp_15pct"], "/", r["n_genes"], "FP:", r["new_fp_15pct"], "side TP:", r["side_effect_tp_15pct"],
          r["sweep_15pct_before_tp_fn_fp_tn"], "->", r["sweep_15pct_after_tp_fn_fp_tn"], "LPs", r["lp_lens_lps_in_json"])
print(sum(r["flips_at_15pct"] == "yes" for r in rows), "of 29 flip in the best variant")
