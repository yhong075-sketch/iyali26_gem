"""Small report tables from the step outputs (no LP). Usage: python make_tables.py ../tables"""
import csv
import json
import os
import sys
import warnings

import cobra

from common import ART, EIGHT, KO_GENES, LIST, MODEL, NEIGH

outdir = sys.argv[1]
os.makedirs(outdir, exist_ok=True)
S1 = json.load(open(ART / "out/step1_capacity.json"))
S4 = json.load(open(ART / "out/step4_options.json"))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    m = cobra.io.read_sbml_model(MODEL)
listed = {row["gene_id"] for row in csv.DictReader(open(LIST))}


def fmt(v):
    return "" if v is None else f"{v:.6g}"


def write(name, rows, cols):
    with open(os.path.join(outdir, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote", name, len(rows))


# 1. the eight reactions and their neighbours: static facts and step-1 capacities
rows = []
for rid in EIGHT + NEIGH:
    r = m.reactions.get_by_id(rid)
    row = {"reaction": rid, "in_task": rid in EIGHT, "name": r.name, "equation": r.build_reaction_string(),
           "gene_rule": r.gene_reaction_rule or "NONE", "bounds": f"({r.lower_bound:g}, {r.upper_bound:g})",
           "compartments": ";".join(sorted(r.compartments)),
           "mass_charge_balance": json.dumps(r.check_mass_balance()) if r.check_mass_balance() else "balanced"}
    for cond in ["open", "taskA"]:
        for solver in ["gurobi", "glpk"]:
            d = S1["A"][f"{cond}/{solver}"].get(rid)
            if d and "min" in d:
                row[f"{cond}_{solver}_min"], row[f"{cond}_{solver}_max"] = fmt(d["min"]), fmt(d["max"])
    if rid in EIGHT:
        for k in ["fva_999999", "fva_999"]:
            d = S1["B"][k][rid]
            row[f"taskA_{k}_min"], row[f"taskA_{k}_max"] = fmt(d["min"]), fmt(d["max"])
        c = S1["C"][rid]
        row["taskA_wt_change_if_closed"] = fmt(c["change"])
        row["closed_wt_status"] = c["status"]
    rows.append(row)
cols = ["reaction", "in_task", "name", "equation", "gene_rule", "bounds", "compartments", "mass_charge_balance",
        "open_gurobi_min", "open_gurobi_max", "open_glpk_min", "open_glpk_max", "taskA_gurobi_min",
        "taskA_gurobi_max", "taskA_glpk_min", "taskA_glpk_max", "taskA_fva_999999_min", "taskA_fva_999999_max",
        "taskA_fva_999_min", "taskA_fva_999_max", "taskA_wt_change_if_closed", "closed_wt_status"]
write("capacity_step1.csv", rows, cols)

# 2. targeted knockouts: as is and every variant that ran knockouts
variants = [("step1", k) for k in S1["D"]["variants"]] + [("step4", k) for k, v in S4["variants"].items() if "ko" in v]
rows = []
for gid in KO_GENES:
    g = m.genes.get_by_id(gid)
    row = {"gene": gid, "on_consensus_list": gid in listed,
           "reactions": ";".join(sorted(r.id for r in g.reactions))}
    for src, key in variants:
        rec = (S1["D"][key]["ko"] if src == "step1" else S4["variants"][key]["ko"])[gid]
        tag = f"{src}:{key}"
        row[f"{tag} ratio"] = fmt(rec["ratio"])
        row[f"{tag} class15"] = rec["class"]["0.15"]
        row[f"{tag} status"] = rec["status"]
    rows.append(row)
cols = ["gene", "on_consensus_list", "reactions"] + [f"{s}:{k} {x}" for s, k in variants
                                                     for x in ("ratio", "class15", "status")]
write("targeted_knockouts.csv", rows, cols)

# 3. options summary
KEY = ["R2025", "R2113", "R2077", "R2083", "R2084", "R2118", "R1291", "R2117", "R2117_mi", "R2114_mi",
       "R2112_mi", "R1969", "R25", "R834", "R2024", "R2285", "T_OBDHm"]
rows = []
for key, v in S4["variants"].items():
    g = S4["glpk"][key]
    row = {"variant": key, "wt_gurobi": fmt(v["wt"]), "wt_status": v["wt_status"], "wt_glpk": fmt(g["wt"]),
           "energy_atp_gurobi": fmt(v["energy"]["R_NTP1"]["value"]),
           "energy_biomass_gurobi": fmt(v["energy"]["biomass_C"]["value"]),
           "energy_atp_glpk": fmt(g["energy"]["R_NTP1"]["value"]),
           "energy_biomass_glpk": fmt(g["energy"]["biomass_C"]["value"]),
           "knockouts_run": "ko" in v}
    if "ko" in v:
        base = S1["D"]["as_is"]["ko"]
        row["ko_class15_changes_vs_as_is"] = sum(
            1 for gid, r in v["ko"].items() if r["class"]["0.15"] != base[gid]["class"]["0.15"])
        row["ko_max_abs_ratio_change"] = fmt(max(abs((r["ratio"] or 0) - (base[gid]["ratio"] or 0))
                                                 for gid, r in v["ko"].items()))
    fva = v.get("fva", S1["A"]["taskA/gurobi"] if key == "B0/K" else {})
    for rid in KEY:
        d = fva.get(rid)
        row[rid] = "" if (not d or "min" not in d) else f"[{d['min']:.4g}, {d['max']:.4g}]"
    rows.append(row)
cols = ["variant", "wt_gurobi", "wt_status", "wt_glpk", "energy_atp_gurobi", "energy_biomass_gurobi",
        "energy_atp_glpk", "energy_biomass_glpk", "knockouts_run", "ko_class15_changes_vs_as_is",
        "ko_max_abs_ratio_change"] + KEY
write("options_summary.csv", rows, cols)

# 4. in-memory edits of the options
rows = []
seen = set()
for key, v in S4["variants"].items():
    e = v["edits"]
    for rid, d in e.get("moved", {}).items():
        if rid not in seen:
            seen.add(rid)
            rows.append({"edit": f"{rid} moved to C_mi as {rid}_mi", "equation": d["equation"], "gene_rule": d["rule"],
                         "H_added": d["H_added"], "balance": json.dumps(d["balance"]) or "{}"})
    if "R2024" in e and "R2024" not in seen:
        seen.add("R2024")
        rows.append({"edit": "m1875 merged into m1970 (R2024 rewritten; R2025 duplicate set to 0)",
                     "equation": e["R2024"]["equation"], "gene_rule": "YALI1F08845g",
                     "H_added": e["R2024_H_added"], "balance": json.dumps(e["R2024"]["balance"])})
    if "T_OBDHm" in e and "T_OBDHm" not in seen:
        seen.add("T_OBDHm")
        rows.append({"edit": "B1 base: T_OBDHm added (inferred PDH rule), R2075 closed",
                     "equation": e["T_OBDHm"]["equation"], "gene_rule": e["T_OBDHm"]["rule"], "H_added": 4,
                     "balance": json.dumps(e["T_OBDHm"]["balance"])})
write("option_edits.csv", rows, ["edit", "equation", "gene_rule", "H_added", "balance"])

# 5. TargetP
tp = {}
for line in open(ART / "targetp/output_protein_type.txt"):
    if line.startswith("#"):
        continue
    p = line.rstrip("\n").split("\t")
    tp[p[0]] = p
REACTIONS = {"YALI1D19252g": "R2117, R2282", "YALI1B14044g": "R2116, R2141, R1978, R1546-R1553",
             "YALI1F03883g": "R2114", "YALI1D07981g": "R2140", "YALI1C02550g": "R2112, R2119, R2121",
             "YALI1E15555g": "R2082", "YALI1F13202g": "R265 (mi), R1092 (cy)",
             "YALI1E31947g": "R693 (PDH; positive control)"}
YALI0 = {"YALI1D19252g": "YALI0D15708g", "YALI1B14044g": "YALI0B10406g", "YALI1F03883g": "YALI0F02607g",
         "YALI1D07981g": "YALI0D06215g", "YALI1C02550g": "YALI0C01859g", "YALI1E15555g": "YALI0E12573g",
         "YALI1F13202g": "YALI0F09625g", "YALI1E31947g": "YALI0E27005g"}
rows = []
for name, p in tp.items():
    gid, acc = name.split("_", 1)
    rows.append({"gene": gid, "yali0_tag (W29 audit)": YALI0[gid], "refseq": acc, "model_reactions": REACTIONS[gid],
                 "on_consensus_list": gid in listed, "prediction": p[1], "other": p[2], "sp": p[3], "mtp": p[4],
                 "cleavage_site": p[5] if len(p) > 5 else "",
                 "job": "6AFBA2A60032297C72156452 (Non-Plant, stamp 20261117014111)"})
for line in open(ART / "targetp/output_b14044_m28_protein_type.txt"):
    if not line.startswith("#"):
        p = line.rstrip("\n").split("\t")
        rows.append({"gene": "YALI1B14044g", "yali0_tag (W29 audit)": "YALI0B10406g",
                     "refseq": "XP_500719.3 residues 28-346 (319 aa, the CLIB122/S2 length)",
                     "model_reactions": REACTIONS["YALI1B14044g"], "on_consensus_list": "YALI1B14044g" in listed,
                     "prediction": p[1], "other": p[2], "sp": p[3], "mtp": p[4], "cleavage_site": p[5],
                     "job": "6AFBA4A800322DDEA260BFEA (Non-Plant, stamp 20261117014945)"})
rows.append({"gene": "YALI1D24649g", "yali0_tag (W29 audit)": "YALI0D19426g", "refseq": "XP_503031.2",
             "model_reactions": "R615 (cy PPAT)", "on_consensus_list": "YALI1D24649g" in listed, "prediction": "OTHER",
             "other": "0.632686", "sp": "0.364989", "mtp": "0.002325", "cleavage_site": "",
             "job": "reused: 6AD665190020334089046654 (orphan_gene_search_20261019, stamp 20261019204442)"})
write("targetp.csv", rows, ["gene", "yali0_tag (W29 audit)", "refseq", "model_reactions", "on_consensus_list",
                            "prediction", "other", "sp", "mtp", "cleavage_site", "job"])
