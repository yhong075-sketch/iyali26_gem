"""Assign each gene-less metabolic reaction one recommended next action (first matching rule wins).

Usage: python triage.py orphan_evidence.csv triage.csv
The rules use only model-derived facts in orphan_evidence.csv plus the fixed lists below,
whose basis is given in REPORT.md. A group is a recommendation for review, not a decision.
"""
import sys
import pandas as pd

e = pd.read_csv(sys.argv[1]).fillna("")

COQ9_LINE = {"R19", "R321", "R1043", "R52"}          # CoQ9 approval gate (AGENTS.md section 2)
LIKELY_SPONTANEOUS = {"R1210", "R2197"}               # not flagged by name; needs a cited source
INTERNAL_LOOP = {"R1964", "R2041", "R2042"}           # loop_check.txt: flux with all boundaries closed
DIRECTION_AT_OPTIMUM = {"R1901", "R1308"}             # carry reverse flux in WT pFBA
COMPARTMENT_EXCLUDE = {"R321"}                        # R218 match is cis-prenyl (dolichol), different stereochemistry


def group(r):
    if r.spontaneous_name_flag:
        return "0 excluded: spontaneous by name"
    if r.selection_cleared_gpr:
        return "1 restore reviewed GPR cleared by metadata selection (user decision)"
    if r.id in COQ9_LINE:
        return "2 CoQ9 line, approval gate in force"
    if r.id in LIKELY_SPONTANEOUS:
        return "3 likely spontaneous, cite before reclassifying"
    if r.exact_duplicates:
        return "4 exact duplicate pair, merge then curate one GPR"
    if r.id in INTERNAL_LOOP:
        return "5 on a net-zero internal loop, review direction and pathway"
    if r.id in DIRECTION_AT_OPTIMUM:
        return "6 reverse flux at optimum, review direction"
    if any("[no GPR]" not in x for x in r.compartment_free_duplicates.split(";") if x.strip()) \
            and r.id not in COMPARTMENT_EXCLUDE:
        return "7 compartment copy of a gene-bearing reaction, check targeting"
    if r.missing_formula or not r.mass_balanced:
        return "8 chemistry incomplete, fix formula or balance first"
    if r.flux_status.startswith("blocked"):
        return "9 blocked, run gap analysis before gene work"
    return "10 can carry flux, gene search by BLAST and reciprocal check"


e["recommended_group"] = e.apply(group, axis=1)
flags = []
for r in e.itertuples():
    f = []
    if r.missing_formula:
        f.append("missing formula")
    elif not r.mass_balanced:
        f.append("mass imbalance")
    if not r.charge_balanced:
        f.append("charge imbalance")
    if r.flux_status.startswith("required"):
        f.append(f"required at optimum (KO/WT {r.ko_wt_ratio:.4f})")
    if r.compartment_free_duplicates:
        f.append("compartment-free match: " + r.compartment_free_duplicates)
    flags.append("; ".join(f))
e["flags"] = flags
cols = ["id", "name", "recommended_group", "flags", "flux_status", "ko_wt_ratio", "compartments", "ec",
        "kegg_reaction", "selection_cleared_gpr", "exact_duplicates", "no_other_source", "no_other_sink", "equation"]
e[cols].sort_values(["recommended_group", "id"],
                    key=lambda s: s.str.split(" ").str[0].astype(int) if s.name == "recommended_group" else s
                    ).to_csv(sys.argv[2], index=False)
print(e.recommended_group.value_counts().sort_index(key=lambda s: s.str.split(" ").str[0].astype(int)).to_string())
for g, d in e.groupby("recommended_group"):
    print(g, "->", " ".join(sorted(d.id, key=lambda x: int(x[1:]))))
