"""Step 3 (no LP): gene rules of the nine reactions as recorded in 7ccbb568, each gene's consensus-list status,
W29 tag-protein audit flag, other reactions in the model, and its recorded certified knockout result.

The knockout results are read from the certified sweep on 353cefdc
($IYALI26_RESEARCH_ROOT/artifacts/tasks/folinic_acid_module_20261112/sweep684_candidate.tsv, Task A + PO1f,
Gurobi 1e-9); 353cefdc and 7ccbb568 have identical bounds, stoichiometry, gene rules and objective (static check,
2026-11-16), so the same LPs apply. Writes ../tables/genes_of_the_nine.csv.
"""
import csv
import warnings
from pathlib import Path

import cobra

from common import MODEL, NINE, RR, REPO, check_inputs, sha256

SWEEP = RR / "artifacts/tasks/folinic_acid_module_20261112/sweep684_candidate.tsv"
SWEEP_SHA = "49b4f27be05beffe32b206b75e7fe571b69d07af6d13a48b568f8e39f07e91db"
LIST = RR / "state/essentiality/repository/consensus_essential_genes.csv"
LIST_SHA = "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b"
AUDIT = REPO / "model/curation/annotation/w29_tag_protein_audit.csv"
OUT = Path(__file__).resolve().parents[1] / "tables" / "genes_of_the_nine.csv"

check_inputs()
assert sha256(SWEEP) == SWEEP_SHA and sha256(LIST) == LIST_SHA
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    model = cobra.io.read_sbml_model(MODEL)
positives = {row["gene_id"]: row for row in csv.DictReader(open(LIST))}
sweep = {row["gene"]: row for row in csv.DictReader(open(SWEEP), delimiter="\t")}
audit = {row["gene"]: row for row in csv.DictReader(open(AUDIT))}
wt = 0.3111005359259494  # WT of 353cefdc / 7ccbb568 on Task A + PO1f (step 1 reproduces 0.311100535925948)

rows = []
for rid in NINE:
    r = model.reactions.get_by_id(rid)
    genes = sorted(g.id for g in r.genes)
    if not genes:
        rows.append([rid, r.name, r.compartments and ";".join(sorted(r.compartments)), r.gene_reaction_rule or "",
                     "", "", "", "", "", "", "", ""])
    for gid in genes:
        g = model.genes.get_by_id(gid)
        s = sweep.get(gid, {})
        pos = gid in positives
        cls = s.get("class_B_0.15", "")
        label = ("TP" if pos else "FP") if cls == "essential" else ("FN" if pos else "TN")
        if s.get("mtdna_out_of_scope") == "True" or s.get("untargeted_out_of_scope") == "True":
            label = "out of scope"
        ko = float(s["ko"]) if s.get("ko") not in (None, "") else None
        rows.append([rid, r.name, ";".join(sorted(r.compartments)), r.gene_reaction_rule, gid, pos,
                     audit.get(gid, {}).get("flag", "not in audit"),
                     ";".join(sorted(x.id for x in g.reactions)),
                     s.get("status", ""), "" if ko is None else round(ko / wt, 6), cls, label])
with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["reaction", "reaction_name", "compartment", "gene_rule", "gene", "on_consensus_list",
                "w29_tag_audit_flag", "gene_reactions_in_model", "ko_status", "ko_over_wt", "class_B_15pct",
                "confusion_label_15pct"])
    w.writerows(rows)
for row in rows:
    print(row)
