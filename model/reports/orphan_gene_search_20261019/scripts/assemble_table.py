"""Assemble gene_search.csv: one row per reaction.

Sequence columns (best W29 hit, reciprocal result, targeting) are computed from the BLAST and
predictor tables. Judgment columns (proposal, tier, (a)/(b)/(c), gene naming, paralogue note,
sources) come from curation_calls.json, written by hand for this task. Counts are printed from the
assembled rows.

Usage: python assemble_table.py REPORT_DIR
"""
import collections
import csv
import json
import os
import sys

import pandas as pd


def best_hit(h, rid):
    g = h[h.reaction == rid].sort_values("bitscore", ascending=False)
    if not len(g):
        return "none at E<=1e-5"
    x = g.iloc[0]
    return (f"{x.model_gene_id} ({x.protein_id}) vs {x.query_acc}: {x.pident:.1f}% id, qcov {x.qcov:.0f}%, "
            f"scov {x.scov:.0f}%, E={x.evalue:.0e}; in model: {bool(x.in_model)}")


def reciprocal(rc, rid, genes):
    out = []
    for gene in genes:
        g = rc[(rc.reaction == rid) & (rc.model_gene_id == gene)]
        if not len(g):
            continue
        tops = "; ".join(f"{t.split(' OS=')[0]} [{t.split('OS=')[1].split(' OX=')[0] if 'OS=' in t else ''}] {p:.0f}%"
                         for t, p in zip(g.sp_title.head(3), g.pident.head(3)))
        out.append(f"{gene}: {tops}")
    return " | ".join(out) if out else "not run (no W29 candidate)"


def main(rdir):
    h = pd.read_csv(os.path.join(rdir, "blast_w29_hits.csv"))
    rc = pd.read_csv(os.path.join(rdir, "reciprocal_hits.csv"))
    ctx = {r["id"]: r for r in csv.DictReader(open(os.path.join(rdir, "reaction_context.csv")))}
    calls = json.load(open(os.path.join(rdir, "curation_calls.json")))
    q = pd.read_csv(os.path.join(rdir, "queries.csv")).fillna("")
    rows = []
    for rid, c in calls.items():
        x = ctx[rid]
        qs = q[(q.reaction == rid) & (q.accession != "")]
        rows.append(dict(
            id=rid, group=x["group"], name=x["name"], compartment=x["compartments"], equation=x["equation_ids"],
            flux_status_sdleu=x["flux_status"],
            proposal=c["proposal"], proposed_gpr=c.get("proposed_gpr", ""), evidence_tier=c["tier"],
            a_enzyme_exists=c["a"], b_compartment_direction=c["b"], c_substrate_source=c["c"],
            query_proteins="; ".join(f"{a} ({b.split(';')[-1].strip()})" for a, b in zip(qs.accession, qs.basis)) or "none: " + q[q.reaction == rid].basis.iloc[0],
            best_w29_hit=best_hit(h, rid) + (("; " + c["extra_hit"]) if c.get("extra_hit") else ""),
            reciprocal_result=reciprocal(rc, rid, c.get("reciprocal_genes", [])) + (("; " + c["extra_reciprocal"]) if c.get("extra_reciprocal") else ""),
            targeting=c.get("targeting", ""), paralog_note=c.get("paralog_note", ""), genes_named=c.get("genes_named", ""),
            sources_opened=c.get("sources", ""), flags_outside_scope=c.get("flags", "")))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(rdir, "gene_search.csv"), index=False)
    print("rows", len(df), "groups", dict(collections.Counter(df.group)))
    print("proposal", dict(collections.Counter(df.proposal)))
    print("tier", dict(collections.Counter(df.evidence_tier)))
    print(df.groupby(["proposal", "evidence_tier"]).id.apply(lambda s: ", ".join(s)).to_string())


if __name__ == "__main__":
    main(sys.argv[1])
