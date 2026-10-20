"""Read-only model context for the gene search.

Writes model_genes.csv (gene, rxns) for build_curation_dbs, and reaction_context.csv with, for each
target reaction: equation, annotations, other reactions that share an EC number (with GPRs), and the
static producers of each substrate (criterion (c), no flux run).

Usage: python model_context.py MODEL.xml TRIAGE.csv OUTDIR
"""
import csv
import os
import sys

os.environ.setdefault("HOME", os.path.expanduser("~"))
import cobra  # noqa: E402

GROUPS = ("10 ", "7 ", "6 ")


def as_list(v):
    if v is None:
        return []
    return list(v) if isinstance(v, (list, tuple, set)) else [v]


def producers(model, met, exclude):
    """Reactions other than `exclude` that can make `met` under their current bounds."""
    out = []
    for r in met.reactions:
        if r.id == exclude:
            continue
        c = r.metabolites[met]
        if (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0):
            out.append(r.id + ("[" + r.gene_reaction_rule + "]" if r.gene_reaction_rule else "[noGPR]"))
    return out


def main(model_path, triage_path, outdir):
    model = cobra.io.read_sbml_model(model_path)
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "model_genes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene", "rxns"])
        for g in model.genes:
            w.writerow([g.id, ";".join(sorted(r.id for r in g.reactions))])

    by_ec = {}
    for r in model.reactions:
        for ec in as_list(r.annotation.get("ec-code")):
            by_ec.setdefault(ec, []).append(r)

    rows = list(csv.DictReader(open(triage_path)))
    targets = [r for r in rows if r["recommended_group"].startswith(GROUPS)]
    out = []
    for t in targets:
        r = model.reactions.get_by_id(t["id"])
        ecs = as_list(r.annotation.get("ec-code"))
        same_ec = sorted({x.id + "(" + ",".join(sorted(x.compartments)) + ")[" + (x.gene_reaction_rule or "noGPR") + "]"
                          for ec in ecs for x in by_ec.get(ec, []) if x.id != r.id})
        # Substrates in the direction(s) the bounds allow.
        subs = [m for m, c in r.metabolites.items() if (c < 0 and r.upper_bound > 0) or (c > 0 and r.lower_bound < 0)]
        src = []
        for m in subs:
            p = producers(model, m, r.id)
            src.append(m.id + " " + m.name + ": " + (str(len(p)) + " producers: " + "; ".join(p[:6]) + (" ..." if len(p) > 6 else "") if p else "NO OTHER PRODUCER"))
        out.append(dict(
            id=r.id, group=t["recommended_group"].split()[0], name=r.name,
            equation_ids=r.reaction, bounds=f"[{r.lower_bound}, {r.upper_bound}]",
            compartments=",".join(sorted(r.compartments)),
            ec=";".join(ecs), kegg=";".join(as_list(r.annotation.get("kegg.reaction"))),
            rhea=";".join(as_list(r.annotation.get("rhea"))),
            metanetx=";".join(as_list(r.annotation.get("metanetx.reaction"))),
            flux_status=t["flux_status"], flags=t["flags"],
            same_ec_reactions=" | ".join(same_ec), substrate_sources=" || ".join(src)))
    with open(os.path.join(outdir, "reaction_context.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(len(model.genes), "genes;", len(out), "target reactions")


if __name__ == "__main__":
    main(*sys.argv[1:4])
