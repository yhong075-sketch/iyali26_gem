"""Static model evidence for gene-less metabolic reactions (read-only).

Usage: python orphan_evidence.py --model M.xml --flux flux_evidence.csv --ids ids.csv
       --selection metadata_reaction_selection.json --out evidence.csv
Per reaction: annotations, mass/charge balance, exact and compartment-free duplicates
(with their GPRs), GPR cleared by the metadata selection, and for each metabolite whether
another reaction can produce / consume it (dead-end check).
"""
import argparse, json, os
os.environ.setdefault("HOME", os.path.abspath("cache"))
import pandas as pd
import cobra

ap = argparse.ArgumentParser()
for a in ("--model", "--flux", "--ids", "--selection", "--out"):
    ap.add_argument(a, required=True)
a = ap.parse_args()

M = cobra.io.read_sbml_model(a.model)
flux = pd.read_csv(a.flux).set_index("id")
ids = pd.read_csv(a.ids)
sel = json.load(open(a.selection))["reactions"]


def ann(r, key):
    v = r.annotation.get(key, "")
    return ";".join(v) if isinstance(v, list) else str(v)


def sig(r, by):
    return frozenset((by(m), round(c, 6)) for m, c in r.metabolites.items())


def neg(r, by):
    return frozenset((by(m), -round(c, 6)) for m, c in r.metabolites.items())


exact = lambda m: m.id
nocomp = lambda m: (m.name, m.formula)  # compartment-free identity


def can_make(m, skip):
    """Some other reaction can produce m (respecting bounds)."""
    for r in m.reactions:
        if r.id == skip:
            continue
        c = r.metabolites[m]
        if (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0):
            return True
    return False


def can_use(m, skip):
    for r in m.reactions:
        if r.id == skip:
            continue
        c = r.metabolites[m]
        if (c < 0 and r.upper_bound > 0) or (c > 0 and r.lower_bound < 0):
            return True
    return False


index_exact, index_nocomp = {}, {}
for r in M.reactions:
    index_exact.setdefault(sig(r, exact), []).append(r.id)
    index_nocomp.setdefault(sig(r, nocomp), []).append(r.id)

rows = []
for rid in ids.id:
    r = M.reactions.get_by_id(rid)
    bal = r.check_mass_balance()
    missing = [m.id for m in r.metabolites if not m.formula]
    dup_exact = [x for x in index_exact.get(sig(r, exact), []) + index_exact.get(neg(r, exact), []) if x != rid]
    dup_nocomp = [x for x in index_nocomp.get(sig(r, nocomp), []) + index_nocomp.get(neg(r, nocomp), [])
                  if x != rid and x not in dup_exact]
    fmt = lambda xs: "; ".join(f"{x} [{M.reactions.get_by_id(x).gene_reaction_rule or 'no GPR'}]" for x in xs)
    # dead-end check, in the direction(s) the reaction's bounds allow
    no_source, no_sink = [], []
    for m, c in r.metabolites.items():
        consumed = (c < 0 and r.upper_bound > 0) or (c > 0 and r.lower_bound < 0)
        produced = (c > 0 and r.upper_bound > 0) or (c < 0 and r.lower_bound < 0)
        if consumed and not produced and not can_make(m, rid):
            no_source.append(f"{m.id} {m.name}")
        if produced and not consumed and not can_use(m, rid):
            no_sink.append(f"{m.id} {m.name}")
    s = sel.get(rid, {})
    cleared = bool("gpr" in s.get("fields", []) and s["before"].get("gpr") and not s["after"].get("gpr"))
    f = flux.loc[rid]
    rows.append(dict(
        id=rid, name=r.name, spontaneous_name_flag=bool(ids.set_index("id").loc[rid, "spontaneous_note"]),
        compartments=";".join(sorted(r.compartments)), lb=r.lower_bound, ub=r.upper_bound,
        ec=ann(r, "ec-code"), kegg_reaction=ann(r, "kegg.reaction"), rhea=ann(r, "rhea"),
        metanetx=ann(r, "metanetx.reaction"),
        equation=r.build_reaction_string(use_metabolite_names=True),
        mass_balanced=not {k: v for k, v in bal.items() if k != "charge" and abs(v) > 1e-6} and not missing,
        charge_balanced=abs(bal.get("charge", 0)) < 1e-6,
        imbalance=json.dumps({k: round(v, 4) for k, v in bal.items() if abs(v) > 1e-6}),
        missing_formula=";".join(missing),
        exact_duplicates=fmt(dup_exact), compartment_free_duplicates=fmt(dup_nocomp),
        selection_cleared_gpr=s["before"]["gpr"] if cleared else "",
        selection_changed_fields=";".join(s.get("fields", [])),
        no_other_source=" | ".join(no_source), no_other_sink=" | ".join(no_sink),
        flux_status=f.flux_status, pfba_flux=f.pfba_flux,
        fva_min_0pct_opt=f.fva_min_0pct_opt, fva_max_0pct_opt=f.fva_max_0pct_opt,
        **{"fva_min_99.9pct_opt": f["fva_min_99.9pct_opt"], "fva_max_99.9pct_opt": f["fva_max_99.9pct_opt"]},
        ko_growth=f.ko_growth, ko_wt_ratio=f.ko_wt_ratio,
    ))
out = pd.DataFrame(rows)
out.to_csv(a.out, index=False)
print("rows", len(out))
print("mass-imbalanced or missing formula:", out[~out.mass_balanced].id.tolist())
print("charge-imbalanced:", out[~out.charge_balanced].id.tolist())
print("exact duplicates:", out[out.exact_duplicates != ""][["id", "exact_duplicates"]].values.tolist())
print("compartment-free duplicates:", out[out.compartment_free_duplicates != ""][["id", "compartment_free_duplicates"]].values.tolist())
print("selection cleared GPR:", out[out.selection_cleared_gpr != ""][["id", "selection_cleared_gpr"]].values.tolist())
print("blocked with no other source/sink:", (out.flux_status.str.startswith("blocked") & ((out.no_other_source != "") | (out.no_other_sink != ""))).sum())
