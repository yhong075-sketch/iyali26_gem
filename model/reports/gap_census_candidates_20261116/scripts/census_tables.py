"""Cause table for the inactive reactions of candidate 7ccbb568, and the changes since the 2026-10-26 census.

Usage:
    python census_tables.py ART_DIR REPORT_DIR OCTOBER_ART_DIR OCTOBER_REPORT_DIR

Reads   ../input_hashes.csv  SHA-256 of the model, medium, profile and census script copies (checked first;
                          needs IYALI26_RESEARCH_ROOT)
        ART_DIR/out/      outputs of the unchanged census scripts (gapfind.py, analyze_roots.py, fixtests.py,
                          fva_gurobi.py) for 7ccbb568 (open, taskA_po1f) and a1237821 (open, regression)
        OCTOBER_*         the 2026-10-26 census outputs (regression reference and October causes)
Writes  REPORT_DIR/*.csv  small tables
        ART_DIR/tables/   full tables
Exits non-zero when a gate fails; nothing is written in that case.
"""
import collections
import hashlib
import json
import os
import sys

import pandas as pd

ART, REP, OCT_ART, OCT_REP = sys.argv[1:5]
sys.path.insert(0, os.path.join(ART, "scripts"))
from kernel import assign_primary_cause, load_gap_tables, medium_dependent_blocked, module_table  # noqa: E402

OUT, TAG, REF = os.path.join(ART, "out"), "7ccbb568", "a1237821"

# Structural root count made at plan time on 7ccbb568 (TASK.md); the census root table must reproduce it.
EXPECTED_ROOTS = {"RNP": 98, "RNC": 101, "RNP+RNC (single reversible reaction)": 137,
                  "isolated (RNP+RNC: in no open reaction)": 28}

# The user's two questions, plus the two groups they do not cover. A "coupling" reaction has no dead-end
# metabolite in its module; whether a step is missing there is not decided by the census.
GROUP = {
    "closed in model (0,0)": "closed by bounds",
    "blocked by a closed (0,0) reaction": "closed by bounds",
    "cofactor pseudo-gap (cofactor not drained by biomass)": "steady state: cofactor never drained",
    "dead-end root metabolite": "dead end: something is missing",
    "missing transport": "dead end: something is missing",
    "missing exchange": "dead end: something is missing",
    "duplicate metabolite id": "dead end: something is missing",
    "chain disconnected at several roots (needs >1 fix)": "dead end: something is missing",
    "stoichiometric coupling / pseudo-gap (no root)": "no dead end found (coupling)",
}


def fail(message):
    print("GATE FAILED:", message)
    sys.exit(1)


def blocked_ids(table):
    return set(table[table.blocked].id)


def check_input_hashes(who, artifacts=""):
    """Stop unless every input listed for this script in input_hashes.csv has the recorded SHA-256."""
    here = os.path.dirname(os.path.abspath(__file__))
    bases = {"repo": os.path.abspath(os.path.join(here, "..", "..", "..", "..")),
             "research": os.environ.get("IYALI26_RESEARCH_ROOT", ""), "artifacts": artifacts}
    table = pd.read_csv(os.path.join(here, "..", "input_hashes.csv"))
    paths = {}
    for row in table[table.checked_by.isin([who, "both"])].itertuples():
        if not bases[row.base]:
            fail(f"no {row.base} folder to locate {row.key}")
        path = os.path.join(bases[row.base], row.path)
        with open(path, "rb") as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
        if digest != row.sha256:
            fail(f"{row.key} has SHA-256 {digest}, expected {row.sha256}")
        paths[row.key] = path
    print(f"gate passed: {len(paths)} input hashes equal input_hashes.csv")
    return paths


# ---------------------------------------------------------------- gates
check_input_hashes("census_tables", ART)
summary_today = json.load(open(f"{OUT}/{REF}_open_summary.json"))
summary_october = json.load(open(f"{OCT_ART}/raw_outputs/{REF}_open_summary.json"))
if summary_today != summary_october:
    fail(f"regression summary differs from October: {summary_today} vs {summary_october}")
reg_today = pd.read_csv(f"{OUT}/{REF}_open_reactions.csv")
reg_october = pd.read_csv(f"{OCT_ART}/raw_outputs/{REF}_open_reactions.csv")
if blocked_ids(reg_today) != blocked_ids(reg_october):
    fail("regression blocked set differs from October")

R, Mt, S, summ = load_gap_tables(OUT, TAG, "open")
Rm, Mm, Sm, summ_m = load_gap_tables(OUT, TAG, "taskA_po1f")
roots = pd.read_csv(f"{OUT}/{TAG}_open_roots.csv").fillna("")
pooled = pd.read_csv(f"{OUT}/{TAG}_open_pooled_roots.csv")
root_counts = roots.root_class.value_counts().to_dict()
if root_counts != EXPECTED_ROOTS:
    fail(f"root table {root_counts} differs from the plan-time count {EXPECTED_ROOTS}")

fix_out = json.load(open(f"{OUT}/{TAG}_open_fixtests_out.json"))
drain = json.load(open(f"{OUT}/{TAG}_open_draintests_out.json"))
unresolved = [t for t in fix_out + drain if "solver failed" in t.get("note", "")]
if unresolved:
    fail(f"{len(unresolved)} fix or drain tests ended with a solver failure: "
         + ", ".join(f"{t['kind']} {t['root']}" for t in unresolved[:10]))
print("gates passed: regression equals October; root table", root_counts, "; no failed fix or drain test")

# ---------------------------------------------------------------- one cause per blocked reaction
drain_out = {t["root"]: {"unblocked": t["unblocked"]} for t in drain}
CA = assign_primary_cause(R, Mt, fix_out, drain_out)
CA["cause_group"] = CA.cause.map(GROUP)
if CA.cause_group.isna().any():
    fail(f"cause without a group: {sorted(CA[CA.cause_group.isna()].cause.unique())}")

met_um = Mt.set_index("id").um.to_dict()
pooled_new = pooled[~pooled.already_root] if len(pooled) else pooled
pooled_modules = {met_um.get(m, "") for copies in pooled_new.copies for m in str(copies).split(";")} - {""}

G = pd.read_csv(f"{OUT}/{TAG}_open_fva_gurobi.csv", index_col=0)
Gm = pd.read_csv(f"{OUT}/{TAG}_taskA_po1f_fva_gurobi.csv", index_col=0)

oct_blocked = pd.read_csv(f"{OCT_REP}/blocked_reactions.csv")
oct_cause = oct_blocked.set_index("id").cause.to_dict()
oct_model = set(reg_today.id)
oct_blocked_ids = blocked_ids(reg_today)

Rb = R[R.blocked].merge(CA, on="id", how="left")
Rb["closed_0_0"] = (Rb.lb == 0) & (Rb.ub == 0)
Rb["has_gpr"] = Rb.gpr != ""
Rb["pooled_root_in_module"] = Rb.um.isin(pooled_modules)
Rb["blocked_on_taskA"] = Rb.id.map(Rm.set_index("id").blocked)
Rb["gurobi_absmax_open"] = Rb.id.map(G.absmax)
Rb["in_october_model"] = Rb.id.isin(oct_model)
Rb["blocked_in_october"] = Rb.id.map(lambda x: (x in oct_blocked_ids) if x in oct_model else "")
Rb["cause_october"] = Rb.id.map(oct_cause).fillna("")
full_cols = ["id", "name", "cause_group", "cause", "detail", "um", "closed_0_0", "has_gpr", "gpr", "subsystem",
             "compartments", "boundary", "lb", "ub", "struct_blocked", "pooled_root_in_module", "blocked_on_taskA",
             "gurobi_absmax_open", "in_october_model", "blocked_in_october", "cause_october", "equation"]
Rb = Rb[full_cols]

# ---------------------------------------------------------------- changes since October
now_model, now_blocked = set(R.id), blocked_ids(R)
now_cause = CA.set_index("id").cause.to_dict()
info = R.set_index("id")
info_oct = reg_today.set_index("id")
changes = []
for rid in sorted(oct_model | now_model):
    if rid not in now_model:
        change = "removed from the model (was blocked)" if rid in oct_blocked_ids else "removed from the model (was active)"
    elif rid not in oct_model:
        change = "added since October, blocked" if rid in now_blocked else "added since October, active"
    elif rid in oct_blocked_ids and rid not in now_blocked:
        change = "became active"
    elif rid not in oct_blocked_ids and rid in now_blocked:
        change = "became blocked"
    elif rid in now_blocked and oct_cause.get(rid, "") != now_cause.get(rid, ""):
        change = "still blocked, cause changed"
    else:
        continue
    src = info if rid in now_model else info_oct
    changes.append(dict(id=rid, change=change, name=src.at[rid, "name"], cause_october=oct_cause.get(rid, ""),
                        cause_now=now_cause.get(rid, ""), lb_now=info.at[rid, "lb"] if rid in now_model else "",
                        ub_now=info.at[rid, "ub"] if rid in now_model else "",
                        lb_october=info_oct.at[rid, "lb"] if rid in oct_model else "",
                        ub_october=info_oct.at[rid, "ub"] if rid in oct_model else ""))
changes = pd.DataFrame(changes)

# ---------------------------------------------------------------- solver check, medium, modules, tests
solver_rows = []
for cond, glpk_table, gurobi_table in (("open", R, G), ("taskA_po1f", Rm, Gm)):
    glpk_set = blocked_ids(glpk_table)
    gurobi_set = set(gurobi_table[gurobi_table["blocked_1e-9"]].index)
    for rid in sorted(glpk_set ^ gurobi_set):
        solver_rows.append(dict(condition=cond, id=rid, blocked_glpk=rid in glpk_set, blocked_gurobi_1e9=rid in gurobi_set,
                                gurobi_absmax=gurobi_table.absmax.get(rid)))
    print(f"solver check {cond}: GLPK {len(glpk_set)} blocked, Gurobi (1e-9) {len(gurobi_set)}, "
          f"differ on {len(glpk_set ^ gurobi_set)}")
solver = pd.DataFrame(solver_rows, columns=["condition", "id", "blocked_glpk", "blocked_gurobi_1e9", "gurobi_absmax"])

medium_only, new_rnp = medium_dependent_blocked(R, Mt, Rm, Mm)
modules = module_table(R, Mt, roots, CA, fix_out)

fix_table = pd.DataFrame([dict(kind=t["kind"], root=t["root"], um=t["um"], op=json.dumps(t["op"]), n_um_rxns=t["n_um_rxns"],
                               n_unblocked=t["n_unblocked"], unblocked=";".join(t["unblocked"]), note=t["note"])
                          for t in fix_out])
fix_summary = fix_table.groupby("kind").agg(tests=("root", "size"), tests_unblocking=("n_unblocked", lambda s: int((s > 0).sum())),
                                            largest=("n_unblocked", "max")).reset_index()
drain_table = pd.DataFrame([dict(cofactor=t["root"], drained=";".join(t["op"]["drain"]), n_blocked_tested=t["n_um_rxns"],
                                 n_unblocked=t["n_unblocked"], unblocked=";".join(t["unblocked"]), note=t["note"])
                            for t in drain])

# ---------------------------------------------------------------- counts
def count_row(name, s, table, mets):
    c = s["classes"]
    return dict(run=name, reactions=s["n_rxn"], blocked=s["blocked"], blocked_pct=round(100 * s["blocked"] / s["n_rxn"], 1),
                metabolites=s["n_met"], gap_metabolites=s["gap"], unconnected_modules=s["UM"],
                largest_modules=" ".join(map(str, s["UM_sizes"][:8])), RNP=c.get("RNP", 0), RNC=c.get("RNC", 0),
                single_reversible=c.get("RNP+RNC (single reversible reaction)", 0), isolated=c.get("isolated", 0),
                blocked_explained_structurally=s["blocked_explained_structurally"],
                closed_0_0=int(((table.cond_lb == 0) & (table.cond_ub == 0)).sum()),
                blocked_with_gpr=int((table.blocked & (table.gpr != "")).sum()),
                blocked_boundary=int((table.blocked & (table.boundary != "")).sum()),
                gap_without_class=int((mets.gap & (mets.cls == "")).sum()))


Mr = pd.read_csv(f"{OUT}/{REF}_open_metabolites.csv")
Mr["cls"] = Mr["cls"].fillna("")
reg_today["gpr"] = reg_today["gpr"].fillna("")
reg_today["boundary"] = reg_today["boundary"].fillna("")
counts = pd.DataFrame([count_row("7ccbb568 open (all exchanges open)", summ, R, Mt),
                       count_row("7ccbb568 Task A medium + PO1f", summ_m, Rm, Mm),
                       count_row("a1237821 open (regression, equals October)", summary_today, reg_today, Mr)])

cause_counts = (Rb.groupby(["cause_group", "cause"]).size().rename("blocked_reactions").reset_index()
                .merge(oct_blocked.groupby("cause").size().rename("october_a1237821").reset_index(), on="cause", how="outer")
                .fillna({"blocked_reactions": 0, "october_a1237821": 0}))
cause_counts["cause_group"] = cause_counts.cause.map(GROUP)
cause_counts = cause_counts.astype({"blocked_reactions": int, "october_a1237821": int}).sort_values(["cause_group", "cause"])

lp = collections.OrderedDict()
lp["regression scan a1237821 open (GLPK)"] = 2 * len(reg_today)
lp["scans 7ccbb568 open + Task A (GLPK)"] = 2 * 2 * len(R)
lp["scans 7ccbb568 open + Task A (Gurobi)"] = 2 * 2 * len(R)
lp["single-fix tests"] = 2 * int(fix_table.n_um_rxns.sum())
lp["cofactor drain tests"] = 2 * int(drain_table.n_blocked_tested.sum())
lp_table = pd.DataFrame([dict(part=k, lps=v) for k, v in lp.items()] + [dict(part="check 1 total", lps=sum(lp.values()))])

# ---------------------------------------------------------------- write
os.makedirs(REP, exist_ok=True)
os.makedirs(f"{ART}/tables", exist_ok=True)
counts.to_csv(f"{REP}/summary_counts.csv", index=False)
cause_counts.to_csv(f"{REP}/cause_counts.csv", index=False)
Rb.drop(columns=["equation", "detail"]).to_csv(f"{REP}/blocked_reactions.csv", index=False)
roots.to_csv(f"{REP}/root_metabolites.csv", index=False)
pooled.to_csv(f"{REP}/pooled_roots.csv", index=False)
changes.to_csv(f"{REP}/cause_changes_vs_a1237821.csv", index=False)
solver.to_csv(f"{REP}/solver_check.csv", index=False)
fix_summary.to_csv(f"{REP}/fix_test_summary.csv", index=False)
drain_table.drop(columns=["unblocked"]).to_csv(f"{REP}/cofactor_drains.csv", index=False)
lp_table.to_csv(f"{REP}/lp_count_check1.csv", index=False)
modules.drop(columns=["reactions"]).to_csv(f"{REP}/unconnected_modules.csv", index=False)
Rb.to_csv(f"{ART}/tables/blocked_reactions_full.csv", index=False)
fix_table.to_csv(f"{ART}/tables/fix_tests.csv", index=False)
drain_table.to_csv(f"{ART}/tables/cofactor_drains_full.csv", index=False)
modules.to_csv(f"{ART}/tables/unconnected_modules_full.csv", index=False)
medium_only.to_csv(f"{ART}/tables/medium_only_blocked.csv", index=False)
Mt[Mt.gap].to_csv(f"{ART}/tables/gap_metabolites.csv", index=False)

print(json.dumps(dict(
    blocked_open=int(R.blocked.sum()), blocked_taskA=int(Rm.blocked.sum()), reactions=len(R),
    by_group=Rb.cause_group.value_counts().to_dict(), by_cause=Rb.cause.value_counts().to_dict(),
    coupling_with_pooled_root_in_module=int(((Rb.cause_group == "no dead end found (coupling)") & Rb.pooled_root_in_module).sum()),
    changes=changes.change.value_counts().to_dict() if len(changes) else {},
    medium_only_blocked=len(medium_only), roots=len(roots), pooled_new=len(pooled_new),
    fix_tests=len(fix_table), drain_tests=len(drain_table), lps_check1=int(sum(lp.values()))), indent=1))
