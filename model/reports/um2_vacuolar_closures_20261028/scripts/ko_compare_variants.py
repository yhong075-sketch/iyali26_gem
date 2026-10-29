"""Checked single-gene knockouts of one model as built and under in-memory UM2 options, compared and scored
(read-only). Same logic as ../um7_nuclear_transports_20261028/scripts/ko_compare_two.py (checked_solve.py:
Gurobi 1e-9, a model reset before each knockout, the targeted check on the knocked reactions' flux; flagged
points get a fresh re-solve, at most FRESH_LIMIT), but the second "model" is the first with option bounds set
in memory (GROUPS of stage1_um2.py), so one file serves every option.

Usage (from platform/): PYTHONPATH=$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts python \
    ko_compare_variants.py MODEL.xml MEDIUM.csv PROFILE.json EXPERIMENTAL.csv OUT_PREFIX VARIANT[,VARIANT...]
Class: essential when KO/WT is strictly below the cutoff. Negatives are the user-defined class (model genes not
on the consensus list); mtDNA genes are out of scope; the list was used in development (not independent).
"""
import csv, json, math, sys, time, warnings
from pathlib import Path
from checked_solve import SOLVER_TOLERANCE, configure, matrix, solve
from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
from scripts.gem_annotate.model_layout import MODEL
from tools.dfba_new_fn_essentiality import load_experimental_reference
sys.path.insert(0, str(Path(__file__).parent))
GROUPS = {
    "pump": {"R795": (0, 1000)}, "water": {"R1363": (0, 1000)},
    "hdiff": {"R1163": (0, 1000)}, "hleak": {"R1163": (-1000, 0)},
    "aaimp": {r: (0, 1000) for r in ("R869", "R878", "R887", "R897")},
    "aaexp": {r: (0, 1000) for r in ("R870", "R879", "R883", "R888", "R898")},
    "aspglu": {r: (0, 1000) for r in ("R871", "R876")},
    "pa": {r: (0, 1000) for r in ("R918", "R925", "R927")},
}
CUTOFFS = (0.01, 0.05, 0.10, 0.15)
CHANGE_TOLERANCE = 1e-6
FRESH_LIMIT = 20
model_path, medium, profile, experimental, prefix, variants = sys.argv[1:7]
variants = ["as_built"] + variants.split(",")
positives = set(load_experimental_reference(Path(experimental))[0])
mtdna = set(json.loads((MODEL.root / "expected/mtdna_genes_out_of_scope.json").read_text())["genes"])
fresh_used = []


def context():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ctx = load_effective_simulation_context(model_path=Path(model_path), media_path=Path(medium),
                                                strain_profile_path=Path(profile))
    configure(ctx.model)
    return ctx


def apply(model, label):
    for part in label.split("+"):
        if part != "as_built":
            for rid, b in GROUPS[part].items():
                model.reactions.get_by_id(rid).bounds = b


def closed(model, gene):
    return {r.id for r in model.genes.get_by_id(gene).reactions if r.bounds == (0, 0)}


def fresh(label, gene):
    ctx = context()
    apply(ctx.model, label)
    with ctx.model:
        ctx.model.genes.get_by_id(gene).knock_out()
        return solve(ctx.model, matrix(ctx.model), closed(ctx.model, gene))


def sweep(label):
    ctx = context()
    model = ctx.model
    apply(model, label)
    stoich = matrix(model)
    wt = solve(model, stoich, set())
    rows = {}
    for gene in sorted(g.id for g in model.genes if g.id not in ctx.excluded_runtime_genes):
        with model:
            model.genes.get_by_id(gene).knock_out()
            row = solve(model, stoich, closed(model, gene))
        row["fresh_resolve"] = None
        if not row["passed"] and row["status"] != "infeasible":
            if len(fresh_used) >= FRESH_LIMIT:
                raise RuntimeError(f"More than {FRESH_LIMIT} knockouts need a fresh re-solve; stopping as declared")
            fresh_used.append([label, gene])
            row["fresh_resolve"] = fresh(label, gene)
        rows[gene] = row
    print(f"{label}: WT {wt['status']} {wt['growth']!r} passed={wt['passed']} genes={len(rows)}", flush=True)
    return ctx, wt, rows


def final(row):
    for point in (row, row.get("fresh_resolve") or {}):
        if point.get("passed"):
            return point["growth"]
    return None


def call(wt, g, cutoff):
    if g is None or not wt:
        return "unresolved"
    return "essential" if g / wt < cutoff else "non-essential"


def score(wt, rows):
    res = {}
    for cutoff in CUTOFFS:
        c = {"TP": 0, "FN": 0, "FP": 0, "TN": 0, "unresolved": 0}
        for gene, row in rows.items():
            if gene in mtdna:
                continue
            k = call(wt, final(row), cutoff)
            c["unresolved" if k == "unresolved" else (("TP" if gene in positives else "FP") if k == "essential"
                                                      else ("FN" if gene in positives else "TN"))] += 1
        tp, fn, fp, tn = c["TP"], c["FN"], c["FP"], c["TN"]
        d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        res[str(cutoff)] = {**c, "recall": tp / (tp + fn), "specificity": tn / (tn + fp),
                            "precision": tp / (tp + fp), "mcc": (tp * tn - fp * fn) / d if d else None}
    return res


started = time.monotonic()
runs = {v: sweep(v) for v in variants}
assert all(wt["passed"] for _, wt, _ in runs.values()), "A WT point failed the check"
base_ctx, base_wt, base_rows = runs["as_built"]
summary = {"condition": f"Task A medium + PO1f, Gurobi {SOLVER_TOLERANCE}, reset per knockout, targeted check",
           "model_sha256": base_ctx.canonical_model_sha256, "medium_sha256": base_ctx.medium_sha256,
           "strain_profile_sha256": base_ctx.strain_profile_sha256,
           "negative_class": "user-defined: model genes not on the consensus list; mtDNA genes out of scope",
           "validation_status": "development reference, not independent validation",
           "groups": {k: {r: list(b) for r, b in v.items()} for k, v in GROUPS.items()},
           "fresh_resolves": fresh_used, "variants": {}}
for v, (ctx, wt, rows) in runs.items():
    vs = {"wt": wt, "scores": score(wt["growth"], rows),
          "non_optimal": sorted(g for g, r in rows.items() if r["status"] != "optimal"),
          "failing_check": sorted(g for g, r in rows.items() if r["status"] == "optimal" and not r["passed"]),
          "unresolved": sorted(g for g, r in rows.items() if final(r) is None)}
    if v != "as_built":
        grow, cls, gain = [], {str(c): [] for c in CUTOFFS}, []
        for gene in sorted(rows):
            gb, gv = final(base_rows[gene]), final(rows[gene])
            if (gb is None) != (gv is None) or (gb is not None and abs(gb - gv) > CHANGE_TOLERANCE):
                grow.append({"gene": gene, "as_built": gb, v: gv})
            if gb is not None and gv is not None and gv - gb > CHANGE_TOLERANCE:
                gain.append({"gene": gene, "as_built": gb, v: gv, "gain": gv - gb})
            for c in CUTOFFS:
                p, q = call(base_wt["growth"], gb, c), call(wt["growth"], gv, c)
                if p != q:
                    cls[str(c)].append({"gene": gene, "experimental_positive": gene in positives, "as_built": p, v: q})
        vs.update({"wt_delta": wt["growth"] - base_wt["growth"], "growth_changed_over_1e-6": grow,
                   "ko_gain_over_1e-6": gain, "class_changed": cls})
    summary["variants"][v] = vs
    print(v, "WT", wt["growth"], "15%:", vs["scores"]["0.15"], "| changed", len(vs.get("growth_changed_over_1e-6", [])),
          "| class changes", {c: len(x) for c, x in vs.get("class_changed", {}).items()}, flush=True)
with open(prefix + ".tsv", "w", newline="") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["gene", "experimental_positive", "mtdna_out_of_scope"] + [f"{v}_{k}" for v in variants for k in ("status", "passed", "ko", "max_knocked_flux")])
    for gene in sorted(base_rows):
        line = [gene, gene in positives, gene in mtdna]
        for v in variants:
            r = runs[v][2][gene]
            line += [r["status"], r["passed"], final(r), r.get("max_knocked_flux")]
        w.writerow(line)
summary["lp_count"] = sum(len(r) + 1 for _, _, r in runs.values()) + len(fresh_used)
summary["seconds"] = time.monotonic() - started
json.dump(summary, open(prefix + "_summary.json", "w"), indent=1)
print("LPs", summary["lp_count"], "seconds", round(summary["seconds"]))
