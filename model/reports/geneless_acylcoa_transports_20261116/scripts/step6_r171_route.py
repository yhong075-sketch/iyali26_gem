"""Step 6, mechanism: why the ATP-synthase knockouts grow 0.7 % more under the option (4 LPs).

All 12 genes whose KO/WT moves in the sweep knock out only R171 (ATP synthase). Here R171 is knocked out through one
of them (YALI1F04433g) on 7ccbb568 as is (B0/K) and with the option (B0/ALL), Task A + PO1f, Gurobi 1e-9. Each gets
one pFBA (growth, then minimum total flux at that growth). The fluxes that differ by more than 1e-6 are written out.
One minimum-flux optimum is one route, not proof that the route is unique.
Usage: python step6_r171_route.py  ->  out/step6_r171_route.json, tables/step6_r171_route.csv
"""
import csv
from pathlib import Path

from common import LP_COUNT, dump, load, pfba
from step4_options import apply_base, apply_variant

GENE = "YALI1F04433g"
model = load("taskA", "gurobi")
res = {}
for key in ["B0/K", "B0/ALL"]:
    with model:
        base, name = key.split("/")
        apply_base(model, base)
        apply_variant(model, name)
        model.genes.get_by_id(GENE).knock_out()
        growth, status, flux = pfba(model)
        res[key] = {"growth": growth, "status": status, "flux": flux}
a, b = res["B0/K"]["flux"], res["B0/ALL"]["flux"]
rows = []
for rid in sorted(set(a) | set(b)):
    x, y = a.get(rid, 0.0), b.get(rid, 0.0)
    if abs(x - y) > 1e-6:
        r = model.reactions.get_by_id(rid) if rid in model.reactions else None
        rows.append({"reaction": rid, "name": r.name if r else "(option only)", "as_is": x, "option": y,
                     "change": y - x, "rule": r.gene_reaction_rule if r else ""})
rows.sort(key=lambda r: -abs(r["change"]))
dest = Path(__file__).resolve().parents[1] / "tables/step6_r171_route.csv"
with open(dest, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["reaction", "name", "as_is", "option", "change", "rule"])
    w.writeheader()
    w.writerows(rows)
out = {k: {"growth": v["growth"], "status": v["status"]} for k, v in res.items()}
out.update({"gene": GENE, "differing": len(rows), "lp_count": LP_COUNT["n"],
            "flux": {k: {r: f for r, f in v["flux"].items() if abs(f) > 1e-9} for k, v in res.items()}})
print({k: out[k] for k in ("B0/K", "B0/ALL", "differing", "lp_count")}, "->", dump(out, "step6_r171_route.json"))
