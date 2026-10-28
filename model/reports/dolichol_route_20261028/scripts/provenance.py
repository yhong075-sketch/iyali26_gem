"""Compare the dolichol-route reactions and species across the source models (read-only)."""
import sys, warnings
import cobra
warnings.simplefilter("ignore")
IDS = [f"R{i}" for i in range(217, 236)] + ["R2177","R2178","R2179","R2180","R2181","R2182","R287","R1101","R288","R2102","R2106","R2105","R2138","R289","R2101","R173","R141"]
METS = ["m10","m32","m35","m82","m109","m203","m406","m500","m1249","m413","m343","m268","m1961","m1962","m1965","m1966","m1967","m389","m102","m952","m2018","m2019","m2020","m2021","m499","m501","m503","m502","m411","m412","m432","m1968","m1969","m1989","m344","m345"]
paths = {"ce9e1112": sys.argv[1], "iyali26_src": "../model/source/iyali26.xml", "iyli21": "../model/source/iyli21.xml", "legacy": "../model/candidates/legacy/model.xml"}
models = {}
for k, p in paths.items():
    try:
        models[k] = cobra.io.read_sbml_model(p)
    except Exception as e:
        print("LOAD FAIL", k, e)
base = models["ce9e1112"]
print("=== species identities (ce9e1112) ===")
for mid in METS:
    for m in base.metabolites.query(lambda x: x.id.split("[")[0] == mid):
        print(f"{m.id:14s} {m.name[:70]:70s} {m.formula} q={m.charge}")
print("\n=== reaction comparison (equation|bounds|gpr) ===")
for rid in IDS:
    rows = {}
    for k, M in models.items():
        if rid in M.reactions:
            r = M.reactions.get_by_id(rid)
            rows[k] = (r.reaction, r.bounds, r.gene_reaction_rule)
        else:
            rows[k] = None
    same = len({str(v) for v in rows.values()}) == 1
    print(f"{rid}: {'IDENTICAL in all 4' if same else 'DIFFERS'}")
    if not same:
        for k, v in rows.items():
            print(f"    {k:12s} {v}")
print("\n=== species formula comparison ===")
for mid in METS:
    rows = {}
    for k, M in models.items():
        hits = M.metabolites.query(lambda x: x.id.split("[")[0] == mid)
        rows[k] = sorted((m.id, m.formula, m.charge) for m in hits)
    same = len({str(v) for v in rows.values()}) == 1
    if not same:
        print(mid, "DIFFERS")
        for k, v in rows.items():
            print("   ", k, v)
print("done")
