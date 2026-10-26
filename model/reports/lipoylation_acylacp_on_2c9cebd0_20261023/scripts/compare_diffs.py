"""Compare two whole-model diffs (model_diff.py output) entity by entity (read-only).

Usage: python -I compare_diffs.py REFERENCE_DIFF.json NEW_DIFF.json OUT.json
REFERENCE = b3b9fe22 vs 164f235f (the two steps on their original base); NEW = candidate vs 2c9cebd0.
For every (section, item, field) it reports whether the entry is in both diffs, and whether the after
value and the before value agree. biomass_C stoichiometry is compared on the delta (after - before),
because the two bases differ in biomass_C by the TPP and heme/iron terms.
"""
import json
import sys

reference_path, new_path, out = sys.argv[1:4]


def key(c):
    return (c["section"], c["item"], c.get("field", c.get("change", "")))


def delta(c):
    before, after = c["before"] or {}, c["after"] or {}
    return {m: after.get(m, 0.0) - before.get(m, 0.0) for m in sorted(set(before) | set(after))
            if abs(after.get(m, 0.0) - before.get(m, 0.0)) > 1e-15}


reference = {key(c): c for c in json.load(open(reference_path))}
new = {key(c): c for c in json.load(open(new_path))}
rows = []
for k in sorted(set(reference) | set(new)):
    a, b = reference.get(k), new.get(k)
    row = {"section": k[0], "item": k[1], "field": k[2], "in_reference": a is not None, "in_new": b is not None}
    if a and b and "before" in a:
        if k == ("reactions", "biomass_C", "stoichiometry"):
            row["delta_equal"] = delta(a) == delta(b)
        else:
            row["after_equal"] = a["after"] == b["after"]
            row["before_equal"] = a["before"] == b["before"]
    rows.append(row)
summary = {
    "reference_entries": len(reference), "new_entries": len(new),
    "only_in_reference": [r for r in rows if not r["in_new"]],
    "only_in_new": [r for r in rows if not r["in_reference"]],
    "after_differs": [r for r in rows if r.get("after_equal") is False or r.get("delta_equal") is False],
    "before_differs": [r for r in rows if r.get("before_equal") is False],
    "rows": rows,
}
with open(out, "w") as handle:
    json.dump(summary, handle, indent=1)
    handle.write("\n")
for name in ("only_in_reference", "only_in_new", "after_differs", "before_differs"):
    print(name, len(summary[name]), [(r["section"], r["item"], r["field"]) for r in summary[name]])
print("entries", len(reference), len(new))
