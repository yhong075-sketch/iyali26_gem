"""Check that a combined candidate's diff from the baseline is exactly the union of its parts' diffs.

Usage: python union_diff.py BASELINE.xml COMBINED.xml OUT.json PART.xml [PART.xml ...]
Each PART is a candidate built from the same baseline with one change set (stacked parts are fine).
For every (section, item, field) a part changes, the expected after-value is that part's value; when
several parts change the same key to different values the key is listed as a conflict. Gene
reaction-membership lists are compared as the union of the reactions each part adds or removes.
Prints and writes: keys only in the combined diff, keys missing from it, and value mismatches.
"""
import json
import os
import sys

os.environ.setdefault("HOME", os.path.abspath("cache"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cobra.io import read_sbml_model
from model_diff import table


def flat(t):
    out = {}
    for section, items in t.items():
        if section == "model":
            for k, v in items.items():
                out[(section, k, "")] = v
            continue
        for item, fields in items.items():
            for field, value in fields.items():
                out[(section, item, field)] = value
    return out


def changes(base, other):
    keys = set(base) | set(other)
    return {k: other.get(k, "<absent>") for k in keys if base.get(k, "<absent>") != other.get(k, "<absent>")}


def key(k):
    return "|".join(k)


if __name__ == "__main__":
    base_path, combined_path, out_path, *parts = sys.argv[1:]
    base = flat(table(read_sbml_model(base_path)))
    combined = changes(base, flat(table(read_sbml_model(combined_path))))
    expected, owners, conflicts = {}, {}, []
    membership = {}  # gene -> (added reactions, removed reactions)
    for part in parts:
        diff = changes(base, flat(table(read_sbml_model(part))))
        for k, v in diff.items():
            if k[0] == "genes" and k[2] == "reactions":
                before = set(base.get(k, []) if base.get(k) != "<absent>" else [])
                after = set(v if v != "<absent>" else [])
                add, rem = membership.setdefault(k, (set(), set()))
                add |= after - before
                rem |= before - after
                owners.setdefault(key(k), []).append(os.path.basename(part))
                continue
            if k in expected and expected[k] != v:
                conflicts.append({"key": key(k), "parts": owners[key(k)] + [os.path.basename(part)]})
            expected.setdefault(k, v)
            owners.setdefault(key(k), []).append(os.path.basename(part))
    for k, (add, rem) in membership.items():
        before = set(base.get(k, []) if base.get(k) != "<absent>" else [])
        expected[k] = sorted((before | add) - rem)
    combined_norm = {k: (sorted(v) if k[0] == "genes" and k[2] == "reactions" and v != "<absent>" else v)
                     for k, v in combined.items()}
    only_combined = sorted(key(k) for k in combined_norm if k not in expected)
    missing = sorted(key(k) for k in expected if k not in combined_norm)
    mismatch = sorted(key(k) for k in expected if k in combined_norm and combined_norm[k] != expected[k])
    report = {"combined_changes": len(combined_norm), "expected_changes": len(expected),
              "only_in_combined": only_combined, "missing_from_combined": missing,
              "value_mismatch": mismatch, "part_conflicts": conflicts,
              "owners": {key(k): owners.get(key(k)) for k in expected}}
    json.dump(report, open(out_path, "w"), indent=1, default=str)
    for name in ("combined_changes", "expected_changes"):
        print(name, report[name])
    for name in ("only_in_combined", "missing_from_combined", "value_mismatch", "part_conflicts"):
        print(name, len(report[name]), report[name][:20])
