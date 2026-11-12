"""Test-by-test comparison of two pytest junit XML files (before = archived HEAD, after = HEAD + this change)."""
import sys, json
import xml.etree.ElementTree as ET

def outcomes(path):
    out = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        key = f"{case.get('classname')}::{case.get('name')}"
        kind = "passed"
        for child in case:
            if child.tag in ("failure", "error", "skipped"):
                kind = {"failure": "failed", "error": "error", "skipped": "skipped"}[child.tag]
                break
        out[key] = kind
    return out

before, after = outcomes(sys.argv[1]), outcomes(sys.argv[2])
count = lambda d: {k: sum(1 for v in d.values() if v == k) for k in ("passed", "failed", "error", "skipped")}
changed = {k: (before.get(k), after.get(k)) for k in sorted(set(before) | set(after)) if before.get(k) != after.get(k)}
res = {"before": count(before), "after": count(after), "changed": changed}
json.dump(res, open(sys.argv[3], "w"), indent=1)
print(json.dumps(res, indent=1))
