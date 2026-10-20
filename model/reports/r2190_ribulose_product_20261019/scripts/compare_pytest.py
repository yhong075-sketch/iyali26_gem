"""Compare two pytest JUnit XML files test by test (outcome per test id)."""
import json
import sys
import xml.etree.ElementTree as ET
from collections import Counter


def outcomes(path):
    result = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        key = f"{case.get('classname')}::{case.get('name')}"
        tags = {child.tag for child in case}
        result[key] = ("error" if "error" in tags else "failed" if "failure" in tags
                       else "skipped" if "skipped" in tags else "passed")
    return result


before, after = outcomes(sys.argv[1]), outcomes(sys.argv[2])
report = {
    "counts_before": dict(Counter(before.values())), "counts_after": dict(Counter(after.values())),
    "only_before": sorted(set(before) - set(after)),
    "only_after": {k: after[k] for k in sorted(set(after) - set(before))},
    "changed": {k: [before[k], after[k]] for k in sorted(set(before) & set(after)) if before[k] != after[k]},
    "non_passing_before": {k: v for k, v in sorted(before.items()) if v in ("failed", "error")},
}
json.dump(report, open(sys.argv[3], "w"), indent=1)
print(json.dumps({k: report[k] for k in ("counts_before", "counts_after", "only_before", "only_after", "changed")}, indent=1))
