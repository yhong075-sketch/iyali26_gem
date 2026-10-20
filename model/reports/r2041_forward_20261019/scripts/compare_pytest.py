"""Compare two pytest JUnit XML files test by test.

Usage: python compare_pytest.py BEFORE.xml AFTER.xml

Copied from model/reports/selection_gpr_fix_20261019/scripts/ (worktree claude/clever-kare-7d2e1a); only this line was added.
"""
import sys
import xml.etree.ElementTree as ET
from collections import Counter


def outcomes(path):
    result = {}
    for case in ET.parse(path).iter("testcase"):
        key = f"{case.get('classname')}::{case.get('name')}"
        tags = {child.tag for child in case}
        result[key] = ("error" if "error" in tags else "failed" if "failure" in tags
                       else "skipped" if "skipped" in tags else "passed")
    return result


before, after = (outcomes(p) for p in sys.argv[1:3])
print("before:", dict(Counter(before.values())))
print("after: ", dict(Counter(after.values())))
for key in sorted(set(before) | set(after)):
    if before.get(key) != after.get(key):
        print(f"{before.get(key, 'absent'):>8} -> {after.get(key, 'absent'):<8} {key}")
