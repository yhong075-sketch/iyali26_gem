"""Compare two pytest junit XML files test by test (outcome per classname::name)."""
import sys
import xml.etree.ElementTree as ET


def outcomes(path):
    out = {}
    for case in ET.parse(path).getroot().iter("testcase"):
        key = f"{case.get('classname')}::{case.get('name')}"
        tags = {child.tag for child in case}
        out[key] = ("failure" if "failure" in tags else "error" if "error" in tags
                    else "skipped" if "skipped" in tags else "passed")
    return out


a, b = outcomes(sys.argv[1]), outcomes(sys.argv[2])
print("before", len(a), "after", len(b))
for key in sorted(set(a) | set(b)):
    if a.get(key) != b.get(key):
        print(f"{key}: {a.get(key)} -> {b.get(key)}")
