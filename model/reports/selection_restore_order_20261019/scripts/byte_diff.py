"""Locate byte differences between two SBML files at the reaction-block level (read-only).

Usage: python byte_diff.py A.xml B.xml
Splits each file into <reaction> blocks and the text outside them, then reports: whether the
outside text is identical, whether the reaction blocks come in the same order, and for each block
that differs, whether its lines are identical once sorted (that is, a pure reordering).
"""
import re
import sys

BLOCK = re.compile(r'<reaction [^>]*id="([^"]+)".*?</reaction>', re.S)


def split(text):
    blocks = {m.group(1): m.group(0) for m in BLOCK.finditer(text)}
    return BLOCK.sub("<reaction/>", text), list(blocks), blocks


a_text, b_text = (open(p, encoding="utf-8").read() for p in sys.argv[1:3])
(a_rest, a_order, a_blocks), (b_rest, b_order, b_blocks) = split(a_text), split(b_text)
print("identical bytes:", a_text == b_text)
print("text outside reaction blocks identical:", a_rest == b_rest)
print("reaction blocks:", len(a_order), len(b_order), "same order:", a_order == b_order)
for rid in a_order:
    if a_blocks[rid] != b_blocks.get(rid):
        a_lines, b_lines = a_blocks[rid].splitlines(), b_blocks.get(rid, "").splitlines()
        print(f"differs: {rid}; identical when sorted: {sorted(a_lines) == sorted(b_lines)}")
        for side, lines in (("A", a_lines), ("B", b_lines)):
            refs = [re.search(r'species="([^"]+)"', l).group(1) for l in lines if "speciesReference" in l]
            print(f"  {side} speciesReference order: {', '.join(refs)}")
