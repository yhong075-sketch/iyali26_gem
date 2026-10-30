"""Every claim level's ECO term exists in the pinned release, is not obsolete and sits in its branch."""

import csv

from scripts.gem_annotate.model_layout import MODEL

from evidence_helpers import W

TABLE = MODEL.curation_file("eco_terms.tsv")


def load_terms():
    with TABLE.open() as handle:
        header = handle.readline()
        terms = {row["id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    return header, terms


def ancestors(term_id, terms):
    seen, stack = set(), [term_id]
    while stack:
        row = terms.get(stack.pop())
        for parent in filter(None, (row["is_a"] if row else "").split(";")):
            if parent not in seen:
                seen.add(parent)
                stack.append(parent)
    return seen


def test_table_matches_pinned_release():
    header, terms = load_terms()
    assert header.startswith(f"# ECO {W.eco['release']};")
    assert len(terms) > 1000


def test_level_eco_terms_exist_are_current_and_in_their_branch():
    _, terms = load_terms()
    for level, info in W.levels.items():
        assert info.eco in terms, f"{level.value}: {info.eco} not in ECO"
        assert terms[info.eco]["obsolete"] == "false", f"{level.value}: {info.eco} is obsolete"
        root = W.eco["branches"][info.branch]
        assert info.eco == root or root in ancestors(info.eco, terms), \
            f"{level.value}: {info.eco} is not under {info.branch} ({root})"


def test_branch_roots_are_distinct_so_a_mislabel_fails():
    _, terms = load_terms()
    experimental = W.eco["branches"]["experimental"]
    similarity = W.eco["branches"]["similarity"]
    assert experimental not in ancestors("ECO:0000250", terms)
    assert similarity not in ancestors("ECO:0000314", terms)
