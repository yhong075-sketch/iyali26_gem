"""Generate the committed ECO term table from a pinned ECO OBO release.

    python -m evidence.eco_table <eco.obo> <out.tsv>

The table keeps id, name, is_a parents and the obsolete flag for every term, with the OBO's
data-version and SHA-256 in a header line, so tests can check ECO IDs without the OBO file.
"""

import argparse
import csv
import hashlib
from pathlib import Path

FIELDS = ["id", "name", "is_a", "obsolete"]


def parse_obo(text: str) -> tuple[str, list[dict]]:
    """(data-version, terms) from OBO text; only [Term] stanzas, is_a without the comment."""
    version = ""
    terms: list[dict] = []
    current: dict | None = None
    for line in text.splitlines():
        if line.startswith("data-version: ") and not terms and current is None:
            version = line.removeprefix("data-version: ").strip()
        if line == "[Term]":
            current = {"id": "", "name": "", "is_a": [], "obsolete": False}
            terms.append(current)
            continue
        if line.startswith("["):
            current = None
            continue
        if current is None:
            continue
        if line.startswith("id: "):
            current["id"] = line[4:].strip()
        elif line.startswith("name: "):
            current["name"] = line[6:].strip()
        elif line.startswith("is_a: "):
            current["is_a"].append(line[6:].split(" ")[0])
        elif line.strip() == "is_obsolete: true":
            current["obsolete"] = True
    return version, sorted(terms, key=lambda t: t["id"])


def write_table(obo_path: Path, out_path: Path) -> Path:
    raw = obo_path.read_bytes()
    version, terms = parse_obo(raw.decode("utf-8"))
    if not version or not terms:
        raise ValueError(f"{obo_path} has no data-version or no terms")
    with out_path.open("w", newline="") as handle:
        handle.write(f"# ECO {version}; source {obo_path.name} "
                     f"sha256 {hashlib.sha256(raw).hexdigest()}; {len(terms)} terms\n")
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for term in terms:
            writer.writerow({"id": term["id"], "name": term["name"], "is_a": ";".join(sorted(term["is_a"])),
                             "obsolete": "true" if term["obsolete"] else "false"})
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("obo", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    print(write_table(args.obo, args.out))


if __name__ == "__main__":
    main()
