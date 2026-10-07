#!/usr/bin/env python3
"""Build a separate E5 vacuole candidate only when explicitly enabled."""

import argparse
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.gem_annotate.execution import execution_limits
from scripts.gem_annotate.vacuole_candidates import build_candidate_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", "--baseline-model", type=Path, required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--output", "--candidate-model", type=Path, required=True)
    parser.add_argument("--enabled", action="store_true")
    args = parser.parse_args()
    with execution_limits(no_solve=True, allow_network=False) as attempts:
        if not args.source.resolve().is_relative_to(ROOT):
            raise ValueError("Source must remain in the project workspace")
        if hashlib.sha256(args.source.read_bytes()).hexdigest() != args.source_sha256:
            raise ValueError("Explicit source SHA256 differs")
        result = build_candidate_file(args.source, args.output, args.enabled)
        print(result["candidate"], result["output_sha256"], "reload verified", attempts)


if __name__ == "__main__":
    main()
