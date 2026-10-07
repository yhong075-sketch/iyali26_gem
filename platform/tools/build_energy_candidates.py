#!/usr/bin/env python3
"""Build separate E0–E3 XML files from the pinned completed local reference."""

import argparse
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("IYALI26_RESEARCH_ROOT", str(ROOT / "artifacts/reference_pipeline_restore_20260909/research"))

from scripts.gem_annotate.energy_candidates import VARIANTS, build_candidate_file
from scripts.gem_annotate.execution import execution_limits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "model_metadata_trna_r1159_leak.xml")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--variants", nargs="+", choices=VARIANTS, default=["E0", "E1", "E2", "E3"])
    args = parser.parse_args()
    if not args.output_dir.resolve().is_relative_to(ROOT):
        raise ValueError("Candidate output must remain in this project workspace")
    with execution_limits(no_solve=True, allow_network=False):
        for variant in args.variants:
            record = build_candidate_file(args.source, args.output_dir / f"{variant}.xml", variant)
            print(variant, record["output_sha256"], "reload verified", flush=True)


if __name__ == "__main__":
    main()
