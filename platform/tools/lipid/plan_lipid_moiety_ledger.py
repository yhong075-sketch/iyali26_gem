#!/usr/bin/env python3
"""CLI for the read-only lipid moiety-ledger manifest compiler."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

try:
    from .lipid_moiety_ledger import (
        LedgerError,
        atomic_write_manifest,
        compile_dry_run,
    )
except ImportError:
    from lipid_moiety_ledger import LedgerError, atomic_write_manifest, compile_dry_run


ROOT = Path(__file__).resolve().parents[1]


def _json_output_path(value: str) -> Path:
    path = Path(value)
    suffix = path.suffix.lower()
    if suffix in {".xml", ".sbml"}:
        raise argparse.ArgumentTypeError(
            "SBML/XML output is forbidden; this compiler emits JSON manifests only"
        )
    if suffix != ".json":
        raise argparse.ArgumentTypeError("--output must have a .json suffix")
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compile a read-only in-memory TAG acyl-moiety ledger manifest."
    )
    parser.add_argument(
        "model_xml", nargs="?", default=str(ROOT / "model.xml"), help="source SBML"
    )
    parser.add_argument(
        "curation_csv",
        nargs="?",
        default=str(ROOT / "data" / "lipid_combo_curation.csv"),
        help="lipid_combo_curation.csv input",
    )
    parser.add_argument(
        "--spec",
        default=str(ROOT / "data" / "lipid_moiety_ledger_spec.json"),
        help="read-only planning contract included in input SHA-256",
    )
    parser.add_argument(
        "--output",
        type=_json_output_path,
        help="optional JSON manifest path; otherwise write canonical JSON to stdout",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = compile_dry_run(args.model_xml, args.curation_csv, spec_path=args.spec)
    except LedgerError as error:
        parser.error(str(error))
    if args.output is None:
        sys.stdout.buffer.write(result.canonical_json)
    else:
        try:
            atomic_write_manifest(
                args.output,
                result.canonical_json,
                {
                    "model.xml": Path(args.model_xml),
                    "data/lipid_combo_curation.csv": Path(args.curation_csv),
                    "data/lipid_moiety_ledger_spec.json": Path(args.spec),
                },
                result.manifest["input_sha256"],
            )
        except LedgerError as error:
            parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
