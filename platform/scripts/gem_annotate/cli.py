"""Configure input paths before importing the reference build chain."""

import argparse
import os
from pathlib import Path


def parse_args(argv=None):
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description="Rebuild iYali26 with the required tRNA-coupled biomass and metadata reaction-field selection (202 stoichiometries, 7 bounds, 9 GPRs)")
    parser.add_argument("--research-root", type=Path, help="Local reference tables and annotation cache")
    parser.add_argument("--starting-model", type=Path, default=root / "data" / "iyali26.xml")
    parser.add_argument("--output-model", type=Path, default=root / "model.xml")
    parser.add_argument("--mnx-dir", type=Path, help="Override the local MetaNetX table directory")
    parser.add_argument("--cache-dir", type=Path, help="Override the annotation cache directory")
    parser.add_argument("--offline", action="store_true", help="Use cached annotations without network requests")
    parser.add_argument("--no-solve", action="store_true", help="Skip all optimization diagnostics; require a new output path")
    parser.add_argument("--coq9-curation", choices=("off", "metadata", "qcycle"), default="metadata",
                        help="metadata (default): annotation/Boolean corrections; qcycle: explicit R305 proton candidate; off: bypass these corrections")
    parser.add_argument("--coq9-functional-gpr", action="store_true",
                        help="Opt-in R695 COQ7 AND COQ9 substrate-access hypothesis; requires separate offline/no-solve metadata build")
    parser.add_argument("--coq-c5-gpr", action="store_true",
                        help="Opt-in mitochondrial NADH-coupled R39 with COQ6 AND YAH1-like AND ARH1-like GPR; offline/no-solve only")
    parser.add_argument("--coq-literature-revision", action="store_true",
                        help="Opt-in R19/R18 quinol revision, including the C5 and COQ9 candidates; partial pathway with unresolved redox gap; offline/no-solve only")
    parser.add_argument("--canonical-copy", action="store_true", help="Build the complete reference chain before final reaction selection (also the ordinary-build default); require a new output path")
    parser.add_argument("--provisional-capacity-profile", type=Path)
    parser.add_argument("--trna-biomass-mode", choices=("split",))
    parser.add_argument("--r608-curation", type=Path, help="Existing explicit R608 candidate; disabled by default")
    parser.add_argument("--vatpase-gpr-hypothesis", action="store_true",
                        help="Test three common AND dependencies in R794/R795; unvalidated hypothesis, requires offline/no-solve and a new output path")
    parser.add_argument("--energy-candidate", choices=("E0", "E1", "E2", "E3", "E4", "E5"), default="E0",
                        help="Optional energy repair: E1 directions, E2 chemistry, E3 both, E4 NDP1/R72, E5 CAT2p/OAADCm; offline/no-solve and new output required")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    for argument, variable in (("research_root", "IYALI26_RESEARCH_ROOT"),
                               ("mnx_dir", "IYALI26_MNX_DIR"),
                               ("cache_dir", "IYALI26_CACHE_DIR")):
        value = getattr(args, argument)
        if value is not None:
            os.environ[variable] = str(value.resolve())
    from .main import build_model
    return build_model(args)
