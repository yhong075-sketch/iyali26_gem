"""Configure input paths before importing the reference build chain."""

import argparse
import os
from pathlib import Path

from .model_layout import load_model_layout


def parse_args(argv=None):
    model = load_model_layout()
    parser = argparse.ArgumentParser(description="Rebuild iYali26 with the required tRNA-coupled biomass and metadata reaction-field selection (202 stoichiometries, 7 bounds, 9 GPRs)")
    parser.add_argument("--research-root", type=Path, help="Local reference tables and annotation cache")
    parser.add_argument("--starting-model", type=Path, default=model.start_model)
    parser.add_argument("--output-model", type=Path, default=model.canonical_model)
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
    parser.add_argument("--coq9-biomass-dilution", action="store_true",
                        help="Add the curated CoQ9 growth-dilution term (alpha is a user-selected, uncalibrated assumption); requires --energy-candidate E5")
    parser.add_argument("--coq9-respiratory-package", action="store_true",
                        help="Apply the curated 2026-10-08 CoQ9/respiratory-chain package; requires --energy-candidate E5 --coq9-curation qcycle --coq9-biomass-dilution")
    parser.add_argument("--r132-mito-adh", action="store_true",
                        help="Remove R132 (matrix alcohol dehydrogenase) per the curated 2026-10-13 decision; requires --coq9-respiratory-package")
    parser.add_argument("--restore-selection-cleared-gprs", action="store_true",
                        help="Re-apply the reviewed R612/R2274/R1302/R1303 GPRs and balanced R2274/R1303 equations that the metadata selection clears (curated 2026-10-19 decision); requires a complete selection")
    parser.add_argument("--restore-selection-overwritten-gprs", action="store_true",
                        help="Re-apply the reviewed R671 GPR without the mis-annotated GatB gene that the metadata selection re-adds (curated 2026-10-19 decision); requires a complete selection")
    parser.add_argument("--restore-selection-cpa-ura2", action="store_true",
                        help="Re-apply EG-GPR-001 on R159/R190 (URA2-only channelled R159, CPA2 and CPA1 on R190) after the metadata selection (curated 2026-10-19 decision); requires a complete selection")
    parser.add_argument("--remove-selection-tkl-fragment-gpr", action="store_true",
                        help="Remove the transketolase N-terminal fragment YALI1D02625g from R765/R766 after the metadata selection, keeping YALI1E07744g (curated 2026-10-19 decision); requires a complete selection")
    parser.add_argument("--remove-r153-peptidase-gpr", action="store_true",
                        help="Remove the A1 peptidase gene YALI1D17462g from R153 after the R153/R2176 merge, leaving R153 gene-less because no W29 argininosuccinate synthase gene is annotated (curated 2026-10-19 decision); the gene is kept unlinked")
    parser.add_argument("--r2041-forward", action="store_true",
                        help="Make R2041 (D-arabinokinase) forward-only after the metadata selection, keeping its balanced 2 H+ equation (curated 2026-10-19)")
    parser.add_argument("--r2190-ribulose-product", action="store_true",
                        help="Write D-ribulose instead of D-arabinose as R2190's product per the curated 2026-10-19 decision (R2190-RIBULOSE-20261019)")
    parser.add_argument("--r1964-gpr", action="store_true",
                        help="Give orphan R1964 (arabinitol dehydrogenase) the provisional GPR YALI1F03452g or YALI1E06616g and D-xylulose cross-references (curated 2026-10-19)")
    parser.add_argument("--gpi-chain-links", action="store_true",
                        help="Connect the blocked GPI-anchor chain R371-R376: gene-less step 1 (GlcNAc-PI synthesis), PE/DAG ER links, ER acetate export and an open GPI demand (curated 2026-10-21, GPI-CHAIN-LINKS-20261021)")
    parser.add_argument("--gpi-coproducts", action="store_true",
                        help="Add the missing co-products to the GPI chain: DAG on R373 and Dol-P on R375 (curated 2026-10-21, GPI-COPRODUCTS-20261021); requires --gpi-chain-links")
    parser.add_argument("--tpp-requirement", action="store_true",
                        help="Require TPP in biomass_C at 1.2e-6 mmol/gDW and close the gene-less TMP kinase R1342 (curated 2026-10-22, TPP-REQUIREMENT-20261022); requires --gpi-coproducts")
    parser.add_argument("--thiamine-phosphate-chemistry", action="store_true",
                        help="Give TMP, TPP, TTP and HET-P their KEGG fully protonated formulas and IDs, replacing the trichloromethoxypyridine mapping of TMP, and drop EC 2.7.4.9 from R1342 (curated 2026-10-22, THIAMINE-P-CHEMISTRY-20261022); requires --tpp-requirement")
    parser.add_argument("--thiamine-proton-lock", action="store_true",
                        help="Lock the pH 7.3 H+ of R752, R408 and R750 against the proton/water balance gate with the SBML-safe note lock_proton_water_stoichiometry: true (curated 2026-10-22, THIAMINE-PROTON-LOCK-20261022); requires --thiamine-phosphate-chemistry")
    parser.add_argument("--heme-iron-biomass", action="store_true",
                        help="Add heme a (1e-6) and iron(2+) (3.04e-5 mmol/gDW) to biomass_C at the Yeast8 consensus values (curated 2026-10-22, HEME-IRON-BIOMASS-20261022); requires --coq9-biomass-dilution")
    parser.add_argument("--r2193-iron-gpr", action="store_true",
                        help="Replace R2193's GPR (YALI1F07747g, a bHLH regulator, kept unlinked) with the two FET3/FTR1 pairs (YALI1D09277g and YALI1D09308g) or (YALI1D08684g and YALI1D08564g) (curated 2026-10-22, R2193-IRON-GPR-20261022)")
    parser.add_argument("--lipoylation-module", action="store_true",
                        help="Add the mitochondrial lipoylation route (LIP2-type, LIP5, LIP3-type), fill the ACP1 and lipoyl-H-protein formulas and add a provisional lipoyl growth-dilution requirement (curated 2026-10-20 decision)")
    parser.add_argument("--acyl-acp-chemistry", action="store_true",
                        help="Fill the acyl-ACP/ACP/fatty-acid formulas of the mtFAS and R2216-R2241 block, remove its surplus H+ and make its reversible cytosolic acyl-ACP hydrolases forward-only (curated 2026-10-20 decision); requires --lipoylation-module")
    parser.add_argument("--cofactor-biomass-terms", action="store_true",
                        help="Add THF (6.34e-5), FAD (1e-5), riboflavin (9.9e-4), NAD (2.65e-3), NADP (5.7e-4) and CoA (1.9e-4 mmol/gDW) to biomass_C at the yeast-GEM 9.1.1 r_4598 values (curated 2026-10-27, COFACTOR-BIOMASS-TERMS-20261027); requires --lipoylation-module")
    parser.add_argument("--r2135-acyl-coa-pool", action="store_true",
                        help="Replace acetyl-CoA by the ER-membrane pooled acyl-CoA m567 in R2135 and write R2145 as 1 PA + H2O -> 1 DAG + Pi (curated 2026-10-27, R2135-ACYL-COA-POOL-20261027); requires --cofactor-biomass-terms and runs after it")
    parser.add_argument("--r1726-cpt-open", action="store_true",
                        help="Open the inherited (0, 0) cholinephosphotransferase R1726 to [0, 1000] and keep R1725 closed (curated 2026-10-27, R1726-CPT-OPEN-20261027); requires --r2135-acyl-coa-pool and runs after it")
    parser.add_argument("--dolichol-route", action="store_true",
                        help="Join the split dolichol route: close the gene-less placeholder R2178, add the lumped gene-less dehydrodolichyl-PP -> dolichol step R_DOLSYN and drain dolichyl phosphate m503[C_er] from biomass_C at 1e-4 mmol/gDW (curated 2026-10-28, DOLICHOL-ROUTE-20261028); requires --r1726-cpt-open and runs after it")
    parser.add_argument("--um7-nuclear-diffusion", action="store_true",
                        help="Add the gene-less nucleus-cytosol diffusions of ADP-ribose (R_ADPRtn) and diphosphate (R_PPItn) in the Transport, nucleus group, joining the nuclear NAD block (gap-census UM7) to the cytosol (curated 2026-10-28, UM7-NUCLEAR-DIFFUSION-20261028); requires --dolichol-route and runs after it")
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
