"""Small regression check for the dFBA new-false-negative definition."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

from cobra import Metabolite, Model, Reaction


REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from scripts.dfba_new_fn_essentiality import (  # noqa: E402
    compare_models,
    experimental_reference_coverage,
    load_dynamic_medium,
    load_experimental_reference,
    sha256,
    simulate_dfba,
)
from scripts.gem_annotate.validate_essential_genes import load_experimental  # noqa: E402


def write_toy_medium(path: Path, *, concentration: str = "10", pool_mode: str = "finite") -> None:
    initial_status = "nominal_formulation" if concentration else "unresolved"
    uptake_status = "inferred_upper_bound" if concentration else "permissive_upper_bound"
    uptake_basis = (
        "initial_concentration_mmol_l_times_10_divided_by_111_rounded"
        if concentration
        else "permissive_sensitivity_assumption_not_measured_kinetics"
    )
    path.write_text(
        "reaction_id,compound,initial_concentration_mmol_l,pool_mode,"
        "initial_concentration_status,max_uptake_mmol_gdw_h,"
        "uptake_evidence_status,uptake_basis,source_locator,source_accessed_on,notes\n"
        f"EX_glc,glucose,{concentration},{pool_mode},{initial_status},10,"
        f"{uptake_status},{uptake_basis},test_fixture,2026-08-20,toy medium\n",
        encoding="utf-8",
    )


def toy_model(with_bypass: bool) -> Model:
    model = Model("toy")
    external = Metabolite("glc_e", compartment="e")
    internal = Metabolite("glc_c", compartment="c")
    exchange = Reaction("EX_glc")
    exchange.add_metabolites({external: -1})
    exchange.bounds = (-10, 1000)

    transport = Reaction("T")
    transport.add_metabolites({external: -1, internal: 1})
    transport.gene_reaction_rule = "g_transport"
    biomass = Reaction("biomass_C")
    biomass.add_metabolites({internal: -1})
    biomass.bounds = (0, 1000)
    model.add_reactions([exchange, transport, biomass])
    if with_bypass:
        bypass = Reaction("BYPASS")
        bypass.add_metabolites({external: -1, internal: 1})
        model.add_reactions([bypass])
    model.objective = biomass
    return model


class DfbaNewFalseNegativeTests(unittest.TestCase):
    def test_gene_becomes_a_new_false_negative_when_candidate_adds_bypass(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            medium_path = root / "medium.csv"
            write_toy_medium(medium_path)
            baseline = toy_model(False)
            candidate = toy_model(True)
            baseline_bounds = {reaction.id: reaction.bounds for reaction in baseline.reactions}
            candidate_bounds = {reaction.id: reaction.bounds for reaction in candidate.reactions}
            rows, _, _ = compare_models(
                baseline,
                candidate,
                ["g_transport"],
                load_dynamic_medium(medium_path),
                {"hours": 1.0, "step_hours": 0.1, "initial_biomass_gdw_l": 0.01},
                0.01,
                root / "results.tsv",
            )

            self.assertEqual(len(rows), 1)
            self.assertTrue(rows[0]["baseline_predicted_essential"])
            self.assertFalse(rows[0]["candidate_predicted_essential"])
            self.assertTrue(rows[0]["new_false_negative"])
            self.assertEqual(
                {reaction.id: reaction.bounds for reaction in baseline.reactions}, baseline_bounds
            )
            self.assertEqual(
                {reaction.id: reaction.bounds for reaction in candidate.reactions}, candidate_bounds
            )
            with (root / "results.tsv").open(newline="") as stream:
                self.assertEqual(len(list(csv.DictReader(stream, delimiter="\t"))), 1)

            same_rows, _, _ = compare_models(
                toy_model(False),
                toy_model(False),
                ["g_transport"],
                load_dynamic_medium(medium_path),
                {"hours": 1.0, "step_hours": 0.1, "initial_biomass_gdw_l": 0.01},
                0.01,
                root / "same.tsv",
            )
            self.assertFalse(same_rows[0]["new_false_negative"])

    def test_conflicting_experimental_labels_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "essential.csv"
            path.write_text("gene_id,essential\ng1,1\ng1,0\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Conflicting duplicate"):
                load_experimental(path)
            path.write_text("gene_id,essential\n,1\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "blank gene_id"):
                load_experimental(path)

    def test_positive_only_reference_never_implies_negative_labels(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "positive_only.csv"
            path.write_text(
                "gene_id,source_gene_id,function,source,confidence\n"
                "g_transport,g_transport,transport protein,"
                "https://doi.org/10.1038/s42003-023-04996-8,"
                "consensus_essential_in_at_least_2_of_3_screens\n"
                "g_outside,g_outside,unmodeled protein,"
                "https://doi.org/10.1038/s42003-023-04996-8,"
                "consensus_essential_in_at_least_2_of_3_screens\n",
                encoding="utf-8",
            )
            positive_ids, contract = load_experimental_reference(path)
            tested_ids, coverage = experimental_reference_coverage(
                positive_ids, toy_model(False), toy_model(True)
            )

            self.assertEqual(positive_ids, ["g_outside", "g_transport"])
            self.assertEqual(tested_ids, ["g_transport"])
            self.assertTrue(contract["positive_only"])
            self.assertEqual(contract["negative_label_count"], 0)
            self.assertEqual(contract["unlisted_gene_semantics"], "unknown_not_nonessential")
            self.assertEqual(coverage["jointly_testable_positive_gene_count"], 1)
            self.assertEqual(coverage["untested_positive_gene_count"], 1)
            self.assertEqual(coverage["absence_semantics"], "unknown_not_nonessential")
            self.assertEqual(
                coverage["partitions"]["reference_gene_ids_in_neither_model"]["gene_ids"],
                ["g_outside"],
            )
            slurm = (REPOSITORY / "scripts/hpcc_dfba_new_fn.slurm").read_text(
                encoding="utf-8"
            )
            self.assertIn(
                "1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b",
                slurm,
            )
            self.assertIn("consensus_essential_genes.csv", slurm)

    def test_model_and_numeric_contracts_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            medium_path = root / "medium.csv"
            write_toy_medium(medium_path)
            medium = load_dynamic_medium(medium_path)
            with self.assertRaisesRegex(ValueError, "finite and positive"):
                simulate_dfba(
                    toy_model(False), medium, hours=float("inf"), step_hours=0.1,
                    initial_biomass_gdw_l=0.01,
                )

            candidate = toy_model(False)
            candidate.objective = candidate.reactions.T
            output = root / "invalid.tsv"
            with self.assertRaisesRegex(ValueError, "biomass_C"):
                compare_models(
                    toy_model(False), candidate, ["g_transport"], medium,
                    {"hours": 1.0, "step_hours": 0.1, "initial_biomass_gdw_l": 0.01},
                    0.01, output,
                )
            self.assertFalse(output.exists())

            duplicate_header = root / "duplicate_header.csv"
            write_toy_medium(duplicate_header)
            duplicate_header.write_text(
                duplicate_header.read_text().replace("reaction_id,compound", "reaction_id,reaction_id"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "exact columns"):
                load_dynamic_medium(duplicate_header)

            conflicting_status = root / "conflicting_status.csv"
            write_toy_medium(conflicting_status)
            conflicting_status.write_text(
                conflicting_status.read_text().replace(
                    ",finite,nominal_formulation,", ",finite,unresolved,"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "conflicts with pool_mode"):
                load_dynamic_medium(conflicting_status)

    def test_po1f_medium_contract_and_nondepleting_pool(self) -> None:
        medium = load_dynamic_medium(
            REPOSITORY / "data/media/po1f_csm_leu_dfba.csv"
        )
        by_id = {row["reaction_id"]: row for row in medium}
        self.assertEqual(len(by_id), 36)
        self.assertEqual(by_id["R1219"]["initial_concentration_mmol_l"], 0)
        self.assertEqual(by_id["R1219"]["max_uptake_mmol_gdw_h"], 0)
        self.assertEqual(by_id["R1219"]["pool_mode"], "closed")
        self.assertEqual(by_id["R1219"]["initial_concentration_status"], "formulation_absent")
        self.assertEqual(by_id["R1219"]["uptake_evidence_status"], "closed")
        inferred = {
            "R1189": (0.000740, 0.0000667),
            "R1003": (0.054300608, 0.00489),
            "R1202": (0.287026406, 0.02586),
            "R1204": (0.601051841, 0.05414),
            "R1215": (0.095406192, 0.00859),
            "R1217": (0.381184722, 0.03434),
            "R1220": (0.273747605, 0.02466),
            "R1222": (0.134039274, 0.01207),
            "R1223": (0.302681760, 0.02727),
            "R1231": (0.839489590, 0.07562),
            "R1232": (0.244822014, 0.02205),
            "R1233": (0.275953419, 0.02486),
            "R1234": (1.195049082, 0.10765),
            "R1354": (0.1784280489, 0.01607),
        }
        self.assertEqual(
            {
                reaction_id: (
                    row["initial_concentration_mmol_l"], row["max_uptake_mmol_gdw_h"]
                )
                for reaction_id, row in by_id.items()
                if row["uptake_evidence_status"] == "inferred_upper_bound"
            },
            inferred,
        )
        self.assertTrue(all(
            by_id[reaction_id]["uptake_basis"]
            == "initial_concentration_mmol_l_times_10_divided_by_111_rounded"
            for reaction_id in inferred
        ))
        self.assertIsNone(by_id["R1016"]["initial_concentration_mmol_l"])
        medium_path = REPOSITORY / "data/media/po1f_csm_leu_dfba.csv"
        self.assertIn(
            sha256(medium_path),
            (REPOSITORY / "scripts/hpcc_dfba_new_fn.slurm").read_text(encoding="utf-8"),
        )

        with tempfile.TemporaryDirectory() as directory:
            medium_path = Path(directory) / "nondepleting.csv"
            write_toy_medium(medium_path, concentration="", pool_mode="nondepleting")
            result = simulate_dfba(
                toy_model(False),
                load_dynamic_medium(medium_path),
                hours=0.1,
                step_hours=0.1,
                initial_biomass_gdw_l=0.01,
            )
            self.assertGreater(result["biomass_gain_gdw_l"], 0)
            self.assertEqual(result["final_concentrations_mmol_l"], {})


if __name__ == "__main__":
    unittest.main()
