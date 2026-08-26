"""Small regression check for the dFBA new-false-negative definition."""

from __future__ import annotations

import argparse
import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from cobra import Metabolite, Model, Reaction


REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

import scripts.dfba_new_fn_essentiality as dfba  # noqa: E402
from scripts.dfba_new_fn_essentiality import (  # noqa: E402
    RESULT_COLUMNS,
    _id_digest,
    _shared_contract,
    aggregate_shards,
    compare_models,
    experimental_reference_coverage,
    load_dynamic_medium,
    load_experimental_reference,
    sha256,
    shard_gene_ids,
    simulate_dfba,
    write_json,
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


def toy_model(
    with_bypass: bool, *, requires_flux: bool = False, include_closed_leucine: bool = False
) -> Model:
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
    if include_closed_leucine:
        leucine = Metabolite("leu_e", compartment="e")
        closed_exchange = Reaction("R1219")
        closed_exchange.add_metabolites({leucine: -1})
        closed_exchange.bounds = (-1000, 1000)
        model.add_reactions([closed_exchange])
    if requires_flux:
        maintenance = Reaction("MAINTENANCE")
        maintenance.add_metabolites({internal: -1})
        maintenance.bounds = (1, 1000)
        model.add_reactions([maintenance])
    if with_bypass:
        bypass = Reaction("BYPASS")
        bypass.add_metabolites({external: -1, internal: 1})
        model.add_reactions([bypass])
    model.objective = biomass
    return model


def write_toy_rescue_medium(path: Path) -> None:
    write_toy_medium(path, concentration="0")
    with path.open("a", encoding="utf-8") as stream:
        stream.write(
            "R1219,L-leucine,0,closed,formulation_absent,0,closed,"
            "zero_by_formulation,test_fixture,2026-08-20,closed leucine\n"
        )


def array_context(gene_ids: list[str]) -> dict:
    return {
        "git_commit": "test-commit",
        "inputs": {"test": {"path": "/test", "sha256": "0" * 64}},
        "medium": [],
        "positive_gene_ids": gene_ids,
        "gene_ids": gene_ids,
        "baseline": None,
        "candidate": None,
        "settings": {
            "hours": 24.0,
            "step_hours": 0.1,
            "initial_biomass_gdw_l": 0.05,
            "growth_cutoff": 0.01,
            "solver": "glpk",
        },
        "model_contracts": {"baseline": {"fingerprint": "baseline"}, "candidate": {"fingerprint": "candidate"}},
        "dynamic_medium_contract": {"row_count": 1},
        "experimental_reference_contract": {"mode": "positive_only_consensus"},
        "experimental_reference_coverage": {"jointly_testable_positive_gene_count": len(gene_ids)},
        "software": {"python": "test", "cobra": "test", "solver_interface": "test"},
    }


def result_row(gene_id: str, *, new_false_negative: bool = False) -> dict:
    baseline_essential = new_false_negative
    candidate_essential = not new_false_negative
    return {
        "gene_id": gene_id,
        "baseline_ko_biomass_gain_gdw_l": "0",
        "candidate_ko_biomass_gain_gdw_l": "0",
        "baseline_ko_to_wt_gain_ratio": "0",
        "candidate_ko_to_wt_gain_ratio": "0",
        "baseline_predicted_essential": str(baseline_essential),
        "candidate_predicted_essential": str(candidate_essential),
        "new_false_negative": str(new_false_negative),
        "baseline_reaction_ids": "R1",
        "candidate_reaction_ids": "R1",
        "gpr_evidence_status": "model/GPR assignment only",
    }


def write_shard_fixture(
    root: Path,
    context: dict,
    *,
    shard_index: int,
    shard_count: int,
    new_false_negative_ids: set[str] | None = None,
) -> None:
    new_false_negative_ids = new_false_negative_ids or set()
    gene_ids = shard_gene_ids(
        context["gene_ids"], shard_index=shard_index, shard_count=shard_count
    )
    shard = root / "shards" / f"{shard_index:03d}"
    shard.mkdir()
    results = shard / "gene_results.tsv"
    rows = [
        result_row(gene_id, new_false_negative=gene_id in new_false_negative_ids)
        for gene_id in gene_ids
    ]
    with results.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=RESULT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    new_ids = [row["gene_id"] for row in rows if row["new_false_negative"] == "True"]
    execution = {
        "mode": "shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "evaluated_gene_ids": gene_ids,
        "evaluated_gene_id_sha256": _id_digest(gene_ids),
        "universe_gene_count": len(context["gene_ids"]),
        "universe_gene_id_sha256": _id_digest(context["gene_ids"]),
        "final_scientific_gate_deferred_to_aggregate": True,
    }
    summary = {
        "schema_version": dfba.SCHEMA_VERSION,
        "status": "complete",
        **_shared_contract(context),
        "execution": execution,
        "wild_type": {
            "baseline": {"status": "optimal", "steps": 240, "biomass_gain_gdw_l": 1.0},
            "candidate": {"status": "optimal", "steps": 240, "biomass_gain_gdw_l": 1.0},
        },
        "evaluated_gene_count": len(rows),
        "new_false_negative_count": len(new_ids),
        "new_false_negative_gene_ids": new_ids,
        "no_new_false_negative_gate_passed": not new_ids,
        "final_scientific_gate_passed": None,
        "gene_results": {"path": str(results.resolve()), "sha256": sha256(results)},
    }
    summary_path = shard / "summary.json"
    write_json(summary_path, summary)
    write_json(shard / "wrapper_status.json", {
        "mode": "shard",
        "array_task_id": str(shard_index),
        "shard_count": str(shard_count),
        "preflight_exit_code": 0,
        "runner_exit_code": 0,
        "postflight_exit_code": 0,
        "technical_execution_completed": True,
        "final_scientific_gate_deferred_to_aggregate": True,
        "wrapper_gate_passed": True,
        "summary_sha256": sha256(summary_path),
        "gene_results_sha256": sha256(results),
    })


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
            self.assertIn(sha256(REPOSITORY / "scripts/dfba_new_fn_essentiality.py"), slurm)
            candidate = REPOSITORY / "model_lipid_unlump_strict_sn_candidate_bc2aac8f_r989_r1521_r39_coq_provisional.xml"
            self.assertIn(candidate.name, slurm)
            self.assertIn(sha256(candidate), slurm)
            self.assertIn(sha256(REPOSITORY / "scripts/lp_sn12_candidate.py"), slurm)
            self.assertIn("DFBA_MODE", slurm)
            self.assertIn("wt_diagnostic", slurm)
            self.assertIn("--wt-diagnostic", slurm)

    def test_wt_diagnostic_records_pfba_infeasibility_without_screening_genes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            medium_path = root / "medium.csv"
            write_toy_medium(medium_path, concentration="0")
            medium = load_dynamic_medium(medium_path)
            with self.assertRaises(dfba.DfbaInfeasibleError) as caught:
                simulate_dfba(
                    toy_model(False, requires_flux=True),
                    medium,
                    hours=0.1,
                    step_hours=0.1,
                    initial_biomass_gdw_l=0.01,
                    model_identity={
                        "role": "baseline",
                        "input_path": "/baseline.xml",
                        "input_sha256": "a" * 64,
                        "semantic_fingerprint": "baseline-fingerprint",
                    },
                )
            diagnostic = caught.exception.diagnostic
            self.assertEqual(diagnostic["time_hours"], 0.0)
            self.assertEqual(diagnostic["steps_completed"], 0)
            self.assertEqual(
                diagnostic["finite_or_closed_pool_concentrations_mmol_l"], {"EX_glc": 0.0}
            )
            self.assertEqual(diagnostic["depleted_finite_or_closed_pool_reaction_ids"], ["EX_glc"])
            exchange = diagnostic["dynamic_exchange_bounds"][0]
            self.assertEqual(exchange["pool_availability_cap_mmol_gdw_h"], 0.0)
            self.assertEqual(exchange["effective_uptake_cap_mmol_gdw_h"], 0.0)
            self.assertEqual(diagnostic["model"]["semantic_fingerprint"], "baseline-fingerprint")
            self.assertEqual(
                diagnostic["ordinary_fba_control"]["purpose"],
                "diagnostic_only_not_a_fallback",
            )

            context = array_context([])
            context.update({
                "medium": medium,
                "positive_gene_ids": ["g_transport"],
                "baseline": toy_model(False, requires_flux=True),
                "candidate": toy_model(False),
                "inputs": {
                    "baseline_model": {"path": "/baseline.xml", "sha256": "a" * 64},
                    "candidate_model": {"path": "/candidate.xml", "sha256": "b" * 64},
                },
                "model_contracts": {
                    "baseline": {"semantic_fingerprint": "baseline-fingerprint"},
                    "candidate": {"semantic_fingerprint": "candidate-fingerprint"},
                },
            })
            write_json(root / "summary.json", {"status": "running"})
            args = argparse.Namespace(output_dir=root)
            with (
                mock.patch.object(dfba, "_prepare_context", return_value=context),
                mock.patch.object(dfba, "_assert_context_unchanged"),
                mock.patch.object(dfba, "compare_models", side_effect=AssertionError),
            ):
                summary = dfba.run_wt_diagnostic(args)
            self.assertEqual(summary["status"], "complete")
            self.assertEqual(summary["diagnostic_outcome"], "pfba_infeasibility_observed")
            self.assertFalse(summary["execution"]["essentiality_screen_performed"])
            self.assertEqual(summary["evaluated_gene_count"], 0)
            self.assertIsNone(summary["final_scientific_gate_passed"])
            self.assertEqual(summary["wild_type"]["baseline"]["status"], "infeasible")

    def test_rescue_diagnostic_finds_the_minimum_snapshot_bound_relaxation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            medium_path = root / "medium.csv"
            write_toy_rescue_medium(medium_path)
            medium = load_dynamic_medium(medium_path)
            baseline = toy_model(False, requires_flux=True, include_closed_leucine=True)
            candidate = toy_model(False, include_closed_leucine=True)
            baseline_contract = dfba._model_contract(baseline, medium)
            candidate_contract = dfba._model_contract(candidate, medium)
            identity = {
                "role": "baseline",
                "input_path": "/baseline.xml",
                "input_sha256": "a" * 64,
                "semantic_fingerprint": baseline_contract["semantic_fingerprint"],
            }
            with self.assertRaises(dfba.DfbaInfeasibleError) as caught:
                simulate_dfba(
                    baseline, medium, hours=0.1, step_hours=0.1,
                    initial_biomass_gdw_l=0.01, model_identity=identity,
                )
            diagnostic = {**caught.exception.diagnostic, "time_hours": 4.3}
            context = array_context([])
            context.update({
                "medium": medium,
                "positive_gene_ids": ["g_transport"],
                "baseline": baseline,
                "candidate": candidate,
                "inputs": {
                    "baseline_model": {"path": "/baseline.xml", "sha256": "a" * 64},
                    "candidate_model": {"path": "/candidate.xml", "sha256": "b" * 64},
                    "dynamic_medium": {"path": str(medium_path), "sha256": sha256(medium_path)},
                    "experimental": {"path": "/reference.csv", "sha256": "c" * 64},
                    "runner_script": {"path": "/source-runner.py", "sha256": "d" * 64},
                    "essentiality_loader": {"path": "/loader.py", "sha256": "e" * 64},
                    "fingerprint_helper": {"path": "/helper.py", "sha256": "f" * 64},
                },
                "model_contracts": {
                    "baseline": baseline_contract,
                    "candidate": candidate_contract,
                },
                "dynamic_medium_contract": dfba._medium_contract(medium),
            })
            source = {
                "schema_version": 5,
                "status": "complete",
                "diagnostic_outcome": "pfba_infeasibility_observed",
                "execution": {"mode": "wt_diagnostic"},
                "git_commit": "source-commit",
                "inputs": {
                    **context["inputs"],
                    "runner_script": {
                        "path": "/old-wt-runner.py", "sha256": "0" * 64
                    },
                },
                "settings": context["settings"],
                "model_contracts": context["model_contracts"],
                "dynamic_medium_contract": context["dynamic_medium_contract"],
                "wild_type": {
                    "baseline": {
                        "status": "infeasible",
                        "infeasibility_diagnostic": diagnostic,
                    }
                },
            }
            source_path = root / "wt_summary.json"
            write_json(source_path, source)
            context["inputs"]["rescue_diagnostic_summary"] = {
                "path": str(source_path), "sha256": sha256(source_path)
            }
            write_json(root / "summary.json", {"status": "running"})
            args = argparse.Namespace(
                output_dir=root, rescue_diagnostic_summary=source_path
            )
            original_bounds = {reaction.id: reaction.bounds for reaction in baseline.reactions}
            with (
                mock.patch.object(dfba, "_prepare_context", return_value=context),
                mock.patch.object(dfba, "_assert_context_unchanged"),
            ):
                summary = dfba.run_rescue_diagnostic(args)
            self.assertEqual(summary["status"], "complete")
            self.assertEqual(summary["source_wt_diagnostic"]["source_runner_sha256"], "0" * 64)
            self.assertEqual(
                summary["diagnostic_outcome"],
                "minimum_fba_pfba_feasible_within_depleted_finite_pool_universe_found",
            )
            self.assertEqual(
                summary[
                    "minimum_fba_pfba_feasible_within_depleted_finite_pool_universe_cardinality"
                ],
                1,
            )
            self.assertEqual(
                summary[
                    "minimum_fba_pfba_feasible_within_depleted_finite_pool_universe_sets"
                ][0]["rescue_reaction_ids"],
                ["EX_glc"],
            )
            self.assertEqual(
                summary["snapshot_bound_counterfactual"]["excluded_closed_by_formulation_reaction_ids"],
                ["R1219"],
            )
            self.assertEqual(summary["snapshot_bound_counterfactual"]["no_rescue_control"]["fba"]["status"], "infeasible")
            self.assertEqual(summary["snapshot_bound_counterfactual"]["no_rescue_control"]["pfba"]["status"], "infeasible")
            with Path(summary["rescue_results"]["path"]).open(newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["rescue_reaction_ids"], "EX_glc")
            self.assertEqual(rows[0]["positive_growth_rescue"], "True")
            self.assertEqual(
                {reaction.id: reaction.bounds for reaction in baseline.reactions}, original_bounds
            )

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

    def test_gene_shards_are_deterministic_complete_and_bounded(self) -> None:
        genes = ["g3", "g1", "g2", "g4", "g5"]
        shards = [
            shard_gene_ids(genes, shard_index=index, shard_count=3)
            for index in range(3)
        ]
        self.assertEqual(shards, [["g1", "g4"], ["g2", "g5"], ["g3"]])
        self.assertEqual(sorted(gene_id for shard in shards for gene_id in shard), sorted(genes))
        self.assertEqual(shard_gene_ids(genes, shard_index=0, shard_count=1), sorted(genes))
        with self.assertRaisesRegex(ValueError, "shard_count"):
            shard_gene_ids(genes, shard_index=0, shard_count=0)
        with self.assertRaisesRegex(ValueError, "shard_index"):
            shard_gene_ids(genes, shard_index=3, shard_count=3)
        with self.assertRaisesRegex(ValueError, "shard_count"):
            shard_gene_ids(["g1"], shard_index=0, shard_count=2)
        with self.assertRaisesRegex(ValueError, "exactly True or False"):
            dfba._strict_bool(1, field="test")

    def test_aggregate_shards_defers_and_applies_the_final_gate(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "shards").mkdir()
            context = array_context(["g1", "g2", "g3"])
            write_shard_fixture(root, context, shard_index=0, shard_count=2)
            write_shard_fixture(
                root,
                context,
                shard_index=1,
                shard_count=2,
                new_false_negative_ids={"g2"},
            )
            args = argparse.Namespace(
                output_dir=root,
                aggregate_shards=root / "shards",
                shard_count=2,
            )
            with (
                mock.patch.object(dfba, "_prepare_context", return_value=context),
                mock.patch.object(dfba, "_assert_context_unchanged"),
            ):
                summary = aggregate_shards(args)
            self.assertEqual(summary["evaluated_gene_count"], 3)
            self.assertEqual(summary["new_false_negative_gene_ids"], ["g2"])
            self.assertFalse(summary["final_scientific_gate_passed"])
            with (root / "gene_results.tsv").open(newline="", encoding="utf-8") as stream:
                self.assertEqual(
                    [row["gene_id"] for row in csv.DictReader(stream, delimiter="\t")],
                    ["g1", "g2", "g3"],
                )

    def test_aggregate_shards_rejects_missing_or_wrong_gene_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "shards").mkdir()
            context = array_context(["g1", "g2"])
            write_shard_fixture(root, context, shard_index=0, shard_count=2)
            args = argparse.Namespace(
                output_dir=root,
                aggregate_shards=root / "shards",
                shard_count=2,
            )
            with (
                mock.patch.object(dfba, "_prepare_context", return_value=context),
                mock.patch.object(dfba, "_assert_context_unchanged"),
            ):
                with self.assertRaises(FileNotFoundError):
                    aggregate_shards(args)

            write_shard_fixture(root, context, shard_index=1, shard_count=2)
            path = root / "shards" / "001" / "gene_results.tsv"
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            rows[0]["gene_id"] = "g1"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=RESULT_COLUMNS, delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            with (
                mock.patch.object(dfba, "_prepare_context", return_value=context),
                mock.patch.object(dfba, "_assert_context_unchanged"),
            ):
                with self.assertRaises(ValueError):
                    aggregate_shards(args)

            (root / "shards" / "stale").mkdir()
            with (
                mock.patch.object(dfba, "_prepare_context", return_value=context),
                mock.patch.object(dfba, "_assert_context_unchanged"),
            ):
                with self.assertRaisesRegex(ValueError, "exactly the expected"):
                    aggregate_shards(args)


if __name__ == "__main__":
    unittest.main()
