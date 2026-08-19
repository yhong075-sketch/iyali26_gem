"""Acceptance checks for the read-only strict sn-position lipid-core candidate."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import unittest
from math import fsum, isclose
from pathlib import Path

import cobra
from cobra.flux_analysis import flux_variability_analysis, pfba


REPOSITORY = Path(__file__).resolve().parents[1]
MODEL_PATH = Path(os.environ.get(
    "IYALI26_SOURCE_MODEL", REPOSITORY.parent / "iyali26_gem" / "model.xml"
))
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from scripts.lp_sn12_candidate import (  # noqa: E402
    CARDIOLIPIN_AUDIT_SHA256, CURATION_PATH, ContractError, MARKER,
    build_candidate, cardiolipin_audit, main, report, source_fingerprint,
)


class StrictSnCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.curation = json.loads(CURATION_PATH.read_text())
        cls.source = cobra.io.read_sbml_model(str(MODEL_PATH))
        cls.source_fingerprint = source_fingerprint(cls.source)
        cls.source_objective = cls.source.slim_optimize(error_value=None)
        cls.candidate = build_candidate(cls.source)
        cls.generated = [reaction for reaction in cls.candidate.reactions if reaction.annotation.get(MARKER) == "true"]

    def test_source_contract_is_frozen_and_input_unchanged(self) -> None:
        self.assertEqual(hashlib.sha256(MODEL_PATH.read_bytes()).hexdigest(), self.curation["source"]["model_sha256"])
        self.assertEqual(self.source_fingerprint, self.curation["source"]["model_fingerprint"])
        self.assertTrue(isclose(self.source.slim_optimize(error_value=None), self.source_objective, rel_tol=0.0, abs_tol=1e-9))
        for mutate in (lambda model: setattr(model.metabolites[0], "name", "drift"), lambda model: setattr(model.reactions[0], "annotation", {"drift": "true"})):
            with self.subTest(mutate=mutate):
                drifted = self.source.copy()
                mutate(drifted)
                with self.assertRaises(ContractError):
                    build_candidate(drifted)
        with self.assertRaises(ContractError):
            build_candidate(self.candidate)

    def test_structure_balance_and_template_copying(self) -> None:
        self.assertEqual(len(self.generated), 1134)
        self.assertEqual(source_fingerprint(build_candidate(self.source)), source_fingerprint(self.candidate))
        self.assertEqual(sum(reaction.id.startswith("UL_R353_") for reaction in self.generated), 7)
        self.assertEqual(sum(reaction.id.startswith("UL_R1771_") for reaction in self.generated), 343)
        self.assertEqual(len({met.id for reaction in self.generated for met in reaction.metabolites if met.annotation.get(MARKER) == "true"}), 1134)
        self.assertFalse(any(metabolite_id in self.candidate.metabolites for metabolite_id in self.curation["generic_acyl_coa_ids"]))
        self.assertFalse(any(reaction.id.startswith("xREMIX_") for reaction in self.candidate.reactions))
        for reaction in self.generated:
            template = self.source.reactions.get_by_id(reaction.annotation["template_reaction"])
            self.assertEqual(reaction.bounds, template.bounds)
            self.assertEqual(reaction.gene_reaction_rule, template.gene_reaction_rule)
            self.assertTrue(all(reaction.annotation.get(key) == value for key, value in template.annotation.items()))
            self.assertEqual(reaction.check_mass_balance(), {})

    def test_lipid_exchange_metadata_is_explicit_and_gpr_free(self) -> None:
        exchange_ids = {"R1776", "R1778", "R1786", "R1787", "R1790", "R1792", "R1793"}
        classification = "composite_or_nonenzymatic_lipid_exchange"
        steps = {step["reaction_id"]: step for step in self.curation["steps"]}
        self.assertEqual(self.candidate.compartments["C_em"], "ER membrane")
        self.assertEqual({reaction_id for reaction_id, step in steps.items() if step.get("reaction_classification") == classification}, exchange_ids)
        self.assertTrue(all(steps[reaction_id]["gpr"] == "" for reaction_id in exchange_ids))
        self.assertEqual(
            {key: steps["R1776"][key] for key in ("source_state", "target_state", "template_source_id", "template_target_id")},
            {"source_state": "dag_cy", "target_state": "dag_em", "template_source_id": "m1670[C_cy]", "template_target_id": "m1649[C_em]"},
        )
        self.assertEqual(steps["R1776"]["construction_direction"], "reverse")
        for reaction in self.generated:
            template_id = reaction.annotation["template_reaction"]
            if template_id in exchange_ids:
                self.assertEqual(reaction.annotation.get("reaction_classification"), classification)

    def test_cli_refuses_source_overwrite(self) -> None:
        original_arguments = sys.argv[:]
        self.addCleanup(setattr, sys, "argv", original_arguments)
        sys.argv = ["lp_sn12_candidate.py", str(MODEL_PATH), "--report", "candidate.json", "--candidate-sbml", str(MODEL_PATH)]
        with self.assertRaisesRegex(SystemExit, "cannot overwrite"):
            main()

    def test_tuple_order_and_biomass_weights(self) -> None:
        self.assertIn("UL_R1846_SN1__lauroyl__SN2__myristoyl", self.candidate.reactions)
        self.assertIn("UL_R1846_SN1__myristoyl__SN2__lauroyl", self.candidate.reactions)
        tag = self.candidate.reactions.get_by_id("UL_R1771_SN1__lauroyl__SN2__myristoyl__SN3__oleoyl")
        self.assertIn("ul_dag_lp__sn1__lauroyl__sn2__myristoyl[C_lp]", [met.id for met in tag.metabolites])
        self.assertIn("m1490[C_lp]", [met.id for met in tag.metabolites])
        biomass = self.candidate.reactions.get_by_id("biomass_C")
        for row in self.curation["biomass"]["terms"]:
            terms = [(-coefficient, metabolite) for metabolite, coefficient in biomass.metabolites.items() if metabolite.annotation.get("state") == row["state"]]
            self.assertEqual(len(terms), 7 ** row["arity"])
            self.assertTrue(isclose(fsum(value for value, _ in terms), row["coefficient"], abs_tol=1e-18))

    def test_cardiolipin_audit_is_empty_fail_closed_and_stereospecific(self) -> None:
        before = source_fingerprint(self.source)
        audit = cardiolipin_audit(self.candidate)
        self.assertEqual(cardiolipin_audit(self.source)["state_space"], audit["state_space"])
        contract = self.curation["cardiolipin_four_chain_audit"]
        space = audit["state_space"]
        envelope = space["not_adopted_modeling_envelope"]
        self.assertEqual(hashlib.sha256(json.dumps(contract, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest(), CARDIOLIPIN_AUDIT_SHA256)
        self.assertEqual(audit["contract_sha256"], CARDIOLIPIN_AUDIT_SHA256)
        self.assertEqual(audit["mode"], "audit_only")
        genes = contract["gene_identity_ledger"]
        evidence_classes = {"experimental_support", "database_inference", "model_assignment", "unverified"}
        evidence_categories = {"experimentally verified", "curated annotation", "model/GPR assignment only", "uncharacterized"}
        self.assertEqual(len(genes), 9)
        self.assertTrue(all(row.get("legacy_id") and any(identifier.startswith("GeneID:") for identifier in row["database_ids"]) and row.get("established_name") and row.get("protein_function") and "model_gpr_role" in row and row["evidence_category"] in evidence_categories and row["evidence_status"] and set(row["evidence_status"]) <= evidence_classes and row.get("evidence_detail") for row in genes))
        ledger_ids = {row["systematic_id"] for row in genes}
        self.assertEqual(len(ledger_ids), len(genes))
        self.assertLessEqual({row["gpr"] for row in contract["template_reactions"] if row["gpr"]}, ledger_ids)
        self.assertFalse(any("source-model-assigned" in row["protein_function"] for row in genes))
        self.assertEqual(contract["R197_gpr_decision_contract"]["current_gpr"], "YALI1E17431g")
        self.assertEqual(contract["R197_gpr_decision_contract"]["mitochondrial_R197_evidence_status"], "unverified")
        self.assertEqual((contract["R197_gpr_decision_contract"]["replacement_candidate"], contract["R197_gpr_decision_contract"]["replacement_action"]), ("YALI1C17238g", "forbidden"))
        pgc1 = next(row for row in genes if row["systematic_id"] == "YALI1F11605g")
        self.assertEqual((pgc1["candidate_label"], pgc1["database_ids"], pgc1["deleted_database_ids"]), ("PGC1", ["GeneID:2908510", "UniProt:Q6C2F9"], ["UniProt:A0A1D8NMM1"]))
        self.assertEqual(pgc1["old_aliases"], ["YALI2_F00162g", "RefSeq:XP_505153.1"])
        lclat = next(row for row in genes if row["systematic_id"] == "YALI1B16087g")
        self.assertEqual((lclat["legacy_id"], lclat["old_aliases"], lclat["normalized_aliases"], lclat["database_ids"]), ("YALI0B12254g", ["YALI2_C00842g"], ["YALI2C00842g"], ["GeneID:2907165"]))
        self.assertEqual((lclat["candidate_label"], lclat["protein_function"], lclat["model_gpr_role"], lclat["gpr_replacement_or_addition_action"]), ("LCLAT1-like", "lysocardiolipin acyltransferase candidate", None, "forbidden"))
        self.assertEqual(contract["R1742_annotation_conflict"]["model_ec_match_verdict"], "contradicted")
        self.assertEqual((contract["R1742_annotation_conflict"]["kegg_supported_ec_scope"], contract["R1742_annotation_conflict"]["reaction_expansion_decision"]), ("3.1.4.-", "blocked"))
        sources = contract["evidence_sources"]
        self.assertEqual(len(sources), 16)
        self.assertTrue(all(row["accessed_on"] == "2026-08-18" and row.get("url") and row.get("supports") and row.get("reviewer_verdict") and row["audit_status"] == "audited" and row["decision"] in {"use_as_constraint", "defer", "exclude"} and row.get("conditions") for row in sources))
        conflict = contract["R1742_annotation_conflict"]
        claim_ids = [claim_id for row in sources for claim_id in row["claim_ids"]] + [conflict["claim_id"]]
        self.assertEqual(len(claim_ids), len(set(claim_ids)))
        self.assertEqual(len(claim_ids), contract["evidence_audit_coverage"]["total_claims"])
        verdicts = [row["reviewer_verdict"] for row in sources for _ in row["claim_ids"]] + [conflict["reviewer_verdict"]]
        self.assertLessEqual(set(verdicts), {"supported", "partially_supported", "contradicted"})
        self.assertEqual((verdicts.count("supported"), verdicts.count("partially_supported"), verdicts.count("contradicted")), (10, 6, 1))
        self.assertEqual((conflict["audit_status"], conflict["decision"], set(conflict["source_ids"])), ("audited", "exclude", {"R1742_NCBI_UNIPROT_KEGG", "IUBMB_EC_3_1_4_2"}))
        self.assertEqual(contract["evidence_audit_coverage"], {"total_claims": 17, "audited": 17, "supported": 10, "unresolved": 6, "contradicted": 1, "unchecked": 0})
        self.assertEqual(sum(contract["evidence_audit_coverage"][key] for key in ("supported", "unresolved", "contradicted")), contract["evidence_audit_coverage"]["total_claims"])
        self.assertEqual(audit["gene_identity_ledger"], contract["gene_identity_ledger"])
        self.assertEqual(audit["evidence_sources"], contract["evidence_sources"])
        self.assertEqual(audit["evidence_audit_coverage"], contract["evidence_audit_coverage"])
        self.assertEqual((audit["network_evidence"]["R1723_removed_position_and_chain"], audit["network_evidence"]["R1724_MLCL_vacancy_PC_donor_site_and_chain"], audit["network_evidence"]["CL_MLCL_leaflet_redistribution"]), ("unresolved", "unresolved", "unknown"))
        self.assertEqual(source_fingerprint(self.source), before)
        self.assertEqual(len(self.generated), 1134)
        self.assertEqual(space["digest_convention"], "UTF-8; Python/Unicode ascending sort; LF join; trailing LF")
        self.assertEqual(space["ordered_half_pair_count"], 49)
        self.assertEqual(space["ordered_half_labels_sha256"], "2b7fc64ec0a8a6e0ab46bedca36e7bf48fff46b35032418855d5218a3d937dd5")
        self.assertEqual(
            (envelope["status"], envelope["all_state_count"], envelope["same_half_state_count"], envelope["different_half_state_count"]),
            ("not_adopted", 2401, 49, 2352),
        )
        self.assertEqual(envelope["all_ids_sha256"], "af34ecde5ef4a2823208c2ff7e81733b42d907034ca816fdafb9f9f1288f0ecb")
        self.assertEqual(envelope["same_half_ids_sha256"], "3eb9cba676499fa5daa16077fd36282f22988970aeaa12c2e5f357fec7bb39b6")
        self.assertEqual(envelope["different_half_ids_sha256"], "3db29ba9826a80e5d9e777c98c1b2c635427b6dccf3817402629031601b5fcd9")
        forbidden = space["forbidden_lossy_half_swap_fold"]
        self.assertEqual(forbidden, {"count_formula": "49*50/2", "state_count": 1225, "status": "forbidden"})
        self.assertNotIn("ids_sha256", forbidden)
        self.assertTrue(audit["half_swap_witness"]["ids_are_distinct"])
        self.assertNotEqual(audit["half_swap_witness"]["original_id"], audit["half_swap_witness"]["swapped_id"])
        self.assertEqual({value for key, value in space["current_strict_model"].items() if key.endswith("_count")}, {0})
        self.assertIn("does not prove all 2,401", space["current_strict_model"]["reachability_scope"])
        self.assertFalse(any(audit["gates"].values()))
        self.assertEqual(audit["activation"], {"state": "blocked", "ready_for_activation": False, "production_apply_allowed": False, "human_approval_recorded": False, "human_gate_required": True})
        self.assertEqual(audit["biomass_contract"]["positioned_biomass_term_count"], 0)
        self.assertFalse(audit["biomass_contract"]["xLIPID_weight_adopted"])
        self.assertEqual(audit["biomass_contract"]["xLIPID_generic_cardiolipin_coefficient"], -0.842)
        profile = audit["biomass_contract"]["sekova_profile_fraction"]
        self.assertEqual((profile["reported_percent_range"], profile["quantity_type"]), ([7, 9], "relative cardiolipin profile fraction"))
        self.assertFalse(profile["absolute_mol_per_gDW_available"] or profile["species_response_factor_calibration_available"] or profile["sn_position_weights_available"])
        self.assertEqual((audit["biomass_contract"]["cardiolipin_biomass_weight_decision"], audit["biomass_contract"]["independent_chain_product_weighting"]), ("blocked", "forbidden"))
        self.assertFalse(any(met.id.startswith("ul_cl_mm__") for met in self.candidate.metabolites))
        self.assertFalse(any(any(met.id.startswith("ul_cl_mm__") for met in reaction.metabolites) for reaction in self.candidate.reactions))
        for row in contract["template_reactions"]:
            reaction = self.source.reactions.get_by_id(row["reaction_id"])
            self.assertEqual(list(reaction.bounds), row["bounds"])
            self.assertEqual(reaction.gene_reaction_rule, row["gpr"])
            self.assertEqual(sorted(reaction.compartments), row["compartments"])
            self.assertEqual({met.id: value for met, value in reaction.metabolites.items()}, row["stoichiometry"])
            if "notes" in row:
                self.assertEqual(reaction.notes, row["notes"])
        drifted = self.source.copy()
        drifted.reactions.get_by_id("R1721").bounds = (0.0, 999.0)
        with self.assertRaises(ContractError):
            cardiolipin_audit(drifted)
        notes_drifted = self.source.copy()
        notes_drifted.reactions.get_by_id("R1742").notes = {"PROTEIN_CLASS": "3.1.4.-"}
        with self.assertRaises(ContractError):
            cardiolipin_audit(notes_drifted)
        positioned = self.candidate.copy()
        positioned.add_metabolites([cobra.Metabolite("ul_cl_mm__proS_sn1__lauroyl__proS_sn2__lauroyl__proR_sn1__lauroyl__proR_sn2__lauroyl[C_mm]", compartment="C_mm")])
        with self.assertRaises(ContractError):
            cardiolipin_audit(positioned)

    def test_sbml_round_trip_and_function(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.xml"
            cobra.io.write_sbml_model(self.candidate, str(path))
            reread = cobra.io.read_sbml_model(str(path))
        self.assertEqual(len([reaction for reaction in reread.reactions if reaction.annotation.get(MARKER) == "true"]), 1134)
        tag = reread.reactions.get_by_id("UL_R1771_SN1__lauroyl__SN2__myristoyl__SN3__oleoyl")
        self.assertEqual(tag.gene_reaction_rule, self.source.reactions.get_by_id("R1771").gene_reaction_rule)
        self.assertEqual(tag.bounds, self.source.reactions.get_by_id("R1771").bounds)
        self.assertGreater(self.candidate.slim_optimize(error_value=None), 0)
        self.assertGreater(pfba(self.candidate).fluxes["biomass_C"], 0)
        payload = report(self.source, self.candidate)
        self.assertFalse(payload["ready_for_activation"])
        self.assertEqual(payload["cardiolipin_four_chain_audit"]["activation"]["state"], "blocked")
        runtime = payload["runtime_seconds"]
        self.assertEqual(
            payload["performance_review_required"],
            any(
                runtime[f"candidate_{kind}_median"] > 2 * runtime[f"source_{kind}_median"]
                for kind in ("fba", "pfba")
            ),
        )
        self.assertGreaterEqual(payload["fba_probe"]["growth_ratio"], 0.99)
        self.assertLessEqual(payload["fba_probe"]["growth_ratio"], 1.0)

    def test_required_lipid_routes_are_not_bypassed(self) -> None:
        route_sets = {
            "pa_transport": ("UL_R1790_",), "pi": ("UL_R1839_",),
            "ps_pe_pc": ("UL_R1841_", "UL_R629_", "UL_R1835_", "UL_R1845_", "UL_R1844_"), "tag": ("UL_R1771_",),
        }
        for name, prefixes in route_sets.items():
            with self.candidate:
                for reaction in self.candidate.reactions:
                    if reaction.id.startswith(prefixes):
                        reaction.bounds = (0.0, 0.0)
                self.assertEqual(self.candidate.slim_optimize(error_value=None), 0.0, name)
        closed = self.curation["close_reactions"] + ["newBiom", "xBIOMASS", "xLIPID"]
        variability = flux_variability_analysis(self.candidate, reaction_list=closed, processes=1)
        self.assertTrue((variability["minimum"] == 0).all() and (variability["maximum"] == 0).all())

    def test_each_required_chain_is_not_bypassed(self) -> None:
        for chain in self.curation["chains"]:
            with self.candidate:
                for reaction in self.candidate.reactions:
                    if chain["lp_acyl_coa_id"] in (met.id for met in reaction.metabolites):
                        reaction.bounds = (0.0, 0.0)
                self.assertEqual(self.candidate.slim_optimize(error_value=None), 0.0, chain["id"])


if __name__ == "__main__":
    unittest.main()
