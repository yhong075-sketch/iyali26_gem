"""TPP requirement: preservation, rejection, SBML roundtrip, exact edits and the flag requirement; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.tpp_requirement import CURATION_PATH, apply_tpp_requirement
from tests.test_coq9_curation import annotations, semantics
from tests.test_gpi_chain_links import group_table

BASELINE = "E5_coq9_resp_pkg_r132_combined_gpi_coproducts_20261021.xml"
TARGETS = ("R1342", "biomass_C")
TPP = "m975[C_cy]"


def state(model):
    return semantics(model), annotations(model), group_table(model)


class TppRequirementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        precondition = "precondition differs"
        expected = {"species_identity": "m975\\[C_cy\\] identity differs",
                    "tpp_context": "m975\\[C_cy\\] reaction context differs",
                    "thi80_route": "R752 differs from the reviewed TPP context",
                    "mito_route": "R933 differs from the reviewed TPP context",
                    "mito_source": "m1116\\[C_mi\\] reaction context differs",
                    "objective": "biomass_C is not the sole objective",
                    "biomass_stoichiometry": f"biomass_C {precondition}",
                    "biomass_other_coefficient": f"biomass_C {precondition}",
                    "bounds": f"R1342 {precondition}", "gpr": f"R1342 {precondition}",
                    "annotation": f"R1342 {precondition}",
                    "extra_note": "R1342 notes outside this curation differ",
                    "conflicting_note": "Conflicting biomass_C TPP requirement notes",
                    "after_without_notes": "R1342 already carries the TPP change without this curation",
                    "notes_without_change": "R1342 carries this curation's notes without the change",
                    "partial": "Only one TPP requirement change is present"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            if field == "species_identity":
                bad.metabolites.get_by_id(TPP).formula = "C12H16N4O7P2S"
            elif field == "tpp_context":
                # A second TPP source would be a route the knockout review did not cover.
                source = Reaction("R_TEST_SOURCE", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([source])
                source.add_metabolites({bad.metabolites.get_by_id(TPP): 1})
            elif field == "thi80_route":
                bad.reactions.R752.bounds = (0.0, 0.0)
            elif field == "mito_route":
                # Reversible, R933 could return mitochondrial TPP and bypass THI80.
                bad.reactions.R933.bounds = (-1000.0, 1000.0)
            elif field == "mito_source":
                source = Reaction("R_TEST_MITO_SOURCE", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([source])
                source.add_metabolites({bad.metabolites.get_by_id("m1116[C_mi]"): 1})
            elif field == "objective":
                bad.objective = "R752"
            elif field == "biomass_stoichiometry":
                bad.reactions.biomass_C.add_metabolites({bad.metabolites.get_by_id("m1000[C_cy]"): -0.001})
            elif field == "biomass_other_coefficient":
                bad.reactions.biomass_C.add_metabolites({bad.metabolites.get_by_id(TPP): -1e-5})
            elif field == "bounds":
                bad.reactions.R1342.bounds = (0.0, 500.0)
            elif field == "gpr":
                bad.reactions.R1342.gene_reaction_rule = "YALI1E25443g"
            elif field == "annotation":
                bad.reactions.R1342.annotation = {**bad.reactions.R1342.annotation, "ec-code": "2.7.4.16"}
            elif field == "extra_note":
                bad.reactions.R1342.notes = {"unreviewed_note": "added elsewhere"}
            elif field == "conflicting_note":
                bad.reactions.biomass_C.notes = {**bad.reactions.biomass_C.notes,
                                                 "tpp_requirement_status": "another review"}
            elif field == "after_without_notes":
                bad.reactions.R1342.bounds = (0.0, 0.0)
            elif field == "notes_without_change":
                apply_tpp_requirement(bad)
                bad.reactions.R1342.bounds = (0.0, 1000.0)
            else:
                apply_tpp_requirement(bad)
                biomass = bad.reactions.biomass_C
                biomass.subtract_metabolites({bad.metabolites.get_by_id(TPP): biomass.get_coefficient(TPP)})
                biomass.notes = dict(self.spec["targets"]["biomass_C"]["notes_before"])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_tpp_requirement(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for key, value in (("status", "user_authorized"), ("curation_id", "TPP-REQUIREMENT"),
                           ("targets", {"R1342": self.spec["targets"]["R1342"]}),
                           ("requires_reactions", {}), ("species", {})):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_tpp_requirement(self.model.copy(), spec)

    def test_spec_values_outside_the_reviewed_change_are_rejected_without_edits(self):
        def r1342_open(spec):
            spec["targets"]["R1342"]["bounds_after"] = [0.0, 500.0]

        def r1342_invalid(spec):
            # Would fail inside cobra after biomass_C was already edited, if it got that far.
            spec["targets"]["R1342"]["bounds_after"] = [10.0, 0.0]

        def biomass_extra_term(spec):
            spec["targets"]["biomass_C"]["after_stoichiometry"]["m1000[C_cy]"] = -0.001

        def biomass_bounds(spec):
            spec["targets"]["biomass_C"]["bounds_after"] = [0.0, 10.0]

        def tpp_positive(spec):
            spec["targets"]["biomass_C"]["after_stoichiometry"][TPP] = 1.2e-6

        model = self.model.copy()
        before = state(model)
        for edit in (r1342_open, r1342_invalid, biomass_extra_term, biomass_bounds, tpp_positive):
            spec = copy.deepcopy(self.spec)
            edit(spec)
            with self.subTest(edit=edit.__name__), self.assertRaisesRegex(ValueError, "spec values differ"):
                apply_tpp_requirement(model, spec)
            self.assertEqual(state(model), before)

    def test_changes_only_the_two_targets_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected.reactions.biomass_C.add_metabolites({expected.metabolites.get_by_id(TPP): -1.2e-6})
        expected.reactions.R1342.bounds = (0.0, 0.0)
        for rid in TARGETS:
            reaction, row = expected.reactions.get_by_id(rid), self.spec["targets"][rid]
            reaction.notes = {**reaction.notes, **{
                "tpp_requirement_curation_id": self.spec["curation_id"],
                "tpp_requirement_status": self.spec["status"],
                "tpp_requirement_authorization": self.spec["authorization"],
                "tpp_requirement_evidence": row["evidence"],
                "tpp_requirement_evidence_tier": row["evidence_tier"],
                "tpp_requirement_rationale": row["rationale"],
                "tpp_requirement_previous_state": row["previous_state"]}}
        result = apply_tpp_requirement(model)
        self.assertEqual((result["status"], result["changed_reactions"]), ("applied", list(TARGETS)))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual(model.reactions.biomass_C.get_coefficient(TPP), -1.2e-6)
        self.assertEqual(model.reactions.R1342.bounds, (0.0, 0.0))
        # Closed, not deleted: the equation and GPR stay for review.
        self.assertEqual({m.id: c for m, c in model.reactions.R1342.metabolites.items()},
                         self.spec["targets"]["R1342"]["before_stoichiometry"])
        self.assertEqual(model.reactions.R752.gene_reaction_rule, "YALI1E25443g")
        self.assertEqual(apply_tpp_requirement(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tpp_requirement.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_tpp_requirement(loaded)["status"], "already_correct")

    def test_build_time_notes_are_accepted(self):
        # In a build the biomass_C notes hold -> and "; read back from the SBML file they are escaped.
        model = self.model.copy()
        model.reactions.biomass_C.notes = dict(self.spec["targets"]["biomass_C"]["notes_before"])
        self.assertNotEqual(model.reactions.biomass_C.notes, self.model.reactions.biomass_C.notes)
        self.assertEqual(apply_tpp_requirement(model)["status"], "applied")
        self.assertEqual(apply_tpp_requirement(model)["status"], "already_correct")

    def test_flag_requires_the_gpi_coproducts(self):
        self.assertFalse(parse_args([]).tpp_requirement)
        self.assertTrue(parse_args(["--tpp-requirement"]).tpp_requirement)
        with self.assertRaisesRegex(ValueError, "defined on the GPI co-products candidate"):
            build_reference_chain(tpp_requirement=True)


if __name__ == "__main__":
    unittest.main()
