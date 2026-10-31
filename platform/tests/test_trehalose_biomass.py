"""Trehalose biomass removal: preservation, rejection, SBML roundtrip, context exit and opt-in flag; no LPs."""

import copy
import json
import re
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.trehalose_biomass import BIOMASS, CURATION_PATH, TERM, apply_trehalose_biomass
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_"
            "gpr_r673_prolyl_split_aco_split_20261029.xml")


def state(model):
    return semantics(model), annotations(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class TrehaloseBiomassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        composition = "composition differs from the reviewed baseline"
        expected = {"other_coefficient": composition, "extra_entry": composition, "bounds": composition,
                    "coefficient": re.escape(f"Conflicting {TERM} biomass coefficient"),
                    "identity": "identity differs", "context": "reaction context differs",
                    "notes": "Conflicting", "extra_note": "notes outside this curation",
                    "notes_with_term": "carries this curation's notes with the",
                    "removed_without_notes": "without this curation",
                    "objective": "objective is not", "direction": "objective is not", "missing": "missing"}
        baseline = state(self.model)
        notes = self.spec["biomass_term"]["notes"]
        for field, message in expected.items():
            bad = self.model.copy()
            biomass, term = bad.reactions.get_by_id(BIOMASS), bad.metabolites.get_by_id(TERM)
            if field == "other_coefficient":
                biomass.add_metabolites({bad.metabolites.get_by_id("m1324[C_cy]"): -0.1})
            elif field == "extra_entry":
                biomass.add_metabolites({bad.metabolites.get_by_id("m297[C_cy]"): -0.001})
            elif field == "bounds":
                biomass.bounds = (0.0, 100.0)
            elif field == "coefficient":
                biomass.add_metabolites({term: -0.001})
            elif field == "identity":
                term.formula = "C12H22O12"
            elif field == "context":
                bad.reactions.R1013.add_metabolites({term: -1})
            elif field == "notes":
                biomass.notes = {**biomass.notes, "trehalose_biomass_status": "another review"}
            elif field == "extra_note":
                biomass.notes = {**biomass.notes, "unreviewed_note": "added elsewhere"}
            elif field == "notes_with_term":
                biomass.notes = {**biomass.notes, **copy.deepcopy(notes)}
            elif field == "removed_without_notes":
                biomass.add_metabolites({term: 0.0}, combine=False)
            elif field == "objective":
                bad.objective = "xBIOMASS"
            elif field == "direction":
                bad.objective.direction = "min"
            else:
                bad.remove_reactions([biomass])
            bad.solver.update()  # optlang refreshes the objective lazily; flush before the snapshot
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_trehalose_biomass(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for edit in ("status", "curation_id", "coefficient_before", "reactions_after", "metabolite", "notes"):
            spec = copy.deepcopy(self.spec)
            if edit == "status":
                spec["status"] = "user_authorized"
            elif edit == "curation_id":
                spec["curation_id"] = "TREHALOSE-BIOMASS"
            elif edit == "coefficient_before":
                spec["biomass_term"]["coefficient_before"] = -0.001
            elif edit == "reactions_after":
                spec["biomass_term"]["reactions_after"] = spec["biomass_term"]["reactions_before"]
            elif edit == "metabolite":
                spec["biomass_term"]["metabolite_id"] = "m1123[C_cy]"
            else:
                spec["biomass_term"]["notes"] = {}
            with self.subTest(edit=edit), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_trehalose_biomass(self.model.copy(), spec)

    def test_changes_only_the_term_and_notes_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        target = expected.reactions.get_by_id(BIOMASS)
        target.add_metabolites({expected.metabolites.get_by_id(TERM): 0.0}, combine=False)
        target.notes = {**target.notes, **self.spec["biomass_term"]["notes"]}
        result = apply_trehalose_biomass(model)
        self.assertEqual(result["status"], "applied")
        self.assertEqual(result["coefficient_before"], -0.00686)
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertNotIn(TERM, {m.id for m in model.reactions.get_by_id(BIOMASS).metabolites})
        self.assertEqual(sorted(r.id for r in model.metabolites.get_by_id(TERM).reactions),
                         self.spec["biomass_term"]["reactions_after"])
        self.assertEqual(apply_trehalose_biomass(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "trehalose.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_trehalose_biomass(loaded)["status"], "already_correct")
        # On an applied model the composition guard still holds.
        bad = loaded.copy()
        bad.reactions.get_by_id(BIOMASS).add_metabolites({bad.metabolites.get_by_id("m1324[C_cy]"): -0.1})
        before = state(bad)
        with self.assertRaisesRegex(ValueError, "composition differs"):
            apply_trehalose_biomass(bad)
        self.assertEqual(state(bad), before)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_trehalose_biomass(model)["status"], "applied")
            self.assertNotIn(TERM, {m.id for m in model.reactions.get_by_id(BIOMASS).metabolites})
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_trehalose_biomass(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).trehalose_biomass_removal)
        self.assertTrue(parse_args(["--trehalose-biomass-removal"]).trehalose_biomass_removal)


if __name__ == "__main__":
    unittest.main()
