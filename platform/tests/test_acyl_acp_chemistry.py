"""In-memory checks for the 2026-10-20 opt-in acyl-ACP chemistry correction; no optimization."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.acyl_acp_chemistry import SPEC_PATH, apply_acyl_acp_chemistry
from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.lipoylation_module import apply_lipoylation_module
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model

SPEC = json.loads(SPEC_PATH.read_text())
# fdd00ac3 (164f235f + the lipoylation module) is not on this branch; it is rebuilt in memory and
# round-tripped through SBML so the step sees a saved model, as it did on fdd00ac3.
BASELINE = "E5_coq9_resp_pkg_r132_combined_20261020.xml"


def snapshot(model):
    return (model_definition(model), {r.id: dict(r.notes) for r in model.reactions})


class AcylAcpChemistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        apply_lipoylation_module(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lipoylation.xml"
            write_deterministic_sbml_model(model, path)
            cls.baseline = read_sbml_model(str(path))

    def fresh(self):
        return self.baseline.copy()

    def test_changes_only_the_block_and_balances_elements(self):
        model = self.fresh()
        before = model_definition(model)
        self.assertEqual(apply_acyl_acp_chemistry(model)["status"], "applied")
        after = model_definition(model)
        changed_reactions = sorted(r for r in after["reactions"] if after["reactions"][r] != before["reactions"][r])
        self.assertEqual(changed_reactions, sorted(SPEC["changed_reactions"]))
        changed_species = sorted(m for m in after["metabolites"] if after["metabolites"][m] != before["metabolites"][m])
        self.assertEqual(changed_species, sorted(SPEC["formula_fills"]))
        for rid in SPEC["block_reactions"]:
            residual = model.reactions.get_by_id(rid).check_mass_balance()
            self.assertEqual({k for k in residual if abs(residual[k]) > 1e-9} - {"charge"}, set(), rid)
        for rid in SPEC["made_irreversible"]:
            self.assertEqual(model.reactions.get_by_id(rid).lower_bound, 0.0)

    def test_rerun_and_lipoylation_rerun_are_no_ops_after_sbml_round_trip(self):
        model = self.fresh()
        apply_acyl_acp_chemistry(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "acyl_acp.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        state = snapshot(loaded)
        self.assertEqual(apply_acyl_acp_chemistry(loaded)["status"], "already_correct")
        self.assertEqual(apply_lipoylation_module(loaded)["status"], "already_correct")  # P5
        self.assertEqual(snapshot(loaded), state)
        stripped = loaded.copy()
        stripped.reactions.R2226.notes = {}
        with self.assertRaises(ValueError):
            apply_acyl_acp_chemistry(stripped)
        # P5 negatives: the lipoylation guard must not accept a corrected residual without this
        # change's note, nor a note with another value.
        no_note = loaded.copy()
        del no_note.reactions.R2230.notes["acyl_acp_chemistry"]
        other_note = loaded.copy()
        other_note.reactions.R2230.notes["acyl_acp_chemistry"] = "OTHER"
        for case in (no_note, other_note):
            with self.assertRaises(ValueError):
                apply_lipoylation_module(case)

    def test_stale_or_partial_inputs_fail_closed_without_edits(self):
        def one_fill_only(m): m.metabolites.get_by_id("m1411[C_mi]").formula = "C4H7OSR"
        def hydrolase_bound(m): m.reactions.R2226.lower_bound = 0.0
        def extra_user(m):
            extra = Reaction("EXTRA"); m.add_reactions([extra])
            extra.add_metabolites({m.metabolites.get_by_id("m1920[C_cy]"): -1.0})
        def context_species(m): m.metabolites.get_by_id("m801[C_mi]").formula = "C3H4O3SR"
        def no_lipoylation_acp1(m): m.metabolites.get_by_id("m179[C_mi]").formula = None
        for change in (one_fill_only, hydrolase_bound, extra_user, context_species, no_lipoylation_acp1):
            with self.subTest(change.__name__):
                model = self.fresh()
                change(model)
                state = snapshot(model)
                with self.assertRaises(ValueError):
                    apply_acyl_acp_chemistry(model)
                self.assertEqual(snapshot(model), state)

    def test_unexpected_scope_is_rejected(self):
        spec = copy.deepcopy(SPEC)
        spec["change_id"] = "OTHER"
        model = self.fresh()
        state = snapshot(model)
        with self.assertRaises(ValueError):
            apply_acyl_acp_chemistry(model, spec)
        self.assertEqual(snapshot(model), state)

    def test_flag_is_opt_in_and_requires_the_lipoylation_module(self):
        self.assertFalse(parse_args([]).acyl_acp_chemistry)
        self.assertTrue(parse_args(["--acyl-acp-chemistry"]).acyl_acp_chemistry)
        with self.assertRaisesRegex(ValueError, "requires --lipoylation-module"):
            build_reference_chain(acyl_acp_chemistry=True)


if __name__ == "__main__":
    unittest.main()
