"""In-memory checks for the 2026-10-20 opt-in lipoylation module; no optimization."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.energy_candidates import model_definition
from scripts.gem_annotate.lipoylation_module import SPEC_PATH, apply_lipoylation_module
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model

SPEC = json.loads(SPEC_PATH.read_text())
BASELINE = "E5_coq9_resp_pkg_r132_combined_20261020.xml"


def snapshot(model):
    """Definition plus formulas and the gene set, so a fail-closed check sees any edit."""
    return (model_definition(model), sorted(g.id for g in model.genes),
            {r.id: dict(r.notes) for r in model.reactions})


class LipoylationModuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.baseline = read_sbml_model(str(MODEL.candidate_file(BASELINE)))

    def fresh(self):
        return self.baseline.copy()

    def test_applies_only_the_curated_targets_and_balances(self):
        model = self.fresh()
        before = model_definition(model)
        self.assertEqual(apply_lipoylation_module(model)["status"], "applied")
        after = model_definition(model)
        changed_reactions = sorted(r for r in after["reactions"]
                                   if r not in before["reactions"] or after["reactions"][r] != before["reactions"][r])
        self.assertEqual(changed_reactions, sorted(["R362", "R364", "biomass_C", *SPEC["new_reactions"]]))
        changed_species = sorted(m for m in after["metabolites"]
                                 if m not in before["metabolites"] or after["metabolites"][m] != before["metabolites"][m])
        self.assertEqual(changed_species, sorted(["m179[C_mi]", "m586[C_mi]", *SPEC["new_metabolites"]]))
        self.assertEqual(after["objective"], before["objective"])
        for rid in SPEC["balanced_after"]:
            self.assertEqual(model.reactions.get_by_id(rid).check_mass_balance(), {}, rid)
        for gid in ("YALI1F34729g", "YALI1E13472g", "YALI1C04991g"):
            self.assertIn(gid, model.genes)
        with model:
            model.genes.YALI1E00574g.knock_out()
            self.assertEqual(model.reactions.R2288.bounds, (0, 0))
            self.assertEqual(model.reactions.R2286.bounds, (0.0, 1000.0))  # GCV3 is AND-ed in LA3 only.

    def test_rerun_is_a_no_op_also_after_sbml_round_trip(self):
        model = self.fresh()
        apply_lipoylation_module(model)
        state = snapshot(model)
        self.assertEqual(apply_lipoylation_module(model)["status"], "already_correct")
        self.assertEqual(snapshot(model), state)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lipoylation.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        loaded_state = snapshot(loaded)
        self.assertEqual(apply_lipoylation_module(loaded)["status"], "already_correct")
        self.assertEqual(snapshot(loaded), loaded_state)
        # Notes must survive SBML (an unescaped '&' once dropped a whole notes block silently).
        for rid, row in SPEC["new_reactions"].items():
            self.assertEqual(loaded.reactions.get_by_id(rid).notes.get("lipoylation_module"), SPEC["change_id"], rid)
        self.assertEqual(loaded.metabolites.get_by_id("m586[C_mi]").name, "lipoylprotein_C8H13OS2R")
        stripped = loaded.copy()
        stripped.reactions.R2288.notes = {}
        with self.assertRaises(ValueError):
            apply_lipoylation_module(stripped)

    def test_stale_or_partial_inputs_fail_closed_without_edits(self):
        def r362_bound(m): m.reactions.R362.upper_bound = 999.0
        def m586_formula(m): m.metabolites.get_by_id("m586[C_mi]").formula = "C8H13OS2R"
        def acp1_formula(m): m.metabolites.get_by_id("m179[C_mi]").formula = "HSR"
        def new_acp1_user(m):
            extra = Reaction("EXTRA"); m.add_reactions([extra])
            extra.add_metabolites({m.metabolites.get_by_id("m179[C_mi]"): -1.0})
        def biomass_term(m):
            m.reactions.biomass_C.add_metabolites({m.metabolites.get_by_id("m586[C_mi]"): -1e-4})
        def partial(m):
            m.add_reactions([Reaction("R2286")])
        def species_drift(m): m.metabolites.get_by_id("m60[C_mi]").charge = 1
        def acp1_user_residual(m):  # an ACP1 reaction outside the targets drifts before the edit
            m.reactions.R1394.add_metabolites({m.metabolites.get_by_id("m28[C_mi]"): 1.0})
        for change in (r362_bound, m586_formula, acp1_formula, new_acp1_user, biomass_term, partial, species_drift,
                       acp1_user_residual):
            with self.subTest(change.__name__):
                model = self.fresh()
                change(model)
                state = snapshot(model)
                with self.assertRaises(ValueError):
                    apply_lipoylation_module(model)
                self.assertEqual(snapshot(model), state)

    def test_unexpected_scope_is_rejected(self):
        spec = copy.deepcopy(SPEC)
        spec["change_id"] = "OTHER"
        model = self.fresh()
        state = snapshot(model)
        with self.assertRaises(ValueError):
            apply_lipoylation_module(model, spec)
        self.assertEqual(snapshot(model), state)

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).lipoylation_module)
        self.assertTrue(parse_args(["--lipoylation-module"]).lipoylation_module)


if __name__ == "__main__":
    unittest.main()
