"""Targeted complex-I GPR checks: dependency propagation and preservation, no solve."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra.io import read_sbml_model

from scripts.gem_annotate.patches import R1889_ASSIGNMENT_PATH, apply_r1889_gpr_assignment
from scripts.gem_annotate.energy_candidates import export_candidate
from tests.test_coq9_curation import annotations, semantics

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts/atp_candidate_repair_20260924/candidates/E5.xml"


class R1889AssignmentTests(unittest.TestCase):
    def test_partial_rule_conflicts_ko_and_roundtrip(self):
        model = read_sbml_model(SOURCE)
        spec = json.loads(R1889_ASSIGNMENT_PATH.read_text())
        before = semantics(model), annotations(model)
        genes_before = {g.id: (g.name, copy.deepcopy(g.annotation), copy.deepcopy(g.notes))
                        for g in model.genes}
        self.assertEqual(model.reactions.R1889.gene_reaction_rule, "")
        for gid in spec["genes"]:
            with model:
                model.genes.get_by_id(gid).knock_out()
                self.assertTrue(model.reactions.R1889.functional)
                self.assertEqual(model.reactions.R1889.bounds, (0, 1000))

        # All rejected inputs must leave every scientific field untouched.
        for field in ("gpr", "stoichiometry", "bounds", "species", "name", "gene", "notes", "scope"):
            bad, changed_spec = model.copy(), copy.deepcopy(spec)
            rx = bad.reactions.R1889
            if field == "gpr":
                rx.gene_reaction_rule = "YALI1F32476g"
            elif field == "stoichiometry":
                rx.add_metabolites({bad.metabolites.get_by_id("m28[C_mi]"): 1})
            elif field == "bounds":
                rx.upper_bound = 20
            elif field == "species":
                bad.metabolites.get_by_id("m30[C_mi]").charge = -2
            elif field == "name":
                rx.name = "alternative non-pumping NADH dehydrogenase"
            elif field == "gene":
                bad.genes.YALI1B26679g.annotation = {"ncbigene": "wrong"}
            elif field == "notes":
                rx.notes = {**rx.notes, "gpr_scope_limit": "complete enzyme"}
            else:
                changed_spec["after_gpr"] = " or ".join(sorted(spec["genes"]))
            original = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaises(ValueError):
                apply_r1889_gpr_assignment(bad, changed_spec)
            self.assertEqual((semantics(bad), annotations(bad)), original)

        with patch("optlang.gurobi_interface.Model.optimize", side_effect=AssertionError("No solve")):
            self.assertEqual(apply_r1889_gpr_assignment(model)["status"], "applied")
            self.assertEqual(apply_r1889_gpr_assignment(model)["status"], "already_correct")
            rx = model.reactions.R1889
            self.assertEqual(rx.gene_reaction_rule, spec["after_gpr"])
            after = semantics(model)
            expected = copy.deepcopy(before[0])
            expected["reactions"]["R1889"] = after["reactions"]["R1889"]
            self.assertEqual(after, expected)
            self.assertEqual({k: v for k, v in annotations(model).items() if k != "R1889"},
                             {k: v for k, v in before[1].items() if k != "R1889"})
            self.assertEqual((rx.name, rx.annotation), before[1]["R1889"][:2])
            self.assertEqual(rx.notes, {**before[1]["R1889"][2], **spec["notes"]})
            self.assertEqual({g.id: (g.name, g.annotation, g.notes) for g in model.genes}, genes_before)
            for gid in spec["genes"]:
                with model:
                    model.genes.get_by_id(gid).knock_out()
                    self.assertFalse(rx.functional)
                    self.assertEqual(rx.bounds, (0, 0))
                    self.assertTrue(model.reactions.R570.functional)
            for gid in ("YALI1F32476g", "YALI1F11481g", "YALI1D02444g", "YALI1D00766g"):
                with model:
                    model.genes.get_by_id(gid).knock_out()
                    self.assertTrue(rx.functional)
                    self.assertEqual(rx.bounds, (0, 1000))
            with tempfile.TemporaryDirectory() as folder:
                loaded = export_candidate(model, Path(folder) / "R1889.xml")
            self.assertEqual(semantics(loaded), semantics(model))
            self.assertEqual(annotations(loaded), annotations(model))
            self.assertEqual(apply_r1889_gpr_assignment(loaded)["status"], "already_correct")


if __name__ == "__main__":
    unittest.main()
