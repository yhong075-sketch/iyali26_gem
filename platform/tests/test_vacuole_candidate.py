"""Bounded vacuole candidate regression with optimization and network disabled."""

import copy
import hashlib
from pathlib import Path
import tempfile
import unittest

from cobra.io import read_sbml_model

from scripts.gem_annotate.energy_candidates import model_definition, protected_definitions
from scripts.gem_annotate.execution import execution_limits
from scripts.gem_annotate.vacuole_candidates import (
    TARGET_BOUNDS, apply_vacuole_candidate, build_candidate_file, load_spec,
)

from scripts.gem_annotate.config import resolve_recorded_path
ROOT = Path(__file__).resolve().parents[2]
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.model_layout import SCRATCH_DIR
TASK = MODEL.reports / "vacuole_open_supply_20260924"


class VacuoleCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.guard = execution_limits(no_solve=True, allow_network=False)
        cls.attempts = cls.guard.__enter__()
        cls.spec = load_spec()
        cls.source = resolve_recorded_path(cls.spec["source_path"])

    @classmethod
    def tearDownClass(cls):
        cls.guard.__exit__(None, None, None)
        assert cls.attempts == {"optimization": 0, "network": 0}, cls.attempts

    def test_default_noop_bounded_scope_and_idempotency(self):
        model = read_sbml_model(self.source)
        before = model_definition(model)
        locks = protected_definitions(model)
        notes = {r.id: copy.deepcopy(r.notes) for r in model.reactions}
        self.assertEqual(apply_vacuole_candidate(model)["status"], "disabled")
        self.assertEqual(model_definition(model), before)
        result = apply_vacuole_candidate(model, enabled=True)
        expected = copy.deepcopy(before)
        for rid, bounds in TARGET_BOUNDS.items():
            expected["reactions"][rid]["bounds"] = bounds
        self.assertEqual(model_definition(model), expected)
        self.assertEqual(protected_definitions(model), locks)
        self.assertEqual({r.id: r.notes for r in model.reactions}, notes)
        self.assertEqual(len(result["items"]), 4)
        repeat = apply_vacuole_candidate(model, enabled=True)
        self.assertTrue(all(row["status"] == "already_correct" for row in repeat["items"]))
        self.assertEqual(model_definition(model), expected)

    def test_preflight_mismatches_are_atomic(self):
        for kind in ("stoichiometry", "gpr", "bounds", "species", "energy_lock", "authorization"):
            with self.subTest(kind=kind):
                model = read_sbml_model(self.source)
                spec = copy.deepcopy(self.spec)
                if kind == "stoichiometry":
                    model.reactions.R795.add_metabolites({model.metabolites.get_by_id("m10[C_cy]"): 1})
                elif kind == "gpr":
                    model.reactions.R2039.gene_reaction_rule = "candidate_gene"
                elif kind == "bounds":
                    model.reactions.R2030.lower_bound = 0
                elif kind == "species":
                    model.metabolites.get_by_id("m1384[C_va]").compartment = "C_cy"
                elif kind == "energy_lock":
                    model.reactions.R_CAT2p.notes.clear()
                else:
                    spec["bounds"]["R795"] = [0, 1000]
                before = model_definition(model)
                with self.assertRaises(ValueError):
                    apply_vacuole_candidate(model, enabled=True, spec=spec)
                self.assertEqual(model_definition(model), before)

    def test_source_identity_export_reload_and_no_overwrite(self):
        before_sha = hashlib.sha256(self.source.read_bytes()).hexdigest()
        with tempfile.TemporaryDirectory(dir=SCRATCH_DIR) as directory:
            output = Path(directory) / "candidate.xml"
            record = build_candidate_file(self.source, output, enabled=True)
            self.assertEqual(record["source_sha256"], before_sha)
            loaded = read_sbml_model(output)
            self.assertEqual(len(loaded.reactions), 2314)
            self.assertEqual(len(loaded.metabolites), 1877)
            self.assertEqual(protected_definitions(loaded), self.spec["energy_protected_definitions"])
            self.assertFalse(any(r.id.startswith(("DIAG_", "POOL_", "TEST_")) for r in loaded.reactions))
            self.assertEqual(hashlib.sha256(self.source.read_bytes()).hexdigest(), before_sha)
            with self.assertRaises(FileExistsError):
                build_candidate_file(self.source, output, enabled=True)
            spec = copy.deepcopy(self.spec); spec["source_sha256"] = "0" * 64
            invalid = Path(directory) / "invalid.xml"
            with self.assertRaisesRegex(ValueError, "input SHA"):
                build_candidate_file(self.source, invalid, enabled=True, spec=spec)
            self.assertFalse(invalid.exists())
            with self.assertRaisesRegex(ValueError, "workspace"):
                build_candidate_file(self.source, ROOT / "../forbidden.xml", enabled=True)


if __name__ == "__main__":
    unittest.main()
