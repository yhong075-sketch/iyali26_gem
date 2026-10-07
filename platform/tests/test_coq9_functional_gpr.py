"""Conditional COQ9 dependency: fail closed, preserve chemistry, transmit KO."""
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Gene
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.coq9 import (
    FUNCTIONAL_GPR_PATH, apply_coq9_functional_gpr, exact_residual,
)
from scripts.gem_annotate.execution import execution_limits
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import FIXTURE, annotations, semantics


class CoQ9FunctionalGPRTest(unittest.TestCase):
    def test_scope_conflicts_ko_and_roundtrip(self):
        self.assertFalse(parse_args([]).coq9_functional_gpr)
        self.assertTrue(parse_args(["--coq9-functional-gpr"]).coq9_functional_gpr)
        with self.assertRaises(ValueError):
            build_reference_chain(coq9_functional_gpr=True)
        with execution_limits(no_solve=True, allow_network=False) as attempts:
            model = read_sbml_model(str(FIXTURE))
            spec = json.loads(FUNCTIONAL_GPR_PATH.read_text())
            for gid, identity in spec["genes"].items():
                if gid not in model.genes:
                    model.genes.append(Gene(gid))
                model.genes.get_by_id(gid).annotation["refseq"] = identity["refseq"]
            before, metadata = semantics(model), annotations(model)
            for conflict in ("gpr", "charge", "bounds", "refseq", "notes"):
                bad = model.copy()
                r = bad.reactions.R695
                if conflict == "gpr":
                    r.gene_reaction_rule = "YALI1F34675g"
                elif conflict == "charge":
                    next(iter(r.metabolites)).charge = 20
                elif conflict == "bounds":
                    r.upper_bound = 2
                elif conflict == "refseq":
                    gene = bad.genes.YALI1F34675g
                    gene.annotation = {**gene.annotation, "refseq": "wrong_version"}
                else:
                    r.notes["coq9_functional_assumption"] = "experimentally confirmed"
                state = semantics(bad), annotations(bad)
                with self.subTest(conflict=conflict), self.assertRaises(ValueError):
                    apply_coq9_functional_gpr(bad)
                self.assertEqual((semantics(bad), annotations(bad)), state)
            self.assertEqual(apply_coq9_functional_gpr(model)["status"], "applied")
            self.assertEqual(apply_coq9_functional_gpr(model)["status"], "already_correct")
            after = semantics(model)
            after["reactions"]["R695"] = before["reactions"]["R695"]
            self.assertEqual(after, before)
            self.assertEqual(exact_residual(model.reactions.R695), {})
            self.assertEqual({k: v for k, v in annotations(model).items() if k != "R695"},
                             {k: v for k, v in metadata.items() if k != "R695"})
            for gid in spec["genes"]:
                with model:
                    model.genes.get_by_id(gid).knock_out()
                    self.assertFalse(model.reactions.R695.functional)
                    self.assertEqual(model.reactions.R695.bounds, (0, 0))
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "candidate.xml"
                write_deterministic_sbml_model(model, path)
                loaded = read_sbml_model(str(path))
            self.assertEqual(semantics(loaded), semantics(model))
            self.assertEqual(annotations(loaded), annotations(model))
            self.assertEqual(apply_coq9_functional_gpr(loaded)["status"], "already_correct")
            self.assertEqual(attempts, {"optimization": 0, "network": 0})


if __name__ == "__main__":
    unittest.main()
