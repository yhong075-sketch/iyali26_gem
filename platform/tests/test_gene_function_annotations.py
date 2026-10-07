"""Candidate status survives export without changing the optimization problem."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Model, Metabolite, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.genes import (
    GENE_FUNCTION_CURATION_PATH, apply_curated_gene_function_annotations,
)
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics


class GeneFunctionAnnotationTests(unittest.TestCase):
    def setUp(self):
        self.model = Model("candidate_function")
        a = Metabolite("a", formula="C", charge=0, compartment="c")
        b = Metabolite("b", formula="C", charge=0, compartment="c")
        for rid in ("R794", "R795"):
            reaction = Reaction(rid)
            reaction.add_metabolites({a: -1, b: 1})
            reaction.gene_reaction_rule = "other_gene or YALI1F38820g"
            self.model.add_reactions([reaction])
        self.model.reactions.R795.bounds = (0.0, 0.0)
        self.model.objective = self.model.reactions.R794
        self.model.genes.other_gene.name = "unchanged control"
        self.gene = self.model.genes.YALI1F38820g
        self.gene.name = "COBRAProtein676"
        self.gene.annotation = {"uniprot": ["A0A1H6PMT1"], "ncbigene": ["2907662"]}
        self.gene.notes = {"existing_note": "preserve this"}

    def test_candidate_roundtrip_idempotency_and_mathematics(self):
        before = semantics(self.model)
        identity = copy.deepcopy(self.gene.annotation)
        self.assertEqual(apply_curated_gene_function_annotations(self.model), 1)
        self.assertEqual(semantics(self.model), before)
        self.assertEqual(self.gene.annotation, identity)
        self.assertEqual(self.gene.notes["existing_note"], "preserve this")
        self.assertEqual(self.gene.notes["function_candidate_previous_name"], "COBRAProtein676")
        self.assertIn("candidate; experimental confirmation required", self.gene.name)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.xml"
            write_deterministic_sbml_model(self.model, path)
            loaded = read_sbml_model(str(path))
        self.assertEqual(semantics(loaded), before)
        gene = loaded.genes.YALI1F38820g
        self.assertEqual(gene.name, self.gene.name)
        self.assertEqual(gene.notes, self.gene.notes)
        self.assertEqual(gene.notes["function_experimental_confirmation"], "required")
        self.assertEqual(gene.notes["function_candidate_status"], "provisional")
        self.assertEqual(apply_curated_gene_function_annotations(loaded), 0)
        self.assertEqual(loaded.genes.other_gene.name, "unchanged control")

    def test_wrong_identity_or_missing_caveat_rejected_before_edits(self):
        original = (self.gene.name, copy.deepcopy(self.gene.notes))
        self.gene.annotation["uniprot"] = ["wrong_accession"]
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            apply_curated_gene_function_annotations(self.model)
        self.assertEqual((self.gene.name, self.gene.notes), original)
        self.gene.annotation["uniprot"] = "A0A1H6PMT1"
        spec = json.loads(GENE_FUNCTION_CURATION_PATH.read_text())
        del spec["genes"][self.gene.id]["notes"]["function_experimental_confirmation"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.json"
            path.write_text(json.dumps(spec))
            with self.assertRaisesRegex(ValueError, "evidence limits"):
                apply_curated_gene_function_annotations(self.model, path)
        self.assertEqual((self.gene.name, self.gene.notes), original)


if __name__ == "__main__":
    unittest.main()
