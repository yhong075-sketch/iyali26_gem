"""R2193 iron GPR: preservation, rejection, SBML roundtrip, pair logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r2193_iron_gpr import CURATION_PATH, apply_r2193_iron_gpr
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = "E5_coq9_resp_pkg_r132_combined_gpi_coproducts_heme_iron_20261022.xml"
OLD_GENE = "YALI1F07747g"
D09, D08 = ("YALI1D09277g", "YALI1D09308g"), ("YALI1D08684g", "YALI1D08564g")


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


class R2193IronGprTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        precondition, outside = "precondition differs", "notes outside this curation"
        old_gene, new_gene = f"{OLD_GENE} identity or reaction context", f"{D08[1]} already exists"
        expected = {"name": precondition, "bounds": precondition, "gpr": precondition,
                    "stoichiometry": precondition, "species": precondition, "annotation": precondition,
                    "notes": "Conflicting", "extra_note": outside,
                    "notes_without_gpr": "carries this curation's notes without its GPR",
                    "after_without_notes": "already carries this GPR without this curation",
                    "old_gene_identity": old_gene, "old_gene_context": old_gene,
                    "new_gene_present": new_gene, "new_gene_name": "already exists with other identity",
                    "missing_new_gene": "missing although R2193 carries the curated GPR"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R2193
            # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
            old = bad.genes.get_by_id(OLD_GENE)
            if field == "name":
                reaction.name = "iron (III) transport"
            elif field == "bounds":
                reaction.bounds = (-1000.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = "YALI1D09277g"
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "species":
                bad.metabolites.get_by_id("m948[C_cy]").formula = "Fe"
            elif field == "annotation":
                reaction.annotation = {**reaction.annotation, "metanetx.reaction": "MNXR99502"}
            elif field == "notes":
                reaction.notes = {"r2193_iron_gpr_status": "another review"}
            elif field == "extra_note":
                reaction.notes = {"unreviewed_note": "added elsewhere"}
            elif field == "notes_without_gpr":
                reaction.notes = copy.deepcopy(self.spec["notes"])
            elif field == "after_without_notes":
                reaction.gene_reaction_rule = self.spec["after_gpr"]
                for gene_id, row in self.spec["new_genes"].items():
                    bad.genes.get_by_id(gene_id).annotation = copy.deepcopy(row["annotation"])
            elif field == "old_gene_identity":
                old.annotation = {**old.annotation, "ncbigene": "2910500"}
            elif field == "old_gene_context":
                bad.reactions.R863.gene_reaction_rule = OLD_GENE
            elif field == "new_gene_present":
                # A new gene already on another reaction.
                bad.reactions.R863.gene_reaction_rule = D08[1]
            elif field == "new_gene_name":
                # The full after-state, but with an unverified gene symbol.
                apply_r2193_iron_gpr(bad)
                bad.genes.get_by_id(D08[1]).name = "FTR1"
            else:
                # The curated GPR text with one of its genes dropped from the model's gene list.
                apply_r2193_iron_gpr(bad)
                bad.genes.remove(bad.genes.get_by_id(D09[0]))
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r2193_iron_gpr(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for key, value in (("after_gpr", "YALI1D09277g and YALI1D09308g"), ("status", "user_authorized"),
                           ("new_genes", {}), ("before_gpr", ""), ("old_gene", {"id": "YALI1D17462g"})):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_r2193_iron_gpr(self.model.copy(), spec)

    def test_changes_only_r2193_and_genes_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        target = expected.reactions.R2193
        target.gene_reaction_rule = self.spec["after_gpr"]
        for gene_id, row in self.spec["new_genes"].items():
            expected.genes.get_by_id(gene_id).annotation = copy.deepcopy(row["annotation"])
        target.notes = {**target.notes, **self.spec["notes"]}
        result = apply_r2193_iron_gpr(model)
        self.assertEqual((result["status"], result["added_genes"]), ("applied", sorted(D09 + D08)))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertIn(OLD_GENE, model.genes)
        self.assertEqual(list(model.genes.get_by_id(OLD_GENE).reactions), [])
        self.assertEqual(len(model.genes), len(self.model.genes) + 4)
        # Sorted, not hash-seed order, so builds with different PYTHONHASHSEED write the same bytes.
        self.assertEqual([g.id for g in model.genes][-4:], sorted(D09 + D08))
        self.assertEqual(apply_r2193_iron_gpr(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r2193.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r2193_iron_gpr(loaded)["status"], "already_correct")

    def test_either_pair_keeps_r2193_open(self):
        model = self.model.copy()
        apply_r2193_iron_gpr(model)
        for knocked, closed in (([OLD_GENE], False), ([D09[0]], False), ([D08[1]], False),
                                (list(D09), False), (list(D08), False),
                                ([D09[0], D08[0]], True), ([D09[1], D08[1]], True)):
            with model, self.subTest(knocked=knocked):
                for gene in knocked:
                    model.genes.get_by_id(gene).knock_out()
                self.assertEqual(model.reactions.R2193.bounds == (0, 0), closed)

    def test_context_exit_restores_the_baseline(self):
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        before = state(model)
        with model:
            self.assertEqual(apply_r2193_iron_gpr(model)["status"], "applied")
            self.assertEqual(model.reactions.R2193.gene_reaction_rule, self.spec["after_gpr"])
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_r2193_iron_gpr(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r2193_iron_gpr)
        self.assertTrue(parse_args(["--r2193-iron-gpr"]).r2193_iron_gpr)


if __name__ == "__main__":
    unittest.main()
