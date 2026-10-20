"""R1964 provisional GPR: preservation, rejection, SBML roundtrip, OR logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r1964_gpr import CURATION_PATH, apply_r1964_gpr_assignment
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = "E5_coq9_resp_pkg_r132_20261012.xml"
NEW_GENE = "YALI1E06616g"


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


class R1964GprTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        precondition, outside = "precondition differs", "notes outside this curation"
        old_gene, new_gene = "YALI1F03452g identity or reaction context", f"{NEW_GENE} already exists"
        expected = {"name": precondition, "bounds": precondition, "gpr": precondition,
                    "stoichiometry": precondition, "species": precondition, "compartment": precondition,
                    "annotation": precondition, "notes": "Conflicting", "protein_class": outside,
                    "extra_note": outside, "after_without_notes": "already carries this GPR",
                    "existing_gene_identity": old_gene, "existing_gene_uniprot": old_gene,
                    "existing_gene_context": old_gene, "new_gene_identity": new_gene,
                    "new_gene_name": new_gene, "new_gene_context": new_gene}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R1964
            # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
            old = bad.genes.YALI1F03452g
            if field == "name":
                reaction.name = "different dehydrogenase"
            elif field == "bounds":
                reaction.bounds = (0.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = "YALI1E15452g"
            elif field == "stoichiometry":
                # The source equation's H+ is not this step's input.
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): -1})
            elif field == "species":
                xylulose = bad.metabolites.get_by_id("m1015[C_cy]")
                xylulose.formula, xylulose.charge = "C5H9O5", -1
            elif field == "compartment":
                bad.metabolites.get_by_id("m1823[C_cy]").compartment = "C_mi"
            elif field == "annotation":
                reaction.annotation["metanetx.reaction"] = "MNXR198270"
            elif field == "notes":
                reaction.notes["curated_gpr_assignment"] = "another review"
            elif field == "protein_class":
                reaction.notes["PROTEIN_CLASS"] = "1.1.1.250"
            elif field == "extra_note":
                reaction.notes["unreviewed_note"] = "added elsewhere"
            elif field == "after_without_notes":
                reaction.gene_reaction_rule = self.spec["after_gpr"]
                reaction.annotation = copy.deepcopy(self.spec["after_annotation"])
                bad.genes.get_by_id(NEW_GENE).annotation = copy.deepcopy(
                    self.spec["new_genes"][NEW_GENE]["annotation"])
            elif field == "existing_gene_identity":
                old.annotation = {**old.annotation, "ncbigene": "2912324"}
            elif field == "existing_gene_uniprot":
                old.annotation = {**old.annotation, "uniprot": "Q6C6W8"}
            elif field == "existing_gene_context":
                bad.reactions.R801.gene_reaction_rule = "YALI1E15452g or YALI1F03452g"
            elif field == "new_gene_identity":
                # Create the gene without reactions (cobra keeps it after the GPR is reset).
                bad.reactions.R801.gene_reaction_rule = "YALI1E15452g or " + NEW_GENE
                bad.reactions.R801.gene_reaction_rule = "YALI1E15452g"
                bad.genes.get_by_id(NEW_GENE).annotation = {"ncbigene": "1"}
            elif field == "new_gene_name":
                # The full after-state, but with an unverified gene symbol.
                apply_r1964_gpr_assignment(bad)
                bad.genes.get_by_id(NEW_GENE).name = "SOU1"
            else:
                bad.reactions.R801.gene_reaction_rule = "YALI1E15452g or " + NEW_GENE
                bad.genes.get_by_id(NEW_GENE).annotation = copy.deepcopy(
                    self.spec["new_genes"][NEW_GENE]["annotation"])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r1964_gpr_assignment(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for key, value in (("after_gpr", "YALI1F03452g"), ("status", "user_authorized"),
                           ("new_genes", {})):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                apply_r1964_gpr_assignment(self.model.copy(), spec)

    def test_changes_only_r1964_and_new_gene_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        target = expected.reactions.R1964
        target.gene_reaction_rule = self.spec["after_gpr"]
        expected.genes.get_by_id(NEW_GENE).annotation = copy.deepcopy(self.spec["new_genes"][NEW_GENE]["annotation"])
        target.annotation = copy.deepcopy(self.spec["after_annotation"])
        target.notes.update(self.spec["notes"])
        self.assertNotIn(NEW_GENE, model.genes)
        result = apply_r1964_gpr_assignment(model)
        self.assertEqual((result["status"], result["added_genes"]), ("applied", [NEW_GENE]))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual(model.reactions.R1964.check_mass_balance(), {})
        self.assertEqual(sorted(r.id for r in model.genes.YALI1F03452g.reactions), ["R1964", "R2190"])
        self.assertEqual(apply_r1964_gpr_assignment(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r1964.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r1964_gpr_assignment(loaded)["status"], "already_correct")

    def test_either_gene_keeps_r1964_open(self):
        model = self.model.copy()
        apply_r1964_gpr_assignment(model)
        for knocked, closed in ((["YALI1F03452g"], {"R2190"}), ([NEW_GENE], set()),
                                (["YALI1F03452g", NEW_GENE], {"R2190", "R1964"})):
            with model, self.subTest(knocked=knocked):
                for gene in knocked:
                    model.genes.get_by_id(gene).knock_out()
                self.assertEqual({r for r in ("R1964", "R2190") if model.reactions.get_by_id(r).bounds == (0, 0)},
                                 closed)

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r1964_gpr)
        self.assertTrue(parse_args(["--r1964-gpr"]).r1964_gpr)


if __name__ == "__main__":
    unittest.main()
