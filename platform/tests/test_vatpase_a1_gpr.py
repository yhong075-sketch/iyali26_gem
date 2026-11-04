"""V-ATPase A1 GPR: preservation, rejection, SBML roundtrip, knockout logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.vatpase_a1_gpr import (A_SUBUNITS, AFTER, BEFORE, CORE, CURATION_PATH, NEW_GENE,
                                                 REACTIONS, apply_vatpase_a1_gpr)
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_"
            "r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_20261101.xml")
A1, A2 = A_SUBUNITS


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class VatpaseA1GprTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        precondition = "precondition differs"
        expected = {"name": precondition, "bounds": precondition, "gpr": precondition,
                    "stoichiometry": precondition, "species": precondition, "annotation": precondition,
                    "notes": "Conflicting", "extra_note": "notes outside this curation",
                    "notes_without_gpr": "carries this curation's notes without its GPR",
                    "after_without_notes": "already carries this GPR without this curation",
                    "one_reaction_applied": "Only one of R794 and R795",
                    "gene_identity": "identity or reaction context", "gene_context": "identity or reaction context",
                    "new_gene_present": f"{NEW_GENE} already exists", "missing_reaction": "R795 missing"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R795
            if field == "name":
                reaction.name = "V-ATPase"
            elif field == "bounds":
                reaction.bounds = (0.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = BEFORE.replace(" or (YALI1F38820g and YALI1F13017g)", "")
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): -1})
            elif field == "species":
                bad.metabolites.get_by_id("m1007[C_va]").charge = 0
            elif field == "annotation":
                reaction.annotation = {**reaction.annotation, "ec-code": "7.1.2.2"}
            elif field == "notes":
                reaction.notes = {**reaction.notes, "vatpase_a1_gpr_status": "another review"}
            elif field == "extra_note":
                reaction.notes = {**reaction.notes, "unreviewed_note": "added elsewhere"}
            elif field == "notes_without_gpr":
                reaction.notes = {**reaction.notes, **copy.deepcopy(self.spec["notes"])}
            elif field == "after_without_notes":
                reaction.gene_reaction_rule = AFTER
            elif field == "one_reaction_applied":
                reaction.gene_reaction_rule = AFTER
                reaction.notes = {**reaction.notes, **copy.deepcopy(self.spec["notes"])}
            elif field == "gene_identity":
                # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
                gene = bad.genes.get_by_id(A2)
                gene.annotation = {**gene.annotation, "uniprot": "Q00000"}
            elif field == "gene_context":
                bad.reactions.R171.gene_reaction_rule = f"{bad.reactions.R171.gene_reaction_rule} and {A1}"
            elif field == "new_gene_present":
                bad.reactions.R172.gene_reaction_rule = f"{bad.reactions.R172.gene_reaction_rule} or {NEW_GENE}"
            else:
                bad.remove_reactions([reaction])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_vatpase_a1_gpr(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_applied_state_pins_the_new_gene_and_the_rule_genes(self):
        cases = {"new_gene_annotation": f"{NEW_GENE} already exists", "new_gene_name": f"{NEW_GENE} already exists",
                 "gene_context": "identity or reaction context"}
        for case, message in cases.items():
            model = self.model.copy()
            apply_vatpase_a1_gpr(model)
            gene = model.genes.get_by_id(NEW_GENE)
            if case == "new_gene_annotation":
                gene.annotation = {**gene.annotation, "uniprot": "Q00000"}
            elif case == "new_gene_name":
                gene.name = "VMA9"
            else:
                model.reactions.R172.gene_reaction_rule = f"{model.reactions.R172.gene_reaction_rule} or {A2}"
            before = state(model)
            with self.subTest(case=case), self.assertRaisesRegex(ValueError, message):
                apply_vatpase_a1_gpr(model)
            self.assertEqual(state(model), before)

    def test_unexpected_scope_is_rejected(self):
        cases = (("status", "user_authorized"), ("assignment_id", "VATPASE-A1-GPR"))
        for key, value in cases:
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_vatpase_a1_gpr(self.model.copy(), spec)
        for edit in ("after_gpr", "reaction_notes", "gene_extra", "gene_after", "new_gene"):
            spec = copy.deepcopy(self.spec)
            if edit == "after_gpr":
                spec["reactions"]["R795"]["after_gpr"] = AFTER.replace(f" or {A2}", "")
            elif edit == "reaction_notes":
                spec["reactions"]["R794"]["notes"] = {**spec["notes"], "extra": "x"}
            elif edit == "gene_extra":
                spec["genes"]["YALI1E41893g"] = copy.deepcopy(spec["genes"][A1])
            elif edit == "gene_after":
                spec["genes"][A1]["reactions_after"] = []
            else:
                spec["new_gene"]["id"] = "YALI1E19327g"
            with self.subTest(edit=edit), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_vatpase_a1_gpr(self.model.copy(), spec)

    def test_changes_only_targets_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        for rid in REACTIONS:
            target = expected.reactions.get_by_id(rid)
            target.gene_reaction_rule = AFTER
            target.notes = {**target.notes, **self.spec["notes"]}
        expected.genes.get_by_id(NEW_GENE).annotation = copy.deepcopy(self.spec["new_gene"]["annotation"])
        result = apply_vatpase_a1_gpr(model)
        self.assertEqual((result["status"], result["added_genes"]), ("applied", [NEW_GENE]))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual([g.id for g in model.genes], [g.id for g in self.model.genes] + [NEW_GENE])
        for rid in REACTIONS:
            self.assertEqual({g.id for g in model.reactions.get_by_id(rid).genes},
                             set(CORE) | set(A_SUBUNITS) | {NEW_GENE})
        self.assertEqual(sorted(r.id for r in model.genes.YALI1A09766g.reactions), ["R794", "R795", "R_NTP1"])
        self.assertEqual(apply_vatpase_a1_gpr(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "vatpase.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        for rid in REACTIONS:
            self.assertEqual(loaded.reactions.get_by_id(rid).gene_reaction_rule, AFTER)
        self.assertEqual(apply_vatpase_a1_gpr(loaded)["status"], "already_correct")

    def test_each_core_gene_closes_both_and_single_a_subunits_close_nothing(self):
        for curated in (False, True):
            model = self.model.copy()
            if curated:
                apply_vatpase_a1_gpr(model)
            for rid in REACTIONS:
                model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)  # R795 is closed; test the rule alone
            genes = CORE + A_SUBUNITS + ((NEW_GENE,) if curated else ())
            for gene_id in genes:
                with model, self.subTest(curated=curated, gene=gene_id):
                    model.genes.get_by_id(gene_id).knock_out()
                    for rid in REACTIONS:
                        # Before the edit the second branch keeps both open for every single knockout.
                        self.assertEqual(model.reactions.get_by_id(rid).bounds == (0, 0),
                                         curated and gene_id not in A_SUBUNITS)
            with model:
                for gene_id in A_SUBUNITS:
                    model.genes.get_by_id(gene_id).knock_out()
                self.assertTrue(all(model.reactions.get_by_id(rid).bounds == (0, 0) for rid in REACTIONS))
            self.assertNotEqual(model.reactions.R_NTP1.bounds, (0, 0))

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_vatpase_a1_gpr(model)["status"], "applied")
            self.assertEqual(model.reactions.R794.gene_reaction_rule, AFTER)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertNotIn(NEW_GENE, model.genes)
        self.assertEqual(apply_vatpase_a1_gpr(model)["status"], "applied")

    def test_flag_is_opt_in_and_excludes_the_hypothesis(self):
        self.assertFalse(parse_args([]).vatpase_a1_gpr)
        self.assertTrue(parse_args(["--vatpase-a1-gpr"]).vatpase_a1_gpr)
        with self.assertRaisesRegex(ValueError, "not with --vatpase-gpr-hypothesis"):
            build_reference_chain(vatpase_a1_gpr=True, vatpase_gpr_hypothesis=True, no_solve=True, allow_network=False)


if __name__ == "__main__":
    unittest.main()
