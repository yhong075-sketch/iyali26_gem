"""R24 homocitrate aconitase: preservation, rejection, SBML roundtrip, knockout logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r24_homocitrate_aconitase import (AFTER, CURATION_PATH, TARGET,
                                                            apply_r24_homocitrate_aconitase)
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_"
            "gpr_r673_prolyl_split_aco_split_20261029.xml")
LYS4, ACO1, ACO2 = "YALI1E03387g", "YALI1D11984g", "YALI1E17935g"


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class R24HomocitrateAconitaseTests(unittest.TestCase):
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
                    "gene_identity": "identity or reaction context",
                    "aco2_still_linked": "identity or reaction context",
                    "missing_gene": f"{ACO2} missing",
                    "missing_reaction": f"{TARGET} missing"}
        baseline = state(self.model)
        row = self.spec["reactions"][TARGET]
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.get_by_id(TARGET)
            if field == "name":
                reaction.name = "homocitrate dehydratase"
            elif field == "bounds":
                reaction.bounds = (0.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = ACO1
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m28[C_mi]"): -3})
            elif field == "species":
                bad.metabolites.get_by_id("m75[C_mi]").charge = -3
            elif field == "annotation":
                reaction.annotation = {**reaction.annotation, "ec-code": "4.2.1.-"}
            elif field == "notes":
                reaction.notes = {**reaction.notes, "r24_homocitrate_aco_status": "another review"}
            elif field == "extra_note":
                reaction.notes = {**reaction.notes, "unreviewed_note": "added elsewhere"}
            elif field == "notes_without_gpr":
                reaction.notes = {**reaction.notes, **copy.deepcopy(row["notes"])}
            elif field == "after_without_notes":
                reaction.gene_reaction_rule = row["after_gpr"]
            elif field == "gene_identity":
                # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
                gene = bad.genes.get_by_id(ACO1)
                gene.annotation = {**gene.annotation, "ncbigene": "2900000"}
            elif field == "aco2_still_linked":
                # The state before the aconitase split: the ACO2-type gene still on the aconitase reactions.
                bad.reactions.R216.gene_reaction_rule = f"{ACO1} or {ACO2}"
            elif field == "missing_gene":
                bad.genes.remove(bad.genes.get_by_id(ACO2))
            else:
                bad.remove_reactions([reaction])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r24_homocitrate_aconitase(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for key, value in (("status", "user_authorized"), ("assignment_id", "R24-HOMOCITRATE-ACO")):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_r24_homocitrate_aconitase(self.model.copy(), spec)
        for edit in ("after_gpr", "before_gpr", "lys4_after", "aco2_after", "aco1_after", "extra_reaction",
                     "missing_gene"):
            spec = copy.deepcopy(self.spec)
            if edit == "after_gpr":
                spec["reactions"][TARGET]["after_gpr"] = ACO2
            elif edit == "before_gpr":
                spec["reactions"][TARGET]["before_gpr"] = ACO1
            elif edit == "lys4_after":
                spec["genes"][LYS4]["reactions_after"] = []
            elif edit == "aco2_after":
                spec["genes"][ACO2]["reactions_after"] = []
            elif edit == "aco1_after":
                spec["genes"][ACO1]["reactions_after"] = spec["genes"][ACO1]["reactions_before"]
            elif edit == "extra_reaction":
                spec["reactions"]["R395"] = copy.deepcopy(spec["reactions"][TARGET])
            else:
                del spec["genes"][ACO2]
            with self.subTest(edit=edit), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_r24_homocitrate_aconitase(self.model.copy(), spec)

    def test_changes_only_the_r24_rule_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        target = expected.reactions.get_by_id(TARGET)
        target.gene_reaction_rule = self.spec["reactions"][TARGET]["after_gpr"]
        target.notes = {**target.notes, **self.spec["reactions"][TARGET]["notes"]}
        result = apply_r24_homocitrate_aconitase(model)
        self.assertEqual(result["status"], "applied")
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual([g.id for g in model.genes], [g.id for g in self.model.genes])
        self.assertEqual(model.reactions.R395.gene_reaction_rule, LYS4)
        self.assertEqual(apply_r24_homocitrate_aconitase(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r24.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r24_homocitrate_aconitase(loaded)["status"], "already_correct")
        # On an applied model the guards still hold: an extra R24 member fails the rule guard, and the ACO2-type gene
        # relinked to an aconitase reaction fails the gene guard on the applied path.
        for case, message in (("r24_rule", "precondition differs"), ("aco2_relinked", "identity or reaction context")):
            bad = loaded.copy()
            if case == "r24_rule":
                bad.reactions.get_by_id(TARGET).gene_reaction_rule = f"{AFTER} or {LYS4}"
            else:
                bad.reactions.R216.gene_reaction_rule = f"{ACO1} or {ACO2}"
            before = state(bad)
            with self.subTest(applied=case), self.assertRaisesRegex(ValueError, message):
                apply_r24_homocitrate_aconitase(bad)
            self.assertEqual(state(bad), before)

    def test_only_both_aconitases_together_close_r24(self):
        model = self.model.copy()
        apply_r24_homocitrate_aconitase(model)
        cases = {(LYS4,): False, (ACO1,): False, (ACO2,): False, (ACO1, ACO2): True}
        for genes, closes in cases.items():
            with model, self.subTest(genes=genes):
                for gene_id in genes:
                    model.genes.get_by_id(gene_id).knock_out()
                self.assertEqual(model.reactions.get_by_id(TARGET).bounds == (0, 0), closes)
                # R395 stays on LYS4 alone.
                self.assertEqual(model.reactions.R395.bounds == (0, 0), LYS4 in genes)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_r24_homocitrate_aconitase(model)["status"], "applied")
            self.assertEqual(model.reactions.get_by_id(TARGET).gene_reaction_rule, AFTER)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_r24_homocitrate_aconitase(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r24_homocitrate_aconitase)
        self.assertTrue(parse_args(["--r24-homocitrate-aconitase"]).r24_homocitrate_aconitase)


if __name__ == "__main__":
    unittest.main()
