"""Aconitase ACO2 split: preservation, rejection, SBML roundtrip, knockout logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.aconitase_aco2_split import (AFTER, CURATION_PATH, TARGETS, UNLINKED,
                                                       apply_aconitase_aco2_split)
from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_"
            "gpr_r673_prolyl_split_notes_v3_20261029.xml")
ACO1, ACO2 = "YALI1D11984g", "YALI1E17935g"
ACONITASES = (*TARGETS, "R239")


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class AconitaseAco2SplitTests(unittest.TestCase):
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
                    "partly_applied": "partly applied",
                    "gene_identity": "identity or reaction context", "gene_context": "identity or reaction context",
                    "missing_reaction": "R216 missing"}
        baseline = state(self.model)
        row = self.spec["reactions"]["R216"]
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R216
            if field == "name":
                reaction.name = "aconitase, mitochondrial"
            elif field == "bounds":
                reaction.bounds = (0.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = ACO2
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m26[C_mi]"): 1})
            elif field == "species":
                bad.metabolites.get_by_id("m409[C_mi]").formula = "C6H5O6"
            elif field == "annotation":
                reaction.annotation = {**reaction.annotation, "ec-code": "4.2.1.36"}
            elif field == "notes":
                reaction.notes = {**reaction.notes, "aconitase_aco2_split_status": "another review"}
            elif field == "extra_note":
                reaction.notes = {**reaction.notes, "unreviewed_note": "added elsewhere"}
            elif field == "notes_without_gpr":
                reaction.notes = {**reaction.notes, **copy.deepcopy(row["notes"])}
            elif field == "after_without_notes":
                reaction.gene_reaction_rule = row["after_gpr"]
            elif field == "partly_applied":
                reaction.gene_reaction_rule = row["after_gpr"]
                reaction.notes = {**reaction.notes, **copy.deepcopy(row["notes"])}
            elif field == "gene_identity":
                # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
                gene = bad.genes.get_by_id(ACO2)
                gene.annotation = {**gene.annotation, "ncbigene": "2900000"}
            elif field == "gene_context":
                bad.reactions.R24.gene_reaction_rule = f"YALI1E03387g or {ACO2}"
            else:
                bad.remove_reactions([reaction])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_aconitase_aco2_split(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        cases = (("status", "user_authorized"), ("assignment_id", "ACONITASE-ACO2-SPLIT"),
                 ("unlinked_gene", {"id": ACO1}))
        for key, value in cases:
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_aconitase_aco2_split(self.model.copy(), spec)
        for edit in ("after_gpr", "before_gpr", "aco1_after", "gene_after", "extra_reaction", "missing_reaction"):
            spec = copy.deepcopy(self.spec)
            if edit == "after_gpr":
                spec["reactions"]["R238"]["after_gpr"] = ACO2
            elif edit == "before_gpr":
                spec["reactions"]["R216"]["before_gpr"] = ACO1
            elif edit == "aco1_after":
                spec["genes"][ACO1]["reactions_after"] = ["R239"]
            elif edit == "gene_after":
                spec["genes"][ACO2]["reactions_after"] = ["R238"]
            elif edit == "extra_reaction":
                spec["reactions"]["R239"] = copy.deepcopy(spec["reactions"]["R238"])
            else:
                del spec["reactions"]["R1556"]
            with self.subTest(edit=edit), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_aconitase_aco2_split(self.model.copy(), spec)

    def test_changes_only_the_three_rules_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        for rid in TARGETS:
            target = expected.reactions.get_by_id(rid)
            target.gene_reaction_rule = self.spec["reactions"][rid]["after_gpr"]
            target.notes = {**target.notes, **self.spec["reactions"][rid]["notes"]}
        result = apply_aconitase_aco2_split(model)
        self.assertEqual(result["status"], "applied")
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual([g.id for g in model.genes], [g.id for g in self.model.genes])
        self.assertEqual(list(model.genes.get_by_id(UNLINKED).reactions), [])
        self.assertEqual(model.reactions.R239.gene_reaction_rule, AFTER)
        self.assertEqual(apply_aconitase_aco2_split(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "aconitase.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_aconitase_aco2_split(loaded)["status"], "already_correct")
        # On an applied model the gene guard still holds: a relinked or removed unlinked gene fails closed.
        for case, message in (("relinked", "identity or reaction context"), ("removed", f"{UNLINKED} missing")):
            bad = loaded.copy()
            if case == "relinked":
                bad.reactions.R24.gene_reaction_rule = f"YALI1E03387g or {UNLINKED}"
            else:
                bad.genes.remove(bad.genes.get_by_id(UNLINKED))
            before = state(bad)
            with self.subTest(applied_gene=case), self.assertRaisesRegex(ValueError, message):
                apply_aconitase_aco2_split(bad)
            self.assertEqual(state(bad), before)

    def test_only_the_aco1_gene_closes_the_aconitases(self):
        for curated in (False, True):
            model = self.model.copy()
            if curated:
                apply_aconitase_aco2_split(model)
            for gene_id in (ACO1, ACO2):
                with model, self.subTest(curated=curated, gene=gene_id):
                    model.genes.get_by_id(gene_id).knock_out()
                    closed = {rid for rid in ACONITASES if model.reactions.get_by_id(rid).bounds == (0, 0)}
                    # Before the edit only R239 (ACO1 alone) closes; the ACO2 knockout never closes anything.
                    expected = set() if gene_id == ACO2 else (set(ACONITASES) if curated else {"R239"})
                    self.assertEqual(closed, expected)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_aconitase_aco2_split(model)["status"], "applied")
            self.assertEqual(model.reactions.R216.gene_reaction_rule, AFTER)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_aconitase_aco2_split(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).aconitase_aco2_split)
        self.assertTrue(parse_args(["--aconitase-aco2-split"]).aconitase_aco2_split)


if __name__ == "__main__":
    unittest.main()
