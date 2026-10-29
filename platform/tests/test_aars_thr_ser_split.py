"""aaRS Thr/Ser split: preservation, rejection, partial state, SBML roundtrip, knockout logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.aars_thr_ser_split import CURATION_PATH, apply_aars_thr_ser_split
from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_r2171_1to1_20261028.xml")
THR_CY, THR_MI, SER_CY, SER_MI = "YALI1F14882g", "YALI1A09121g", "YALI1F03908g", "YALI1B07957g"


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class AarsThrSerSplitTests(unittest.TestCase):
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
                    "partial": "partly applied", "gene_identity": "identity or reaction context",
                    "gene_context": "identity or reaction context", "missing_reaction": "R723 missing"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R723
            row = self.spec["reactions"]["R723"]
            if field == "name":
                reaction.name = "seryl-tRNA synthetase, mitochondrial"
            elif field == "bounds":
                reaction.bounds = (-1000.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = SER_MI
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "species":
                bad.metabolites.get_by_id("m441[C_cy]").formula = "C3H6NO3"
            elif field == "annotation":
                reaction.annotation = {**reaction.annotation, "ec-code": "6.1.1.12"}
            elif field == "notes":
                reaction.notes = {**reaction.notes, "aars_thr_ser_split_status": "another review"}
            elif field == "extra_note":
                reaction.notes = {**reaction.notes, "unreviewed_note": "added elsewhere"}
            elif field == "notes_without_gpr":
                reaction.notes = {**reaction.notes, **copy.deepcopy(row["notes"])}
            elif field == "after_without_notes":
                for rid, r in self.spec["reactions"].items():
                    bad.reactions.get_by_id(rid).gene_reaction_rule = r["after_gpr"]
            elif field == "partial":
                reaction.gene_reaction_rule = row["after_gpr"]
                reaction.notes = {**reaction.notes, **copy.deepcopy(row["notes"])}
            elif field == "gene_identity":
                # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
                gene = bad.genes.get_by_id(SER_MI)
                gene.annotation = {**gene.annotation, "ncbigene": "2900000"}
            elif field == "gene_context":
                bad.reactions.R673.gene_reaction_rule = f"YALI1D19110g or {SER_MI}"
            else:
                bad.remove_reactions([reaction])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_aars_thr_ser_split(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        cases = (("status", "user_authorized"), ("assignment_id", "AARS-THR-SER-SPLIT"),
                 ("unlinked_gene", {"id": THR_MI}))
        for key, value in cases:
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_aars_thr_ser_split(self.model.copy(), spec)
        spec = copy.deepcopy(self.spec)
        spec["reactions"]["R759"]["after_gpr"] = THR_MI
        with self.assertRaisesRegex(ValueError, "scope/status"):
            apply_aars_thr_ser_split(self.model.copy(), spec)
        spec = copy.deepcopy(self.spec)
        spec["reactions"]["R673"] = copy.deepcopy(spec["reactions"]["R723"])
        with self.assertRaisesRegex(ValueError, "scope/status"):
            apply_aars_thr_ser_split(self.model.copy(), spec)

    def test_changes_only_the_three_rules_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        for rid, row in self.spec["reactions"].items():
            target = expected.reactions.get_by_id(rid)
            target.gene_reaction_rule = row["after_gpr"]
            target.notes = {**target.notes, **row["notes"]}
        result = apply_aars_thr_ser_split(model)
        self.assertEqual(result["status"], "applied")
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual([g.id for g in model.genes], [g.id for g in self.model.genes])
        self.assertEqual(list(model.genes.get_by_id(SER_MI).reactions), [])
        self.assertEqual(model.reactions.R673.gene_reaction_rule, "YALI1D19110g or YALI1E05929g")
        self.assertEqual(apply_aars_thr_ser_split(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "aars.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_aars_thr_ser_split(loaded)["status"], "already_correct")

    def test_each_gene_closes_only_its_own_reaction(self):
        model = self.model.copy()
        apply_aars_thr_ser_split(model)
        cases = {THR_CY: {"R759"}, THR_MI: {"R760"}, SER_CY: {"R723"}, SER_MI: set()}
        for gene_id, closed in cases.items():
            with model, self.subTest(gene=gene_id):
                model.genes.get_by_id(gene_id).knock_out()
                shut = {rid for rid in ("R759", "R760", "R723") if model.reactions.get_by_id(rid).bounds == (0, 0)}
                self.assertEqual(shut, closed)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_aars_thr_ser_split(model)["status"], "applied")
            self.assertEqual(model.reactions.R759.gene_reaction_rule, THR_CY)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_aars_thr_ser_split(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).aars_thr_ser_split)
        self.assertTrue(parse_args(["--aars-thr-ser-split"]).aars_thr_ser_split)


if __name__ == "__main__":
    unittest.main()
