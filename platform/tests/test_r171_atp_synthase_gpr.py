"""R171 ATP synthase GPR: preservation, rejection, SBML roundtrip, knockout logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r171_atp_synthase_gpr import (AFTER, BEFORE, CORE, CURATION_PATH, UNLINKED,
                                                        apply_r171_atp_synthase_gpr)
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_aars_thr_ser_"
            "split_20261029.xml")
SUBUNIT_E, SUBUNIT_J = "YALI1E38131g", "YALI0D17490g"


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class R171AtpSynthaseGprTests(unittest.TestCase):
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
                    "gene_identity": "identity or reaction context", "gene_context": "identity or reaction context",
                    "missing_reaction": "R171 missing"}
        baseline = state(self.model)
        row = self.spec["reaction"]
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R171
            if field == "name":
                reaction.name = "ATP synthase, mitochondrial"
            elif field == "bounds":
                reaction.bounds = (-1000.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = BEFORE.replace(" or YALI1B15733g", "")
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): -1})
            elif field == "species":
                bad.metabolites.get_by_id("m26[C_mi]").formula = "C10H12N5O13P3"
            elif field == "annotation":
                reaction.annotation = {**reaction.annotation, "ec-code": "7.1.2.1"}
            elif field == "notes":
                reaction.notes = {**reaction.notes, "r171_atp_synthase_gpr_status": "another review"}
            elif field == "extra_note":
                reaction.notes = {**reaction.notes, "unreviewed_note": "added elsewhere"}
            elif field == "notes_without_gpr":
                reaction.notes = {**reaction.notes, **copy.deepcopy(row["notes"])}
            elif field == "after_without_notes":
                reaction.gene_reaction_rule = row["after_gpr"]
            elif field == "gene_identity":
                # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
                gene = bad.genes.get_by_id(SUBUNIT_E)
                gene.annotation = {**gene.annotation, "uniprot": "Q00000"}
            elif field == "gene_context":
                bad.reactions.R172.gene_reaction_rule = SUBUNIT_E
            else:
                bad.remove_reactions([reaction])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r171_atp_synthase_gpr(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        cases = (("status", "user_authorized"), ("assignment_id", "R171-ATP-SYNTHASE-GPR"),
                 ("unlinked_genes", {g: {} for g in UNLINKED[:3]}))
        for key, value in cases:
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_r171_atp_synthase_gpr(self.model.copy(), spec)
        for edit in ("after_gpr", "gene_after", "gene_extra"):
            spec = copy.deepcopy(self.spec)
            if edit == "after_gpr":
                spec["reaction"]["after_gpr"] = f"{AFTER} and {SUBUNIT_J}"
            elif edit == "gene_after":
                spec["genes"][SUBUNIT_J]["reactions_after"] = ["R171"]
            else:
                spec["genes"]["YALI1F07747g"] = copy.deepcopy(spec["genes"][SUBUNIT_E])
            with self.subTest(edit=edit), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_r171_atp_synthase_gpr(self.model.copy(), spec)

    def test_changes_only_r171_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        target = expected.reactions.R171
        target.gene_reaction_rule = self.spec["reaction"]["after_gpr"]
        target.notes = {**target.notes, **self.spec["reaction"]["notes"]}
        result = apply_r171_atp_synthase_gpr(model)
        self.assertEqual(result["status"], "applied")
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual([g.id for g in model.genes], [g.id for g in self.model.genes])
        self.assertEqual({g.id for g in model.reactions.R171.genes}, set(CORE))
        for gene_id in UNLINKED:
            self.assertEqual(list(model.genes.get_by_id(gene_id).reactions), [])
        self.assertEqual(apply_r171_atp_synthase_gpr(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r171.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r171_atp_synthase_gpr(loaded)["status"], "already_correct")

    def test_each_core_gene_closes_r171_and_unlinked_genes_close_nothing(self):
        for curated in (False, True):
            model = self.model.copy()
            if curated:
                apply_r171_atp_synthase_gpr(model)
            for gene_id in CORE + UNLINKED:
                with model, self.subTest(curated=curated, gene=gene_id):
                    model.genes.get_by_id(gene_id).knock_out()
                    # Before the edit the g / k branches keep R171 open for every single knockout.
                    self.assertEqual(model.reactions.R171.bounds == (0, 0), curated and gene_id in CORE)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_r171_atp_synthase_gpr(model)["status"], "applied")
            self.assertEqual(model.reactions.R171.gene_reaction_rule, AFTER)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_r171_atp_synthase_gpr(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r171_atp_synthase_gpr)
        self.assertTrue(parse_args(["--r171-atp-synthase-gpr"]).r171_atp_synthase_gpr)


if __name__ == "__main__":
    unittest.main()
