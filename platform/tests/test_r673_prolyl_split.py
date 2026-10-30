"""R673 prolyl split: preservation, rejection, SBML roundtrip, knockout logic and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r673_prolyl_split import AFTER, BEFORE, CURATION_PATH, UNLINKED, apply_r673_prolyl_split
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_um7_nuclear_diffusion_r2171_1to1_aars_thr_ser_"
            "split_r171_atp_synthase_gpr_20261029.xml")
CYTO, MITO = "YALI1E05929g", "YALI1D19110g"


def gene_table(model):
    return {g.id: (copy.deepcopy(g.annotation), sorted(r.id for r in g.reactions)) for g in model.genes}


def state(model):
    return semantics(model), annotations(model), gene_table(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class R673ProlylSplitTests(unittest.TestCase):
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
                    "missing_reaction": "R673 missing"}
        baseline = state(self.model)
        row = self.spec["reaction"]
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R673
            if field == "name":
                reaction.name = "prolyl-tRNA synthetase, mitochondrial"
            elif field == "bounds":
                reaction.bounds = (-1000.0, 1000.0)
            elif field == "gpr":
                reaction.gene_reaction_rule = MITO
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "species":
                bad.metabolites.get_by_id("m765[C_cy]").formula = "C5H8NO2"
            elif field == "annotation":
                reaction.annotation = {**reaction.annotation, "ec-code": "6.1.1.16"}
            elif field == "notes":
                reaction.notes = {**reaction.notes, "r673_prolyl_split_status": "another review"}
            elif field == "extra_note":
                reaction.notes = {**reaction.notes, "unreviewed_note": "added elsewhere"}
            elif field == "notes_without_gpr":
                reaction.notes = {**reaction.notes, **copy.deepcopy(row["notes"])}
            elif field == "after_without_notes":
                reaction.gene_reaction_rule = row["after_gpr"]
            elif field == "gene_identity":
                # Model.copy() shares Gene.annotation dicts with the original: rebind, never mutate them.
                gene = bad.genes.get_by_id(MITO)
                gene.annotation = {**gene.annotation, "ncbigene": "2900000"}
            elif field == "gene_context":
                bad.reactions.R723.gene_reaction_rule = f"YALI1F03908g or {CYTO}"
            else:
                bad.remove_reactions([reaction])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r673_prolyl_split(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        cases = (("status", "user_authorized"), ("assignment_id", "R673-PROLYL-SPLIT"),
                 ("unlinked_gene", {"id": CYTO}))
        for key, value in cases:
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_r673_prolyl_split(self.model.copy(), spec)
        for edit in ("after_gpr", "gene_after"):
            spec = copy.deepcopy(self.spec)
            if edit == "after_gpr":
                spec["reaction"]["after_gpr"] = MITO
            else:
                spec["genes"][MITO]["reactions_after"] = ["R673"]
            with self.subTest(edit=edit), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_r673_prolyl_split(self.model.copy(), spec)

    def test_changes_only_r673_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        target = expected.reactions.R673
        target.gene_reaction_rule = self.spec["reaction"]["after_gpr"]
        target.notes = {**target.notes, **self.spec["reaction"]["notes"]}
        result = apply_r673_prolyl_split(model)
        self.assertEqual(result["status"], "applied")
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual([g.id for g in model.genes], [g.id for g in self.model.genes])
        self.assertEqual(list(model.genes.get_by_id(UNLINKED).reactions), [])
        self.assertEqual(apply_r673_prolyl_split(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r673.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r673_prolyl_split(loaded)["status"], "already_correct")

    def test_only_the_cytosolic_gene_closes_r673(self):
        for curated in (False, True):
            model = self.model.copy()
            if curated:
                apply_r673_prolyl_split(model)
            for gene_id in (CYTO, MITO):
                with model, self.subTest(curated=curated, gene=gene_id):
                    model.genes.get_by_id(gene_id).knock_out()
                    # Before the edit the OR partner keeps R673 open for either single knockout.
                    self.assertEqual(model.reactions.R673.bounds == (0, 0), curated and gene_id == CYTO)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_r673_prolyl_split(model)["status"], "applied")
            self.assertEqual(model.reactions.R673.gene_reaction_rule, AFTER)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertEqual(apply_r673_prolyl_split(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r673_prolyl_split)
        self.assertTrue(parse_args(["--r673-prolyl-split"]).r673_prolyl_split)


if __name__ == "__main__":
    unittest.main()
