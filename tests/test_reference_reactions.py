"""Focused checks for the two explicitly retained reference-build reactions."""

import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra import Model, Metabolite, Reaction

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.gaps import add_gap_fill_reactions
from scripts.gem_annotate.locus_resolver import LocusCrosswalk
from scripts.gem_annotate.main import annotate_retained_reactions
from scripts.gem_annotate.patches import remove_spurious_transport_reactions
from tests.test_coq9_curation import annotations, semantics


class ReferenceReactionTests(unittest.TestCase):
    def test_public_defaults_keep_reference_input_and_metadata(self):
        args = parse_args([])
        self.assertEqual(args.starting_model.name, "iyali26.xml")
        self.assertEqual(args.coq9_curation, "metadata")
        self.assertIsNone(args.provisional_capacity_profile)
        self.assertIsNone(args.r608_curation)

    def test_retention_survives_removal_and_duplicate_filters_idempotently(self):
        model = Model("retained_reactions")
        a = Metabolite("a", formula="C", charge=0, compartment="C_cy")
        b = Metabolite("b", formula="C", charge=0, compartment="C_cy")
        a.annotation = {"metanetx.chemical": "MNXM100"}
        b.annotation = {"metanetx.chemical": "MNXM101"}
        for rid in ("R1172", "R730", "R663"):
            reaction = Reaction(rid)
            reaction.add_metabolites({a: -1, b: 1})
            model.add_reactions([reaction])
        model.reactions.R730.annotation["bigg.reaction"] = "SPHPL"
        before = semantics(model)
        self.assertEqual(remove_spurious_transport_reactions(model), 0)
        self.assertEqual(before, semantics(model))
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            first, second = root / "crosswalk.csv", root / "metabolic.csv"
            first.write_text("model_gene,yali1_s2,yali0,category,n_reactions\n")
            second.write_text("yali1,yali0,geneid,metab_pathways,ec,verdict,in_model\n")
            resolver = LocusCrosswalk.from_csvs(first, second)
            candidates = root / "candidates.csv"
            candidates.write_text(
                "priority,mnxr_id,bigg_reaction,gene_id,equation,ec_number,kegg_reaction\n"
                'P0,MNXR188844,SPHPL,,"1 MNXM100@MNXD1 = 1 MNXM101@MNXD1",,\n'
                'P0,MNXR146152,R_PSPHPL,,"1 MNXM100@MNXD1 = 1 MNXM101@MNXD1",,\n'
            )
            with patch("scripts.gem_annotate.gaps.load_default_locus_crosswalk", return_value=resolver):
                first_run = add_gap_fill_reactions(model, candidates, cache_dir=root / "cache")
                self.assertIn("SPHPL", model.reactions)
                self.assertEqual(first_run["skipped_curated_existing"], ["R_PSPHPL"])
                state = semantics(model)
                add_gap_fill_reactions(model, candidates, cache_dir=root / "cache")
                self.assertEqual(state, semantics(model))
        annotate_retained_reactions(model)
        state = copy.deepcopy(annotations(model))
        self.assertTrue(all(r["status"] == "already_correct" for r in annotate_retained_reactions(model)))
        self.assertEqual(state, annotations(model))
        self.assertEqual(model.reactions.R1172.gene_reaction_rule, "")
        self.assertIn("unconfirmed", model.reactions.R1172.notes["retention_current_status"])
        self.assertIn("R730", model.reactions.SPHPL.notes["retention_previous_policy"])


if __name__ == "__main__":
    unittest.main()
