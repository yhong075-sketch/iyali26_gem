"""R1172 closure curation: preservation, rejection, SBML roundtrip, context restore and the opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import annotate_retained_reactions, build_model
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r1172_hmgcoa_transport_close import (CURATION_PATH, REQUIRED, ROUTE,
                                                                apply_r1172_hmgcoa_transport_close)
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_"
            "r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_r2121_close_"
            "20261112.xml")


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


def state(model):
    return semantics(model), annotations(model)


class R1172HmgcoaTransportCloseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        differs, missing = "R1172 closure precondition differs", "R1172 or a reaction of the cytosolic route missing"
        route = "is not the cytosolic route reaction this curation reviewed"
        expected = {"name": differs, "gpr": differs, "bounds": differs, "stoichiometry": differs,
                    "species": differs, "compartment": differs, "annotation": differs,
                    "extra_note": "notes outside this curation differ",
                    "missing_retention_note": "notes outside this curation differ",
                    "conflicting_note": "Conflicting R1172 closure notes",
                    "closed_without_notes": "already closed without this curation",
                    "open_with_notes": "open but carries this curation's notes",
                    "r82_gpr": route, "r82_bounds": route, "r82_annotation": route, "r411_gpr": route,
                    "r411_name": route, "r411_stoichiometry": route, "r411_species": route, "r410_name": route,
                    "r410_bounds": route, "r410_annotation": route,
                    "r411_missing": missing, "r1172_missing": missing}
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R1172
            if field == "name":
                reaction.name = "HMG-CoA carrier"
            elif field == "gpr":
                reaction.gene_reaction_rule = "YALI1F20788g"
            elif field == "bounds":
                reaction.bounds = (0.0, 1000.0)
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "species":
                bad.metabolites.get_by_id("m648[C_mi]").formula = "C27H43N7O20P3S"
            elif field == "compartment":
                bad.metabolites.get_by_id("m646[C_cy]").compartment = "C_mi"
            elif field == "annotation":
                reaction.annotation["sbo"] = "SBO:0000176"
            elif field == "extra_note":
                reaction.notes["other_review"] = "something else"
            elif field == "missing_retention_note":
                del reaction.notes["retention_current_status"]
            elif field == "conflicting_note":
                reaction.notes["r1172_hmgcoa_transport_close_curation"] = "conflicting review"
            elif field == "closed_without_notes":
                reaction.bounds = (0.0, 0.0)
            elif field == "open_with_notes":
                reaction.notes.update(self.spec["notes"])
            elif field == "r82_gpr":
                bad.reactions.R82.gene_reaction_rule = "YALI1B11532g or YALI1E13899g"
            elif field == "r82_bounds":
                bad.reactions.R82.bounds = (0.0, 0.0)
            elif field == "r82_annotation":
                bad.reactions.R82.annotation["kegg.reaction"] = "R00391"
            elif field == "r411_gpr":
                bad.reactions.R411.gene_reaction_rule = ""
            elif field == "r411_name":
                bad.reactions.R411.name = "hydroxymethylglutaryl CoA synthase (mitochondrial)"
            elif field == "r411_stoichiometry":
                bad.reactions.R411.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "r411_species":
                # Acetyl-CoA m68 is in R411 and R82 but not in R1172, so this hits the route guard only.
                bad.metabolites.get_by_id("m68[C_cy]").formula = "C23H37N7O17P3S"
            elif field == "r410_name":
                bad.reactions.R410.name = "HMG-CoA reductase (NADH)"
            elif field == "r410_bounds":
                bad.reactions.R410.bounds = (-1000.0, 1000.0)
            elif field == "r410_annotation":
                bad.reactions.R410.annotation["kegg.reaction"] = "R02081"
            elif field == "r411_missing":
                bad.remove_reactions([bad.reactions.R411])
            else:
                bad.remove_reactions([bad.reactions.R1172])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r1172_hmgcoa_transport_close(bad)
            self.assertEqual(state(bad), before)

    def test_unexpected_scope_is_rejected(self):
        for path, value in ((("schema_version",), 2), (("curation_id",), "R1172-HMGCOA-TRANSPORT-CLOSE-20261113"),
                            (("status",), "proposed"), (("reaction_id",), "R412"), (("gpr",), "YALI1F20788g"),
                            (("before_bounds",), [0, 1000]), (("after_bounds",), [0, 1000]),
                            (("required_route",), []), (("required_route", 0, "reaction_id"), "R83"),
                            (("required_route", 1, "gpr"), ""), (("notes",), {}),
                            (("notes", "r1172_previous_bounds"), "[0, 1000]"),
                            (("before_notes",), {"r1172_previous_bounds": "[-1000, 1000]"})):
            spec = copy.deepcopy(self.spec)
            target = spec
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, "Unexpected R1172 closure curation scope"):
                apply_r1172_hmgcoa_transport_close(self.model.copy(), spec)
        for key in REQUIRED:
            spec = copy.deepcopy(self.spec)
            del spec[key]
            with self.subTest(missing=("top", key)), self.assertRaisesRegex(ValueError, "curation scope"):
                apply_r1172_hmgcoa_transport_close(self.model.copy(), spec)
        for index in range(len(ROUTE)):
            for key in (*REQUIRED, "bounds"):
                spec = copy.deepcopy(self.spec)
                del spec["required_route"][index][key]
                with self.subTest(missing=(ROUTE[index], key)), self.assertRaisesRegex(ValueError, "curation scope"):
                    apply_r1172_hmgcoa_transport_close(self.model.copy(), spec)

    def test_sets_only_the_r1172_bounds_and_notes_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        expected.reactions.R1172.bounds = (0.0, 0.0)
        expected.reactions.R1172.notes = {**expected.reactions.R1172.notes, **self.spec["notes"]}
        result = apply_r1172_hmgcoa_transport_close(model)
        self.assertEqual((result["status"], result["bounds"]), ("applied", [0.0, 0.0]))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual(model.reactions.R412.bounds, (0.0, 1000.0))
        self.assertEqual(model.reactions.R1172.gene_reaction_rule, "")
        self.assertIn("unconfirmed", model.reactions.R1172.notes["retention_current_status"])
        self.assertEqual(apply_r1172_hmgcoa_transport_close(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r1172.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(path)
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r1172_hmgcoa_transport_close(loaded)["status"], "already_correct")

    def test_independent_of_the_r2121_closure(self):
        # The step neither reads nor touches R2121: with R2121 open again (the e584151f state) it still applies.
        model = read_baseline()
        model.reactions.R2121.bounds = (-1000.0, 1000.0)
        result = apply_r1172_hmgcoa_transport_close(model)
        self.assertEqual((result["status"], result["bounds"]), ("applied", [0.0, 0.0]))
        self.assertEqual(model.reactions.R2121.bounds, (-1000.0, 1000.0))

    def test_retention_annotation_leaves_the_closure_alone(self):
        # The retention step rewrites only its own notes; the dated reversal field is not written into the model.
        model = read_baseline()
        apply_r1172_hmgcoa_transport_close(model)
        before = state(model)
        records = {r["item"]: r for r in annotate_retained_reactions(model)}
        self.assertEqual(records["R1172"]["status"], "already_correct")
        self.assertIn("later_decision_20261112", records["R1172"])
        self.assertEqual(state(model), before)

    def test_model_context_restores_bounds_and_notes(self):
        # Fresh read: after a context exit on a Model.copy() the solver bound text reads 0 instead of 0.0.
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_r1172_hmgcoa_transport_close(model)["status"], "applied")
            self.assertEqual(model.reactions.R1172.bounds, (0.0, 0.0))
        self.assertEqual(state(model), before)

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r1172_hmgcoa_transport_close)
        self.assertTrue(parse_args(["--r1172-hmgcoa-transport-close"]).r1172_hmgcoa_transport_close)

    def test_build_model_forwards_the_flag(self):
        # A dropped pass-through would make the flag a silent no-op; stop the build at the chain call instead.
        class Stop(Exception):
            pass

        with tempfile.TemporaryDirectory() as directory:
            for argv, expected in (([], False), (["--r1172-hmgcoa-transport-close"], True)):
                args = parse_args([*argv, "--output-model", str(Path(directory) / "never_written.xml")])
                with self.subTest(argv=argv), patch("scripts.gem_annotate.main.build_reference_chain",
                                                    side_effect=Stop) as chain, self.assertRaises(Stop):
                    build_model(args)
                self.assertIs(chain.call_args.kwargs["r1172_hmgcoa_transport_close"], expected)


if __name__ == "__main__":
    unittest.main()
