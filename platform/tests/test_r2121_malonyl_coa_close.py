"""R2121 closure curation: preservation, rejection, SBML roundtrip, context restore and the opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_model
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r2121_malonyl_coa_close import CURATION_PATH, REQUIRED, apply_r2121_malonyl_coa_close
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_"
            "r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_20261103.xml")


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


def state(model):
    return semantics(model), annotations(model)


class R2121MalonylCoaCloseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        differs, kept = "R2121 closure precondition differs", "R2119 is not the retained reaction"
        missing = "R2121 or R2119 missing"
        expected = {"name": differs, "gpr": differs, "bounds": differs, "stoichiometry": differs,
                    "species": differs, "compartment": differs, "annotation": differs,
                    "extra_note": "notes outside this curation differ",
                    "conflicting_note": "Conflicting R2121 closure notes",
                    "closed_without_notes": "already closed without this curation",
                    "open_with_notes": "open but carries this curation's notes",
                    "r2119_bounds": kept, "r2119_gpr": kept, "r2119_name": kept, "r2119_stoichiometry": kept,
                    "r2119_annotation": kept, "r2119_species": kept,
                    "r2119_missing": missing, "r2121_missing": missing}
        for field, message in expected.items():
            bad = self.model.copy()
            reaction = bad.reactions.R2121
            if field == "name":
                reaction.name = "malonate-semialdehyde dehydrogenase"
            elif field == "gpr":
                reaction.gene_reaction_rule = "YALI1C02550g or YALI1E21915g"
            elif field == "bounds":
                reaction.bounds = (0.0, 1000.0)
            elif field == "stoichiometry":
                reaction.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "species":
                bad.metabolites.get_by_id("m1976[C_cy]").formula = "C3H4O3"
            elif field == "compartment":
                bad.metabolites.get_by_id("m200[C_cy]").compartment = "C_mi"
            elif field == "annotation":
                reaction.annotation["kegg.reaction"] = "R00705"
            elif field == "extra_note":
                reaction.notes["other_review"] = "something else"
            elif field == "conflicting_note":
                reaction.notes["r2121_malonyl_coa_close_curation"] = "conflicting review"
            elif field == "closed_without_notes":
                reaction.bounds = (0.0, 0.0)
            elif field == "open_with_notes":
                reaction.notes.update(self.spec["notes"])
            elif field == "r2119_bounds":
                bad.reactions.R2119.bounds = (-1000.0, 1000.0)
            elif field == "r2119_gpr":
                bad.reactions.R2119.gene_reaction_rule = "YALI1E21915g"
            elif field == "r2119_name":
                bad.reactions.R2119.name = "malonate-semialdehyde dehydrogenase"
            elif field == "r2119_stoichiometry":
                bad.reactions.R2119.add_metabolites({bad.metabolites.get_by_id("m10[C_cy]"): 1})
            elif field == "r2119_annotation":
                bad.reactions.R2119.annotation["kegg.reaction"] = "R00740"
            elif field == "r2119_species":
                # NADH m123 is in R2119 but not in R2121, so this hits the R2119 guard only.
                bad.metabolites.get_by_id("m123[C_cy]").formula = "C21H30N7O14P2"
            elif field == "r2119_missing":
                bad.remove_reactions([bad.reactions.R2119])
            else:
                bad.remove_reactions([bad.reactions.R2121])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r2121_malonyl_coa_close(bad)
            self.assertEqual(state(bad), before)

    def test_unexpected_scope_is_rejected(self):
        for path, value in ((("schema_version",), 2), (("curation_id",), "R2121-MALONYL-COA-CLOSE-20261113"),
                            (("status",), "proposed"), (("reaction_id",), "R2119"),
                            (("before_bounds",), [0, 1000]), (("after_bounds",), [0, 1000]),
                            (("retained_reaction", "reaction_id"), "R2112"),
                            (("retained_reaction", "gpr"), "YALI1E21915g"),
                            (("retained_reaction", "bounds"), [-1000, 1000]), (("notes",), {}),
                            (("notes", "r2121_previous_bounds"), "[0, 1000]"),
                            (("before_notes",), {"r2121_previous_bounds": "[-1000, 1000]"})):
            spec = copy.deepcopy(self.spec)
            target = spec
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, "Unexpected R2121 closure curation scope"):
                apply_r2121_malonyl_coa_close(self.model.copy(), spec)
        for part in ("top", "retained_reaction"):
            for key in REQUIRED:
                spec = copy.deepcopy(self.spec)
                del (spec if part == "top" else spec["retained_reaction"])[key]
                with self.subTest(missing=(part, key)), self.assertRaisesRegex(ValueError, "curation scope"):
                    apply_r2121_malonyl_coa_close(self.model.copy(), spec)

    def test_sets_only_the_r2121_bounds_and_notes_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model, expected = read_baseline(), read_baseline()
        expected.reactions.R2121.bounds = (0.0, 0.0)
        expected.reactions.R2121.notes = {**expected.reactions.R2121.notes, **self.spec["notes"]}
        result = apply_r2121_malonyl_coa_close(model)
        self.assertEqual((result["status"], result["bounds"]), ("applied", [0.0, 0.0]))
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual(model.reactions.R2119.bounds, (0.0, 1000.0))
        self.assertEqual(model.reactions.R2121.gene_reaction_rule, "YALI1C02550g")
        self.assertEqual(apply_r2121_malonyl_coa_close(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r2121.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(path)
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r2121_malonyl_coa_close(loaded)["status"], "already_correct")

    def test_model_context_restores_bounds_and_notes(self):
        # Fresh read: after a context exit on a Model.copy() the solver bound text reads 0 instead of 0.0.
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_r2121_malonyl_coa_close(model)["status"], "applied")
            self.assertEqual(model.reactions.R2121.bounds, (0.0, 0.0))
        self.assertEqual(state(model), before)

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).r2121_malonyl_coa_close)
        self.assertTrue(parse_args(["--r2121-malonyl-coa-close"]).r2121_malonyl_coa_close)

    def test_build_model_forwards_the_flag(self):
        # A dropped pass-through would make the flag a silent no-op; stop the build at the chain call instead.
        class Stop(Exception):
            pass

        with tempfile.TemporaryDirectory() as directory:
            for argv, expected in (([], False), (["--r2121-malonyl-coa-close"], True)):
                args = parse_args([*argv, "--output-model", str(Path(directory) / "never_written.xml")])
                with self.subTest(argv=argv), patch("scripts.gem_annotate.main.build_reference_chain",
                                                    side_effect=Stop) as chain, self.assertRaises(Stop):
                    build_model(args)
                self.assertIs(chain.call_args.kwargs["r2121_malonyl_coa_close"], expected)


if __name__ == "__main__":
    unittest.main()
