"""BiGG-ID export: compartment codes, skip list, collisions, IDs-only roundtrip and output guards; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Metabolite, Model, Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tools import export_bigg_ids as ex

REPLACEMENTS = {"C_em": {"from": "em", "to": "erm"}, "C_en": {"from": "en", "to": "ce"}}


def toy_model():
    """Three metabolites, two reactions: one in C_em, one transport C_em -> C_en."""
    model = Model("toy")
    a = Metabolite("m1[C_em]", formula="C2H4O2", charge=0, compartment="C_em", name="acetate")
    b = Metabolite("m2[C_em]", formula="C3H6O3", charge=0, compartment="C_em", name="lactate")
    c = Metabolite("m3[C_en]", formula="C3H6O3", charge=0, compartment="C_en", name="lactate")
    a.annotation = {"bigg.metabolite": "ac", "chebi": "CHEBI:30089"}
    b.annotation = {"bigg.metabolite": "lac__L"}
    r1, r2 = Reaction("R1", name="convert"), Reaction("R2", name="move")
    r1.add_metabolites({a: -1, b: 1})
    r2.add_metabolites({b: -1, c: 1})
    r1.bounds, r2.bounds = (-1000, 1000), (0, 1000)
    r1.gene_reaction_rule = "YALI1A00001g"
    r1.annotation = {"bigg.reaction": ["X1", "R_X1"]}
    r1.notes = {"source": "toy"}
    model.add_reactions([r1, r2])
    # An unnamed in-memory gene comes back from SBML as "G_<id>"; models read from SBML already carry that name.
    model.genes.YALI1A00001g.name = "G_YALI1A00001g"
    model.compartments = {"C_em": "endoplasmic reticulum membrane", "C_en": "cell envelope"}
    return model


def toy_rows(model):
    met = [{"upstream_id": mid, "mcnaugh_id": new, "upstream_name": model.metabolites.get_by_id(mid).name,
            "upstream_formula": model.metabolites.get_by_id(mid).formula,
            "upstream_charge": str(model.metabolites.get_by_id(mid).charge)}
           for mid, new in (("m1[C_em]", "ac_em"), ("m2[C_em]", "lac__L_em"), ("m3[C_en]", "lac__L_en"))]
    rxn = [{"upstream_id": r.id, "mcnaugh_id": new, "status": "identical", "upstream_equation": r.build_reaction_string()}
           for r, new in ((model.reactions.R1, "CONVem"), (model.reactions.R2, "LACt_em_en"))]
    return {"metabolite_crosswalk.csv": met, "reaction_crosswalk.csv": rxn}


POLICY = {"renamed_reaction_statuses": ["identical"], "compartment_code_replacements": REPLACEMENTS, "skip_rename": []}


class ExportBiggIdTests(unittest.TestCase):
    def test_compartment_codes_follow_the_policy(self):
        cases = [("pa_em", {"C_em"}, True, "pa_erm"), ("chol_en", {"C_en"}, True, "chol_ce"),
                 ("pa_p", {"C_pe"}, True, "pa_p"), ("ELON160em", {"C_em"}, False, "ELON160erm"),
                 ("PSt_em_c", {"C_em", "C_cy"}, False, "PSt_erm_c"), ("PCt_em_en", {"C_em", "C_en"}, False, "PCt_erm_ce"),
                 ("FACOAL200en", {"C_en"}, False, "FACOAL200ce"), ("DESAT16", {"C_em"}, False, "DESAT16")]
        for old, compartments, metabolite, new in cases:
            with self.subTest(old=old):
                self.assertEqual(ex.yeast_gem_id(old, compartments, REPLACEMENTS, metabolite), new)
        with self.assertRaisesRegex(ValueError, "without the _em suffix"):
            ex.yeast_gem_id("pa_er", {"C_em"}, REPLACEMENTS, True)

    def test_renames_ids_only_and_records_legacy_ids(self):
        model = toy_model()
        maps, report = ex.rename_plan(model, POLICY, toy_rows(model))
        self.assertEqual(maps["metabolite"], {"m1[C_em]": "ac_erm", "m2[C_em]": "lac__L_erm", "m3[C_en]": "lac__L_ce"})
        self.assertEqual(maps["reaction"], {"R1": "CONVerm", "R2": "LACt_erm_ce"})
        expected = ex.definition(model, maps)
        ex.apply_renames(model, maps)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "toy.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(path)
        self.assertEqual(ex.definition(loaded), expected)
        self.assertEqual(loaded.metabolites.get_by_id("ac_erm").notes["legacy_id"], "m1[C_em]")
        self.assertEqual(loaded.reactions.CONVerm.notes, {"source": "toy", "legacy_id": "R1"})
        self.assertEqual(loaded.reactions.CONVerm.annotation, {"bigg.reaction": ["X1", "R_X1"]})

    def test_skip_list_keeps_legacy_ids(self):
        model = toy_model()
        policy = {**POLICY, "skip_rename": [{"entity": "reaction", "entity_id": "R2"},
                                            {"entity": "metabolite", "entity_id": "m1[C_em]",
                                             "reason": "reviewed BiGG annotation row not written"}]}
        maps, report = ex.rename_plan(model, policy, toy_rows(model))
        self.assertNotIn("R2", maps["reaction"])
        self.assertNotIn("m1[C_em]", maps["metabolite"])
        self.assertIn(("reaction", "R2", "", "kept: reviewed BiGG ID names another entity"), report)
        self.assertIn(("metabolite", "m1[C_em]", "", "kept: reviewed BiGG annotation row not written"), report)

    def test_colliding_target_keeps_legacy_id(self):
        model = toy_model()
        rows = toy_rows(model)
        rows["reaction_crosswalk.csv"][1]["mcnaugh_id"] = "CONVem"
        maps, report = ex.rename_plan(model, POLICY, rows)
        self.assertEqual(list(maps["reaction"]), ["R1"])
        self.assertIn(("reaction", "R2", "", "kept: target CONVerm already used"), report)

    def test_metabolite_suffix_must_match_its_compartment(self):
        model = toy_model()
        rows = toy_rows(model)
        # m3 sits in C_en, so an _em suffix is a crosswalk error, not something to rename.
        rows["metabolite_crosswalk.csv"][2]["mcnaugh_id"] = "lac__L_em"
        with self.assertRaisesRegex(ValueError, "without the _en suffix"):
            ex.rename_plan(model, POLICY, rows)

    def test_stale_crosswalk_or_skip_list_fails_closed(self):
        model = toy_model()
        rows = toy_rows(model)
        rows["reaction_crosswalk.csv"][0]["upstream_equation"] = "m1[C_em] --> m2[C_em]"
        with self.assertRaisesRegex(ValueError, r"to be renamed changed since the crosswalk was built: \['R1'\]"):
            ex.rename_plan(model, POLICY, rows)
        for name, row in (("reaction_crosswalk.csv", {"upstream_id": "R9", "mcnaugh_id": "", "status": "differs",
                                                      "upstream_equation": ""}),
                          ("metabolite_crosswalk.csv", {"upstream_id": "m9[C_em]", "mcnaugh_id": ""})):
            rows = toy_rows(model)
            rows[name].append(row)
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "Crosswalk entities not in the model"):
                ex.rename_plan(model, POLICY, rows)
        policy = {**POLICY, "skip_rename": [{"entity": "reaction", "entity_id": "R9"}]}
        with self.assertRaisesRegex(ValueError, "Skip-list entities not in the crosswalk"):
            ex.rename_plan(model, policy, toy_rows(model))

    def test_entities_newer_than_the_crosswalk_keep_legacy_ids(self):
        model = toy_model()
        rows = toy_rows(model)
        # R2 changed since the crosswalk and is not renamed; R3 and m4 came later and are not in it at all.
        rows["reaction_crosswalk.csv"][1].update(status="differs", upstream_equation="m2[C_em] --> 2.0 m3[C_en]")
        r3 = Reaction("R3")
        r3.add_metabolites({model.metabolites.get_by_id("m3[C_en]"): -1,
                            Metabolite("m4[C_en]", formula="C3H6O3", charge=0, compartment="C_en"): 1})
        model.add_reactions([r3])
        maps, report = ex.rename_plan(model, POLICY, rows)
        self.assertEqual(maps["reaction"], {"R1": "CONVerm"})
        self.assertNotIn("m4[C_en]", maps["metabolite"])
        for row in (("reaction", "R2", "", "kept: equation changed since the crosswalk"),
                    ("reaction", "R3", "", "kept: not in the crosswalk"),
                    ("metabolite", "m4[C_en]", "", "kept: not in the crosswalk")):
            self.assertIn(row, report)

    def test_curated_policy_is_the_reviewed_one(self):
        policy = json.loads(ex.POLICY_PATH.read_text())
        self.assertEqual((policy["policy_id"], policy["status"]), ("BIGG-EXPORT-20261103", "user_authorized"))
        self.assertEqual(policy["compartment_code_replacements"], REPLACEMENTS)
        self.assertEqual(len(policy["skip_rename"]), 31)
        for field, value in (("status", "draft"), ("policy_id", "OTHER")):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "scope/status"):
                ex.load_inputs(Path("."), {**copy.deepcopy(policy), field: value})

    def test_reidentified_metabolite_fails_closed(self):
        for field, value in (("formula", "C3H5O3"), ("charge", -1), ("name", "pyruvate")):
            model = toy_model()
            rows = toy_rows(model)
            setattr(model.metabolites.get_by_id("m2[C_em]"), field, value)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, r"changed name, formula or charge.*m2"):
                ex.rename_plan(model, POLICY, rows)

    def test_refuses_existing_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "export.xml"
            for existing in (output, output.with_name("export.unverified.xml"), output.with_suffix(".rename_report.csv"),
                             output.with_suffix(".export.json")):
                existing.write_text("")
                with self.subTest(existing=existing.name), self.assertRaises(FileExistsError):
                    ex.export(Path("unused.xml"), Path("."), output)
                existing.unlink()


if __name__ == "__main__":
    unittest.main()
