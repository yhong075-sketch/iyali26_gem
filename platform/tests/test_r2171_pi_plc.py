"""R2171 PI phospholipase C written 1:1: preservation, rejection, SBML roundtrip and the opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.r2171_pi_plc import AFTER, CURATION_PATH, apply_r2171_pi_plc_1to1
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_20261028.xml")


def load():
    # A Path, not a str: cobra applies the 260-character path test on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class R2171PiPlcTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = load()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def expected(self):
        model = load()
        reaction = model.reactions.R2171
        reaction.add_metabolites({model.metabolites.get_by_id(m): c for m, c in AFTER.items()}, combine=False)
        reaction.notes.update(self.spec["notes"])
        return model

    def test_stale_inputs_fail_closed_without_edits(self):
        species = "A species of the R2171 curation differs"
        expected = {"name": "R2171 precondition differs", "gpr": "R2171 precondition differs",
                    "bounds": "R2171 precondition differs", "kegg_link": "R2171 precondition differs",
                    "stoichiometry": "R2171 stoichiometry is neither", "after_without_notes": "R2171 stoichiometry is neither",
                    "partial_notes": "R2171 stoichiometry is neither", "notes": "Conflicting R2171 curation notes",
                    "after_conflicting_notes": "Conflicting R2171 curation notes", "water_formula": species,
                    "dag_identity": species,
                    "notes_without_equation": "R2171 carries this curation's notes without its equation",
                    "ins1p_charge": species, "ins1p_formula": species, "ins1p_identity": species, "pi_charge": species,
                    "pi_identity": species, "dag_compartment": species, "proton_charge": species,
                    "missing_species": "R2171 or one of its species missing", "missing_reaction": "R2171 or one of its species missing"}
        for field, message in expected.items():
            bad = self.model.copy()
            r2171 = bad.reactions.R2171
            if field == "name":
                r2171.name = "phosphatidylinositol phospholipase C"
            elif field == "gpr":
                r2171.gene_reaction_rule = "YALI1D13510g or YALI1C05570g"
            elif field == "bounds":
                r2171.bounds = (-1000.0, 1000.0)
            elif field == "kegg_link":
                r2171.annotation = {**r2171.annotation, "kegg.reaction": "R03364"}
            elif field == "stoichiometry":
                r2171.add_metabolites({bad.metabolites.get_by_id("m813[C_cy]"): 50.0}, combine=False)
            elif field in ("after_without_notes", "partial_notes"):
                r2171.add_metabolites({bad.metabolites.get_by_id(m): c for m, c in AFTER.items()}, combine=False)
                if field == "partial_notes":
                    r2171.notes.update(dict(list(self.spec["notes"].items())[:2]))
            elif field == "notes":
                r2171.notes["r2171_pi_plc_curation"] = "conflicting review"
            elif field == "after_conflicting_notes":
                r2171.add_metabolites({bad.metabolites.get_by_id(m): c for m, c in AFTER.items()}, combine=False)
                r2171.notes.update(self.spec["notes"])
                r2171.notes["r2171_pi_plc_limits"] = "conflicting review"
            elif field == "water_formula":
                bad.metabolites.get_by_id("m32[C_cy]").formula = "H2O2"
            elif field == "dag_identity":
                met = bad.metabolites.get_by_id("m1994[C_cy]")
                met.annotation = {**met.annotation, "metanetx.chemical": "MNXM1364164"}
            elif field == "notes_without_equation":
                r2171.notes.update(self.spec["notes"])
            elif field == "ins1p_charge":
                bad.metabolites.get_by_id("m813[C_cy]").charge = -1
            elif field == "ins1p_formula":
                bad.metabolites.get_by_id("m813[C_cy]").formula = "C6H13O9P"
            elif field == "ins1p_identity":
                # Replace, not mutate: Model.copy() shares each metabolite's annotation dict with the fixture.
                met = bad.metabolites.get_by_id("m813[C_cy]")
                met.annotation = {**met.annotation, "kegg.compound": "C04299"}
            elif field == "pi_charge":
                bad.metabolites.get_by_id("m1393[C_cy]").charge = -1
            elif field == "pi_identity":
                met = bad.metabolites.get_by_id("m1393[C_cy]")
                met.annotation = {**met.annotation, "metanetx.chemical": "MNXM1091588"}
            elif field == "dag_compartment":
                bad.metabolites.get_by_id("m1994[C_cy]").compartment = "C_er"
            elif field == "proton_charge":
                bad.metabolites.get_by_id("m10[C_cy]").charge = 0
            elif field == "missing_species":
                bad.remove_metabolites([bad.metabolites.get_by_id("m813[C_cy]")])
            else:
                bad.remove_reactions([r2171])
            before = semantics(bad), annotations(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_r2171_pi_plc_1to1(bad)
            self.assertEqual((semantics(bad), annotations(bad)), before)

    def test_unexpected_scope_is_rejected(self):
        after = dict(self.spec["after_stoichiometry"], **{"m10[C_cy]": 1.0})
        before = dict(self.spec["before_stoichiometry"], **{"m35[C_cy]": 1.0})
        species = {k: v for k, v in self.spec["species"].items() if k != "m10[C_cy]"}
        arrow = dict(self.spec["notes"], r2171_previous_equation="PI + 100 H2O -> DAG")
        for key, value in (("schema_version", 2), ("curation_id", "R2171-PI-PLC-1TO1-20261029"),
                           ("status", "proposed"), ("reaction_id", "R2145"), ("after_stoichiometry", after),
                           ("before_stoichiometry", before), ("species", species), ("notes", {}), ("notes", arrow)):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "Unexpected R2171 PI phospholipase C curation scope"):
                apply_r2171_pi_plc_1to1(self.model.copy(), spec)

    def test_changes_only_r2171_balances_charge_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = load()
        expected = self.expected()
        self.assertEqual(apply_r2171_pi_plc_1to1(model)["status"], "applied")
        self.assertTrue((semantics(model), annotations(model)) ==
                        (semantics(expected), annotations(expected)), "Unexpected model change")
        self.assertEqual(model.reactions.R2171.build_reaction_string(),
                         "m1393[C_cy] + m32[C_cy] --> 2.0 m10[C_cy] + m1994[C_cy] + m813[C_cy]")
        self.assertEqual(sum(m.charge * c for m, c in model.reactions.R2171.metabolites.items()), 0)
        self.assertEqual(apply_r2171_pi_plc_1to1(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r2171.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue((semantics(loaded), annotations(loaded)) ==
                        (semantics(expected), annotations(expected)), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_r2171_pi_plc_1to1(loaded)["status"], "already_correct")

    def test_flag_is_opt_in_and_requires_the_dolichol_route(self):
        self.assertFalse(parse_args([]).r2171_pi_plc_1to1)
        self.assertTrue(parse_args(["--r2171-pi-plc-1to1"]).r2171_pi_plc_1to1)
        with self.assertRaisesRegex(ValueError, "requires --dolichol-route"):
            build_reference_chain(r2171_pi_plc_1to1=True)
        with self.assertRaisesRegex(ValueError, "requires --dolichol-route"):
            build_reference_chain(r2171_pi_plc_1to1=True, lipoylation_module=True, cofactor_biomass_terms=True,
                                  r2135_acyl_coa_pool=True, r1726_cpt_open=True)


if __name__ == "__main__":
    unittest.main()
