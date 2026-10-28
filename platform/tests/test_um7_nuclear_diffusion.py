"""UM7 nuclear diffusion curation: preservation, rejection, partial states, SBML roundtrip, exact edits and the flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.core import Group
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from scripts.gem_annotate.um7_nuclear_diffusion import CURATION_PATH, apply_um7_nuclear_diffusion
from tests.test_coq9_curation import annotations, semantics
from tests.test_gpi_chain_links import group_table

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_r2135_pool_r1726_cpt_dolichol_route_v2_20261028.xml")
# Written out here, independently of the curation file, so a changed file cannot pass silently.
LINKS = {"R_ADPRtn": {"m1725[C_cy]": -1.0, "m818[C_nu]": 1.0},
         "R_PPItn": {"m203[C_cy]": -1.0, "m210[C_nu]": 1.0}}


def read_baseline():
    # A Path, not a str: cobra rejects str paths of 260+ characters on any platform whose name contains "win",
    # which includes darwin.
    return read_sbml_model(MODEL.candidate_file(BASELINE))


def state(model):
    return semantics(model), annotations(model), group_table(model)


def expected_model(spec):
    """The baseline with the two reactions added by hand."""
    model = read_baseline()
    group = next(g for g in model.groups if g.name == "Transport, nucleus")
    rows = {row["reaction_id"]: row for row in spec["new_reactions"]}
    for rid, stoichiometry in LINKS.items():
        reaction = Reaction(rid, name=rows[rid]["reaction_name"], subsystem="Transport, nucleus",
                            lower_bound=-1000.0, upper_bound=1000.0)
        model.add_reactions([reaction])
        reaction.add_metabolites({model.metabolites.get_by_id(m): c for m, c in sorted(stoichiometry.items())})
        reaction.annotation = {"sbo": "SBO:0000185"}
        reaction.notes = dict(rows[rid]["notes"])
        group.add_members([reaction])
    return model


class UM7NuclearDiffusionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_spec_matches_the_task(self):
        spec = self.spec
        self.assertEqual({row["reaction_id"]: row["stoichiometry"] for row in spec["new_reactions"]}, LINKS)
        self.assertTrue(all(row["gpr"] == "" and row["bounds"] == [-1000.0, 1000.0] for row in spec["new_reactions"]))
        self.assertEqual(spec["baseline_model_sha256"],
                         "9d655c2bdbf82721a38f148a9a4ee0a385bb509ea7b591111107113c4aa8335c")
        # Notes are written to SBML; these characters are escaped or dropped there.
        text = json.dumps([row["notes"] for row in spec["new_reactions"]])
        self.assertFalse(set("<>&") & set(text))

    def test_stale_inputs_fail_closed_without_edits(self):
        cases = {
            "adpr_cy_charge": (lambda m: setattr(m.metabolites.get_by_id("m1725[C_cy]"), "charge", -3), "m1725\\[C_cy\\] identity differs"),
            "adpr_nu_formula": (lambda m: setattr(m.metabolites.get_by_id("m818[C_nu]"), "formula", "C15H23N5O14P2"), "m818\\[C_nu\\] identity differs"),
            "ppi_cy_compartment": (lambda m: setattr(m.metabolites.get_by_id("m203[C_cy]"), "compartment", "C_nu"), "m203\\[C_cy\\] identity differs"),
            "ppi_nu_name": (lambda m: setattr(m.metabolites.get_by_id("m210[C_nu]"), "name", "pyrophosphate"), "m210\\[C_nu\\] identity differs"),
            "adpr_nu_new_consumer": (self._add_consumer("m818[C_nu]"), "m818\\[C_nu\\] reaction context differs"),
            "adpr_cy_new_producer": (self._add_consumer("m1725[C_cy]", 1), "m1725\\[C_cy\\] reaction context differs"),
            "ppi_nu_new_consumer": (self._add_consumer("m210[C_nu]"), "m210\\[C_nu\\] reaction context differs"),
            "species_missing": (lambda m: m.remove_metabolites([m.metabolites.get_by_id("m818[C_nu]")]), "m818\\[C_nu\\] missing"),
            "group_missing": (lambda m: m.remove_groups([next(g for g in m.groups if g.name == "Transport, nucleus")]), "Group 'Transport, nucleus' is missing"),
            "group_not_unique": (lambda m: m.add_groups([Group("G_TEST_DUP", name="Transport, nucleus")]), "Group 'Transport, nucleus' is missing or not unique"),
            # Any other link of a pair touches a nuclear end, so the context guard catches it.
            "ppi_pair_already_linked": (self._link("m203[C_cy]", "m210[C_nu]"), "m210\\[C_nu\\] reaction context differs"),
            "new_id_taken": (self._other_definition("R_PPItn"), "R_PPItn already exists with another definition"),
        }
        for field, (mutate, message) in cases.items():
            bad = self.model.copy()
            mutate(bad)
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_um7_nuclear_diffusion(bad)
            self.assertEqual(state(bad), before, field)

    @staticmethod
    def _add_consumer(mid, coefficient=-1):
        def mutate(model):
            reaction = Reaction("R_TEST_SINK", lower_bound=0.0, upper_bound=1000.0)
            model.add_reactions([reaction])
            reaction.add_metabolites({model.metabolites.get_by_id(mid): coefficient})
        return mutate

    @staticmethod
    def _link(cytosol, nucleus):
        def mutate(model):
            reaction = Reaction("R_TEST_LINK", lower_bound=-1000.0, upper_bound=1000.0)
            model.add_reactions([reaction])
            reaction.add_metabolites({model.metabolites.get_by_id(cytosol): -1, model.metabolites.get_by_id(nucleus): 1})
        return mutate

    @staticmethod
    def _other_definition(rid):
        def mutate(model):
            reaction = Reaction(rid, lower_bound=0.0, upper_bound=1000.0)
            model.add_reactions([reaction])
            reaction.add_metabolites({model.metabolites.get_by_id("m203[C_cy]"): -1,
                                      model.metabolites.get_by_id("m210[C_nu]"): 1})
        return mutate

    def test_partial_states_are_rejected(self):
        full = expected_model(self.spec)
        for rid in LINKS:
            model = full.copy()
            model.remove_reactions([model.reactions.get_by_id(rid)])
            before = state(model)
            with self.subTest(removed=rid), self.assertRaisesRegex(ValueError, "Only part of the UM7 nuclear diffusion"):
                apply_um7_nuclear_diffusion(model)
            self.assertEqual(state(model), before, rid)

    def test_unexpected_scope_is_rejected(self):
        for path, value in ((("schema_version",), 2), (("curation_id",), "UM7-NUCLEAR-DIFFUSION-20261029"),
                            (("status",), "proposed"), (("group",), "Transport, mitochondrion"),
                            (("new_reactions", 0, "reaction_id"), "R_ADPRtm"),
                            (("new_reactions", 1, "gpr"), "YALI1F18105g"),
                            (("new_reactions", 0, "bounds"), [0.0, 1000.0]),
                            (("new_reactions", 1, "notes"), {}),
                            (("new_reactions", 0, "stoichiometry", "m1725[C_cy]"), -2.0),
                            # A species without an identity row would be looked up only during the edit.
                            (("new_reactions", 1, "stoichiometry", "m1[C_cy]"), 1.0)):
            spec = copy.deepcopy(self.spec)
            target = spec
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, "Unexpected UM7 nuclear diffusion curation scope"):
                apply_um7_nuclear_diffusion(self.model.copy(), spec)

    def test_makes_exactly_the_two_reactions_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_baseline()
        expected = expected_model(self.spec)
        self.assertEqual(apply_um7_nuclear_diffusion(model)["status"], "applied")
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        for rid, stoichiometry in LINKS.items():
            reaction = model.reactions.get_by_id(rid)
            self.assertEqual({m.id: c for m, c in reaction.metabolites.items()}, stoichiometry)
            self.assertEqual((reaction.bounds, reaction.gene_reaction_rule), ((-1000.0, 1000.0), ""))
            self.assertEqual(reaction.check_mass_balance(), {})
        self.assertEqual(apply_um7_nuclear_diffusion(model)["status"], "already_correct")
        self.assertTrue(state(model) == state(expected), "Rerun changed the model")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "um7.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_um7_nuclear_diffusion(loaded)["status"], "already_correct")

    def test_flag_is_opt_in_and_requires_the_dolichol_route(self):
        self.assertFalse(parse_args([]).um7_nuclear_diffusion)
        self.assertTrue(parse_args(["--um7-nuclear-diffusion"]).um7_nuclear_diffusion)
        with self.assertRaisesRegex(ValueError, "requires --dolichol-route"):
            build_reference_chain(um7_nuclear_diffusion=True)


if __name__ == "__main__":
    unittest.main()
