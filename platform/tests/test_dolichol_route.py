"""Dolichol route curation: preservation, rejection, partial states, SBML roundtrip, exact edits and the flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.core import Group
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.dolichol_route import CURATION_PATH, apply_dolichol_route
from scripts.gem_annotate.main import build_reference_chain
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics
from tests.test_gpi_chain_links import group_table

BASELINE = ("E5_coq9_resp_pkg_r132_combined_gpi_tpp_thiamine_p_lock_heme_iron_r2193_gpr_lipoylation_acylacp_"
            "cofactor_biomass_r2135_pool_r1726_cpt_20261027.xml")
NEW, TERM = "R_DOLSYN", "m503[C_er]"
# Written out here, independently of the curation file, so a changed file cannot pass silently.
NEW_STOICHIOMETRY = {"m2019[C_cy]": -1.0, "m32[C_cy]": -1.0, "m122[C_cy]": -1.0, "m40[C_cy]": -2.0, "m10[C_cy]": -2.0,
                     "m499[C_cy]": 1.0, "m203[C_cy]": 1.0, "m123[C_cy]": 1.0, "m41[C_cy]": 2.0}


def state(model):
    return semantics(model), annotations(model), group_table(model)


def expected_model(spec):
    """The baseline with the three edits made by hand."""
    model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
    model.reactions.R2178.bounds = (0.0, 0.0)
    model.reactions.R2178.notes.update(spec["closed_reaction"]["notes"])
    reaction = Reaction(NEW, name=spec["new_reaction"]["reaction_name"], subsystem="N-Glycan biosynthesis",
                        lower_bound=0.0, upper_bound=1000.0)
    model.add_reactions([reaction])
    reaction.add_metabolites({model.metabolites.get_by_id(m): c for m, c in sorted(NEW_STOICHIOMETRY.items())})
    reaction.annotation = {"sbo": "SBO:0000176"}
    reaction.notes = dict(spec["new_reaction"]["notes"])
    next(g for g in model.groups if g.name == "N-Glycan biosynthesis").add_members([reaction])
    model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id(TERM): -1e-4})
    model.reactions.biomass_C.notes.update(spec["biomass_term"]["notes"])
    return model


class DolicholRouteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_spec_matches_the_task(self):
        spec = self.spec
        self.assertEqual(spec["new_reaction"]["stoichiometry"], NEW_STOICHIOMETRY)
        self.assertEqual(spec["biomass_term"]["coefficient_mmol_per_gDW"], 1e-4)
        self.assertEqual(spec["closed_reaction"]["after_bounds"], [0.0, 0.0])
        self.assertEqual(spec["new_reaction"]["gpr"], "")
        self.assertEqual(spec["baseline_model_sha256"],
                         "ce9e1112ee940531469457009234a2fd11de7f589c5948a7dec6b91118005294")
        # Notes are written to SBML; these characters are escaped or dropped there.
        text = json.dumps([spec[k]["notes"] for k in ("closed_reaction", "new_reaction", "biomass_term")])
        self.assertFalse(set("<>&") & set(text))

    def test_stale_inputs_fail_closed_without_edits(self):
        differs = "R2178 differs from the reviewed reaction"
        kinase = "R287 is not the reviewed dolichol kinase"
        biomass = "biomass_C composition differs"
        cases = {
            "r2178_name": (lambda m: setattr(m.reactions.R2178, "name", "dehydrodolichol reductase"), differs),
            "r2178_gpr": (lambda m: setattr(m.reactions.R2178, "gene_reaction_rule", "YALI1D26639g"), differs),
            "r2178_bounds": (lambda m: setattr(m.reactions.R2178, "bounds", (0.0, 500.0)), differs),
            "r2178_annotation": (lambda m: m.reactions.R2178.annotation.update({"ec-code": "1.3.1.94"}), differs),
            "r2178_stoichiometry": (lambda m: m.reactions.R2178.add_metabolites({m.metabolites.get_by_id("m10[C_cy]"): 1}), differs),
            "r2178_group": (lambda m: next(g for g in m.groups if m.reactions.R2178 in g.members).remove_members([m.reactions.R2178]), differs),
            "r287_gpr": (lambda m: setattr(m.reactions.R287, "gene_reaction_rule", "YALI1E22718g or YALI1D26639g"), kinase),
            "r287_bounds": (lambda m: setattr(m.reactions.R287, "bounds", (0.0, 0.0)), kinase),
            "r287_stoichiometry": (lambda m: m.reactions.R287.add_metabolites({m.metabolites.get_by_id("m32[C_cy]"): 1}), kinase),
            "r287_name": (lambda m: setattr(m.reactions.R287, "name", "polyprenol kinase"), kinase),
            "m2019_formula": (lambda m: setattr(m.metabolites.get_by_id("m2019[C_cy]"), "formula", "C80H132O7P2"), "m2019\\[C_cy\\] identity differs"),
            "m499_formula": (lambda m: setattr(m.metabolites.get_by_id("m499[C_cy]"), "formula", "C80H132O"), "m499\\[C_cy\\] identity differs"),
            "m503_charge": (lambda m: setattr(m.metabolites.get_by_id(TERM), "charge", -2), "m503\\[C_er\\] identity differs"),
            "nadph_compartment": (lambda m: setattr(m.metabolites.get_by_id("m40[C_cy]"), "compartment", "C_er"), "m40\\[C_cy\\] identity differs"),
            "dolichol_new_producer": (self._add_dolichol_producer, "m499\\[C_cy\\] reaction context differs"),
            "biomass_coefficient": (lambda m: m.reactions.biomass_C.add_metabolites({m.metabolites.get_by_id("m69[C_cy]"): -1e-5}), biomass),
            "biomass_bounds": (lambda m: setattr(m.reactions.biomass_C, "bounds", (0.0, 10.0)), biomass),
            "biomass_extra_note": (lambda m: m.reactions.biomass_C.notes.update({"other_curation": "x"}), "notes outside this curation"),
            "objective": (lambda m: setattr(m, "objective", "R287"), "objective is not biomass_C"),
            # A missing group trips the R287 guard first (R287 is a member); a duplicated name reaches this one.
            "group_not_unique": (lambda m: m.add_groups([Group("G_TEST_DUP", name="N-Glycan biosynthesis")]), "Group 'N-Glycan biosynthesis' is missing or not unique"),
            "r2178_missing": (lambda m: m.remove_reactions([m.reactions.R2178]), "R2178 missing"),
            "r287_missing": (lambda m: m.remove_reactions([m.reactions.R287]), "R287 missing"),
            "r2178_conflicting_notes": (lambda m: m.reactions.R2178.notes.update({"dolichol_route_curation_id": "other"}), "Conflicting R2178"),
            "biomass_conflicting_notes": (lambda m: m.reactions.biomass_C.notes.update({"dolichol_route_curation_id": "other"}), "Conflicting biomass_C"),
            "r2178_closed_elsewhere": (lambda m: setattr(m.reactions.R2178, "bounds", (0.0, 0.0)), "already closed without this curation"),
            "r2178_open_with_notes": (lambda m: m.reactions.R2178.notes.update(self.spec["closed_reaction"]["notes"]), "open but carries this curation's notes"),
            "term_elsewhere": (lambda m: m.reactions.biomass_C.add_metabolites({m.metabolites.get_by_id(TERM): -1e-4}), "already carries m503\\[C_er\\] without"),
            "term_other_coefficient": (lambda m: m.reactions.biomass_C.add_metabolites({m.metabolites.get_by_id(TERM): -5e-5}), "Conflicting m503\\[C_er\\] biomass coefficient"),
            "notes_without_term": (lambda m: m.reactions.biomass_C.notes.update(self.spec["biomass_term"]["notes"]), "notes without the m503"),
            "new_id_taken": (self._add_other_new_reaction, "R_DOLSYN already exists with another definition"),
            "biomass_missing": (self._remove_biomass, "biomass_C missing"),
            "species_missing": (lambda m: m.remove_metabolites([m.metabolites.get_by_id("m123[C_cy]")]), "m123\\[C_cy\\] missing"),
            "dehydrodolichol_new_consumer": (self._add_consumer("m2019[C_cy]"), "m2019\\[C_cy\\] reaction context differs"),
            "dolp_new_consumer": (self._add_consumer(TERM), "m503\\[C_er\\] reaction context differs"),
        }
        for field, (mutate, message) in cases.items():
            bad = self.model.copy()
            mutate(bad)
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_dolichol_route(bad)
            self.assertEqual(state(bad), before, field)

    @staticmethod
    def _add_dolichol_producer(model):
        reaction = Reaction("R_TEST_DOL", lower_bound=0.0, upper_bound=1000.0)
        model.add_reactions([reaction])
        reaction.add_metabolites({model.metabolites.get_by_id("m499[C_cy]"): 1})

    @staticmethod
    def _remove_biomass(model):
        model.remove_reactions([model.reactions.biomass_C])
        # Flush the lazy solver update, so the snapshot taken before the call is the state the step sees.
        model.solver.update()
        model.objective.expression

    @staticmethod
    def _add_consumer(mid):
        def mutate(model):
            reaction = Reaction("R_TEST_SINK", lower_bound=0.0, upper_bound=1000.0)
            model.add_reactions([reaction])
            reaction.add_metabolites({model.metabolites.get_by_id(mid): -1})
        return mutate

    @staticmethod
    def _add_other_new_reaction(model):
        reaction = Reaction(NEW, lower_bound=0.0, upper_bound=1000.0)
        model.add_reactions([reaction])
        reaction.add_metabolites({model.metabolites.get_by_id("m2019[C_cy]"): -1, model.metabolites.get_by_id("m499[C_cy]"): 1})

    def test_partial_states_are_rejected(self):
        full = expected_model(self.spec)
        for label in ("without_closure", "without_term", "without_reaction"):
            model = full.copy()
            if label == "without_closure":
                model.reactions.R2178.bounds = (0.0, 1000.0)
                for key in self.spec["closed_reaction"]["notes"]:
                    del model.reactions.R2178.notes[key]
            elif label == "without_term":
                model.reactions.biomass_C.add_metabolites({model.metabolites.get_by_id(TERM): 1e-4})
                for key in self.spec["biomass_term"]["notes"]:
                    del model.reactions.biomass_C.notes[key]
            else:
                model.remove_reactions([model.reactions.get_by_id(NEW)])
            before = state(model)
            with self.subTest(label=label), self.assertRaisesRegex(ValueError, "Only part of the dolichol route"):
                apply_dolichol_route(model)
            self.assertEqual(state(model), before, label)

    def test_unexpected_scope_is_rejected(self):
        for path, value in ((("schema_version",), 2), (("curation_id",), "DOLICHOL-ROUTE-20261029"),
                            (("status",), "proposed"), (("requires_objective",), "xBIOMASS"),
                            (("closed_reaction", "reaction_id"), "R2177"), (("closed_reaction", "after_bounds"), [0.0, 1000.0]),
                            (("kept_reaction", "reaction_id"), "R288"), (("new_reaction", "reaction_id"), "R_DOLK2"),
                            (("new_reaction", "gpr"), "YALI1E09357g"), (("new_reaction", "bounds"), [-1000.0, 1000.0]),
                            (("biomass_term", "metabolite_id"), "m501[C_cy]"),
                            (("biomass_term", "coefficient_mmol_per_gDW"), 0),
                            (("biomass_term", "accepted_before_state", "id"), "a1237821_biomass_C"),
                            (("new_reaction", "notes"), {}),
                            # A species without an identity row would be looked up only after R2178 is closed.
                            (("new_reaction", "stoichiometry", "m1[C_cy]"), 1.0)):
            spec = copy.deepcopy(self.spec)
            target = spec
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaisesRegex(ValueError, "Unexpected dolichol route curation scope"):
                apply_dolichol_route(self.model.copy(), spec)

    def test_makes_exactly_the_three_edits_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        expected = expected_model(self.spec)
        self.assertEqual(apply_dolichol_route(model)["status"], "applied")
        self.assertTrue(state(model) == state(expected), "Unexpected model change")
        self.assertEqual(model.reactions.R2178.bounds, (0.0, 0.0))
        self.assertEqual(model.reactions.biomass_C.metabolites[model.metabolites.get_by_id(TERM)], -1e-4)
        self.assertEqual(model.reactions.get_by_id(NEW).gene_reaction_rule, "")
        self.assertEqual(apply_dolichol_route(model)["status"], "already_correct")
        self.assertTrue(state(model) == state(expected), "Rerun changed the model")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dolichol.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(expected), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_dolichol_route(loaded)["status"], "already_correct")

    def test_flag_is_opt_in_and_requires_the_r1726_opening(self):
        self.assertFalse(parse_args([]).dolichol_route)
        self.assertTrue(parse_args(["--dolichol-route"]).dolichol_route)
        with self.assertRaisesRegex(ValueError, "requires --r1726-cpt-open"):
            build_reference_chain(dolichol_route=True)
        with self.assertRaisesRegex(ValueError, "requires --r1726-cpt-open"):
            build_reference_chain(dolichol_route=True, lipoylation_module=True, cofactor_biomass_terms=True,
                                  r2135_acyl_coa_pool=True)


if __name__ == "__main__":
    unittest.main()
