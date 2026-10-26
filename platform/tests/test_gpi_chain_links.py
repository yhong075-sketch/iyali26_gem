"""GPI chain links: preservation, rejection, SBML roundtrip, balance, dead-end closure and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.gpi_chain_links import CURATION_PATH, REACTION_IDS, apply_gpi_chain_links
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = "E5_coq9_resp_pkg_r132_combined_20261020.xml"
DEAD_ENDS = ("m596[C_er]", "m597[C_er]", "m599[C_er]", "m603[C_er]", "m604[C_er]")


def group_table(model):
    return {g.id: (g.name, sorted(r.id for r in g.members)) for g in model.groups}


def state(model):
    return semantics(model), annotations(model), group_table(model)


class GpiChainLinksTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        expected = {"dead_end_consumer": "m604\\[C_er\\] reaction context differs",
                    "intermediate_consumer": "m602\\[C_er\\] reaction context differs",
                    "species_charge": "m403\\[C_cy\\] identity differs",
                    "species_compartment": "m1648\\[C_em\\] identity differs",
                    "chain_bounds": "R376 differs from the recorded GPI chain",
                    "chain_stoichiometry": "R371 differs from the recorded GPI chain",
                    "group_missing": "Group 'Transport, ER to ER membrane' is missing",
                    "other_definition": "R_DM_GPI already exists with another definition",
                    "changed_note": "R_GPIGNT already exists with another definition",
                    "partial": "Only some GPI chain links are present"}
        baseline = state(self.model)
        for field, message in expected.items():
            bad = self.model.copy()
            if field == "dead_end_consumer":
                # E.g. the later R373 co-product change, applied without updating this curation.
                sink = Reaction("R_TEST_SINK", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([sink])
                sink.add_metabolites({bad.metabolites.get_by_id("m604[C_er]"): -1})
            elif field == "intermediate_consumer":
                # E.g. a later GPI7 step on Man3; the chain would no longer be linear.
                sink = Reaction("R_TEST_SINK", lower_bound=0.0, upper_bound=1000.0)
                bad.add_reactions([sink])
                sink.add_metabolites({bad.metabolites.get_by_id("m602[C_er]"): -1})
            elif field == "species_charge":
                bad.metabolites.get_by_id("m403[C_cy]").charge = -1
            elif field == "species_compartment":
                bad.metabolites.get_by_id("m1648[C_em]").compartment = "C_er"
            elif field == "chain_bounds":
                bad.reactions.R376.bounds = (0.0, 1000.0)
            elif field == "chain_stoichiometry":
                # The extra H+ is a separate candidate; this curation must notice if R371 changes.
                bad.reactions.R371.add_metabolites({bad.metabolites.get_by_id("m102[C_er]"): -1})
            elif field == "group_missing":
                next(g for g in bad.groups if g.name == "Transport, ER to ER membrane").name = "renamed"
            elif field == "other_definition":
                demand = Reaction("R_DM_GPI", lower_bound=1e-5, upper_bound=1000.0)
                bad.add_reactions([demand])
                demand.add_metabolites({bad.metabolites.get_by_id("m603[C_er]"): -1})
            elif field == "changed_note":
                apply_gpi_chain_links(bad)
                bad.reactions.R_GPIGNT.notes = {**bad.reactions.R_GPIGNT.notes, "gpi_chain_links_status": "other"}
            else:
                apply_gpi_chain_links(bad)
                bad.remove_reactions([bad.reactions.R_ACt_er_cy])
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_gpi_chain_links(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for key, value in (("status", "user_authorized"), ("curation_id", "GPI-CHAIN-LINKS"),
                           ("reactions", {k: v for k, v in self.spec["reactions"].items() if k != "R_DM_GPI"})):
            spec = copy.deepcopy(self.spec)
            spec[key] = value
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_gpi_chain_links(self.model.copy(), spec)

    def test_adds_only_the_five_links_and_survives_roundtrip(self):
        # Fresh reads, not copies: Model.copy() formats solver bounds differently (0 vs 0.0).
        model = read_sbml_model(str(MODEL.candidate_file(BASELINE)))
        before_species, before_annotations, before_groups = state(model)
        result = apply_gpi_chain_links(model)
        self.assertEqual((result["status"], set(result["added_reactions"])), ("applied", REACTION_IDS))
        after_species, after_annotations, after_groups = state(model)
        self.assertEqual(after_species["species"], before_species["species"])
        self.assertEqual(after_species["genes"], before_species["genes"])
        self.assertEqual(after_species["objective"], before_species["objective"])
        old = before_species["reactions"]
        new = after_species["reactions"]
        self.assertEqual(set(new) - set(old), REACTION_IDS)
        self.assertEqual({rid: new[rid] for rid in old}, old)
        self.assertEqual({rid: after_annotations[rid] for rid in before_annotations}, before_annotations)
        added = {gid: sorted(set(after_groups[gid][1]) - set(before_groups[gid][1])) for gid in before_groups}
        self.assertEqual({gid: ids for gid, ids in added.items() if ids},
                         {"group37": ["R_DM_GPI", "R_GPIGNT"], "group79": ["R_ACt_er_cy"],
                          "group87": ["R_DAGt_er_em", "R_PEt_em_er"]})
        self.assertEqual({gid: after_groups[gid][1] for gid in before_groups if not added[gid]},
                         {gid: before_groups[gid][1] for gid in before_groups if not added[gid]})
        for rid, row in self.spec["reactions"].items():
            reaction = model.reactions.get_by_id(rid)
            with self.subTest(reaction=rid):
                self.assertEqual({m.id: c for m, c in reaction.metabolites.items()}, row["stoichiometry"])
                self.assertEqual((list(reaction.bounds), reaction.gene_reaction_rule), (row["bounds"], ""))
                if rid != "R_DM_GPI":
                    self.assertEqual(reaction.check_mass_balance(), {})
        self.assertEqual(apply_gpi_chain_links(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "gpi.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(model), "Unexpected SBML roundtrip change")
        self.assertEqual({rid: loaded.reactions.get_by_id(rid).subsystem for rid in REACTION_IDS},
                         {rid: row["group"] for rid, row in self.spec["reactions"].items()})
        self.assertEqual(apply_gpi_chain_links(loaded)["status"], "already_correct")

    def test_chain_dead_ends_gain_a_producer_and_a_consumer(self):
        model = self.model.copy()
        apply_gpi_chain_links(model)
        for mid in DEAD_ENDS:
            metabolite = model.metabolites.get_by_id(mid)
            makes = [r.id for r in metabolite.reactions
                     if (r.metabolites[metabolite] > 0 and r.upper_bound > 0)
                     or (r.metabolites[metabolite] < 0 and r.lower_bound < 0)]
            uses = [r.id for r in metabolite.reactions
                    if (r.metabolites[metabolite] < 0 and r.upper_bound > 0)
                    or (r.metabolites[metabolite] > 0 and r.lower_bound < 0)]
            with self.subTest(metabolite=mid):
                self.assertTrue(set(makes) & REACTION_IDS or set(uses) & REACTION_IDS)
                self.assertTrue(makes and uses)

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).gpi_chain_links)
        self.assertTrue(parse_args(["--gpi-chain-links"]).gpi_chain_links)


if __name__ == "__main__":
    unittest.main()
