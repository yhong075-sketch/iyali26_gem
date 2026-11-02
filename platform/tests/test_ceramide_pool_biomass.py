"""Ceramide pool biomass term: preservation, rejection, SBML roundtrip, context exit and opt-in flag; no LPs."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from cobra import Reaction
from cobra.io import read_sbml_model

from scripts.gem_annotate.ceramide_pool_biomass import (BIOMASS, CERAMIDES, CURATION_PATH, POOL, REACTION_PREFIX,
                                                         apply_ceramide_pool_biomass)
from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import annotations, semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_"
            "gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_20261031.xml")
REACTIONS = [f"{REACTION_PREFIX}{c}" for c in CERAMIDES]


def state(model):
    return semantics(model), annotations(model)


def read_baseline():
    # Path object: cobra rejects str paths of 260+ characters on macOS ("win" in "darwin").
    return read_sbml_model(MODEL.candidate_file(BASELINE))


class CeramidePoolBiomassTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = read_baseline()
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_stale_inputs_fail_closed_without_edits(self):
        composition = "composition differs from the reviewed baseline"
        expected = {"other_coefficient": composition, "extra_entry": composition, "bounds": composition,
                    "identity": "identity differs", "missing_ceramide": "missing",
                    "notes": "Conflicting", "extra_note": "notes outside this curation",
                    "stray_reaction": "Partial ceramide pool state", "notes_only": "Partial ceramide pool state",
                    "objective": "objective is not", "direction": "objective is not", "missing": "missing"}
        baseline = state(self.model)
        notes = self.spec["biomass_notes"]
        for field, message in expected.items():
            bad = self.model.copy()
            biomass = bad.reactions.get_by_id(BIOMASS)
            if field == "other_coefficient":
                biomass.add_metabolites({bad.metabolites.get_by_id("m1324[C_cy]"): -0.1})
            elif field == "extra_entry":
                biomass.add_metabolites({bad.metabolites.get_by_id("m386[C_er]"): -0.001})
            elif field == "bounds":
                biomass.bounds = (0.0, 100.0)
            elif field == "identity":
                bad.metabolites.get_by_id("m397[C_er]").formula = "C42H85NO4"
            elif field == "missing_ceramide":
                bad.remove_metabolites([bad.metabolites.get_by_id("m400[C_er]")])
            elif field == "notes":
                biomass.notes = {**biomass.notes, "ceramide_pool_biomass_status": "another review"}
            elif field == "extra_note":
                biomass.notes = {**biomass.notes, "unreviewed_note": "added elsewhere"}
            elif field == "stray_reaction":
                stray = Reaction(REACTIONS[0], lower_bound=0, upper_bound=1000)
                stray.add_metabolites({bad.metabolites.get_by_id("m386[C_er]"): -1})
                bad.add_reactions([stray])
            elif field == "notes_only":
                biomass.notes = {**biomass.notes, **copy.deepcopy(notes)}
            elif field == "objective":
                bad.objective = "xBIOMASS"
            elif field == "direction":
                bad.objective.direction = "min"
            else:
                bad.remove_reactions([biomass])
            bad.solver.update()  # optlang refreshes the objective lazily; flush before the snapshot
            before = state(bad)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, message):
                apply_ceramide_pool_biomass(bad)
            self.assertEqual(state(bad), before)
        self.assertEqual(state(self.model), baseline, "A subcase leaked into the shared model")

    def test_unexpected_scope_is_rejected(self):
        for edit in ("status", "curation_id", "coefficient", "pool", "feeds", "bounds", "notes", "before_with_pool"):
            spec = copy.deepcopy(self.spec)
            if edit == "status":
                spec["status"] = "user_authorized"
            elif edit == "curation_id":
                spec["curation_id"] = "CERAMIDE-POOL"
            elif edit == "coefficient":
                spec["pool"]["coefficient"] = -5e-4
            elif edit == "pool":
                spec["pool"]["metabolite_id"] = "m386[C_er]"
            elif edit == "feeds":
                spec["pool"]["feeds"] = spec["pool"]["feeds"][:-1]
            elif edit == "bounds":
                spec["pool"]["feeds"][0]["bounds"] = [-1000.0, 1000.0]
            elif edit == "notes":
                spec["biomass_notes"] = {}
            else:
                spec["accepted_before_state"]["biomass_stoichiometry"][POOL] = -5.4e-4
            with self.subTest(edit=edit), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_ceramide_pool_biomass(self.model.copy(), spec)

    def test_adds_only_the_pool_and_survives_roundtrip(self):
        model, baseline = read_baseline(), read_baseline()
        result = apply_ceramide_pool_biomass(model)
        self.assertEqual(result["status"], "applied")
        self.assertEqual(result["reactions_added"], REACTIONS)
        biomass = model.reactions.get_by_id(BIOMASS)
        self.assertEqual(biomass.get_coefficient(POOL), -5.4e-4)
        self.assertEqual(sorted(r.id for r in model.metabolites.get_by_id(POOL).reactions), sorted([*REACTIONS, BIOMASS]))
        for rid, ceramide in zip(REACTIONS, CERAMIDES):
            reaction = model.reactions.get_by_id(rid)
            self.assertEqual((reaction.bounds, reaction.gene_reaction_rule), ((0.0, 1000.0), ""))
            self.assertEqual({m.id: c for m, c in reaction.metabolites.items()}, {f"{ceramide}[C_er]": -1.0, POOL: 1.0})
        # Everything outside the targets is unchanged.
        old, new = baseline.reactions.get_by_id(BIOMASS), biomass
        self.assertEqual({m.id: c for m, c in old.metabolites.items()},
                         {m.id: c for m, c in new.metabolites.items() if m.id != POOL})
        self.assertEqual(set(new.notes) - set(old.notes), set(self.spec["biomass_notes"]))
        untouched = [r.id for r in baseline.reactions if r.id != BIOMASS]
        self.assertEqual({r: (baseline.reactions.get_by_id(r).bounds, baseline.reactions.get_by_id(r).gene_reaction_rule,
                              {m.id: c for m, c in baseline.reactions.get_by_id(r).metabolites.items()}) for r in untouched},
                         {r: (model.reactions.get_by_id(r).bounds, model.reactions.get_by_id(r).gene_reaction_rule,
                              {m.id: c for m, c in model.reactions.get_by_id(r).metabolites.items()}) for r in untouched})
        self.assertEqual(len(model.reactions), len(baseline.reactions) + 10)
        self.assertEqual(len(model.metabolites), len(baseline.metabolites) + 1)
        self.assertEqual(len(model.genes), len(baseline.genes))
        self.assertEqual(apply_ceramide_pool_biomass(model)["status"], "already_correct")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ceramide_pool.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(str(path))
        self.assertTrue(state(loaded) == state(model), "Unexpected SBML roundtrip change")
        self.assertEqual(apply_ceramide_pool_biomass(loaded)["status"], "already_correct")
        # On an applied model the guards still hold.
        bad = loaded.copy()
        bad.reactions.get_by_id(REACTIONS[3]).bounds = (-1000.0, 1000.0)
        before = state(bad)
        with self.assertRaisesRegex(ValueError, "differs from the curated state"):
            apply_ceramide_pool_biomass(bad)
        self.assertEqual(state(bad), before)

    def test_context_exit_restores_the_baseline(self):
        model = read_baseline()
        before = state(model)
        with model:
            self.assertEqual(apply_ceramide_pool_biomass(model)["status"], "applied")
            self.assertIn(POOL, model.metabolites)
        self.assertTrue(state(model) == before, "Context exit left curation edits behind")
        self.assertNotIn(POOL, model.metabolites)
        self.assertEqual(apply_ceramide_pool_biomass(model)["status"], "applied")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).ceramide_pool_biomass)
        self.assertTrue(parse_args(["--ceramide-pool-biomass"]).ceramide_pool_biomass)


if __name__ == "__main__":
    unittest.main()
