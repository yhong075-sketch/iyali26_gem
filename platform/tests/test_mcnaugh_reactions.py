"""mcnaugh_pass reactions: exact additions and closures, fail-closed guards, idempotence, SBML roundtrip, opt-in flag."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

import cobra
from cobra.io import read_sbml_model

from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.mcnaugh_reactions import CURATION_PATH, apply_mcnaugh_reactions
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_"
            "r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_20261103.xml")


def network(model):
    """Reaction definitions and species, the parts this change may touch."""
    return {"reactions": {r.id: ({m.id: c for m, c in r.metabolites.items()}, r.bounds, r.gene_reaction_rule,
                                 r.objective_coefficient) for r in model.reactions},
            "species": {m.id: (m.name, m.formula, m.charge, m.compartment) for m in model.metabolites},
            "genes": sorted(g.id for g in model.genes)}


def metadata(model):
    """The reaction fields that are not part of the LP but are written to the file."""
    return {r.id: (r.name, r.subsystem, r.annotation, r.notes, sorted(g.name for g in model.groups if r in g.members))
            for r in model.reactions}


def mismatches(got, expected):
    return sorted(k for k in set(got) | set(expected) if got.get(k) != expected.get(k))


class McnaughReactionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # A Path, not a str: cobra rejects str paths of 260+ characters on macOS.
        cls.model = read_sbml_model(MODEL.candidate_file(BASELINE))
        cls.spec = json.loads(CURATION_PATH.read_text())

    def test_adds_and_closes_exactly_the_reviewed_set(self):
        model = self.model.copy()
        before, meta_before = network(model), metadata(model)
        result = apply_mcnaugh_reactions(model)
        self.assertEqual(result["status"], "applied")
        after, meta_after = network(model), metadata(model)
        added = {row["id"]: row for row in self.spec["reactions"]}
        closed = {row["id"]: row for row in self.spec["closures"]}
        self.assertEqual(mismatches(after["reactions"], before["reactions"]), sorted(set(added) | set(closed)))
        self.assertEqual(mismatches(meta_after, meta_before), sorted(added))
        for rid, row in added.items():
            self.assertEqual(after["reactions"][rid], (row["stoichiometry"], tuple(row["bounds"]), row["gpr"], 0.0))
            self.assertEqual(meta_after[rid], (row["name"], row["groups"][0] if row["groups"] else "", row["annotation"],
                                               row["notes"], row["groups"]))
            self.assertEqual(model.reactions.get_by_id(rid).notes["curation"], "MCNAUGH-REACTIONS-20261103")
        for rid, row in closed.items():
            self.assertEqual(after["reactions"][rid][1], (0.0, 0.0))
            self.assertEqual(after["reactions"][rid][0], before["reactions"][rid][0])
        self.assertEqual(after["species"], before["species"])
        self.assertEqual(after["genes"], before["genes"])
        self.assertEqual([r.id for r in model.reactions if r.objective_coefficient], ["biomass_C"])

    def test_lipid_scope_and_documented_exceptions(self):
        ids = [row["id"] for row in self.spec["reactions"]]
        self.assertEqual(len([i for i in ids if i.endswith("_lp")]), 15)
        self.assertEqual(sorted(i for i in ids if i.startswith("BIOMASS")), ["BIOMASS_C_core", "BIOMASS_C_full"])
        self.assertEqual(len(self.spec["not_added"]), 12)
        charged = sorted(row["id"] for row in self.spec["reactions"] if "charge_exception" in row["notes"])
        self.assertEqual(charged, ["AGPAT182_lp", "DGAT182_lp", "G3PAT182_lp"])
        # SBML escapes these characters, so notes containing them would not compare equal after a roundtrip.
        self.assertEqual([v for row in self.spec["reactions"] for v in row["notes"].values() if set(v) & set("<>&'\"")], [])
        self.assertEqual({row["id"]: row["groups"] for row in self.spec["reactions"]},
                         {i: [] if i.startswith("BIOMASS") else ["Glycerolipid metabolism"] for i in ids})

    def test_second_run_is_already_correct(self):
        model = self.model.copy()
        apply_mcnaugh_reactions(model)
        applied = network(model)
        self.assertEqual(apply_mcnaugh_reactions(model)["status"], "already_correct")
        self.assertEqual(mismatches(network(model)["reactions"], applied["reactions"]), [])

    def test_stale_or_partial_state_fails_closed_without_edits(self):
        cases = {}
        stale = self.model.copy()
        stale.reactions.R1771.bounds = (0.0, 500.0)
        cases["stale closure"] = (stale, "partly applied or stale")
        partial = self.model.copy()
        partial.reactions.R353.bounds = (0.0, 0.0)
        cases["closure already applied, reactions not"] = (partial, "partly applied or stale")
        rewritten = self.model.copy()
        rewritten.reactions.R1846.gene_reaction_rule = "YALI1E22736g"
        cases["closure GPR differs"] = (rewritten, "stoichiometry or GPR differs")
        restoich = self.model.copy()
        restoich.reactions.R1771.add_metabolites({restoich.metabolites.get_by_id("m573[C_lp]"): 1.0})
        cases["closure stoichiometry differs"] = (restoich, "stoichiometry or GPR differs")
        objective = self.model.copy()
        objective.objective = "R353"
        cases["other objective"] = (objective, "Objective is not biomass_C")
        no_group = self.model.copy()
        next(g for g in no_group.groups if g.name == "Glycerolipid metabolism").name = "renamed"
        cases["group missing"] = (no_group, "Group 'Glycerolipid metabolism' is missing")
        duplicate = self.model.copy()
        row = self.spec["reactions"][0]
        other = cobra.Reaction("OTHER_" + row["id"])
        duplicate.add_reactions([other])
        other.add_metabolites({duplicate.metabolites.get_by_id(m): c for m, c in row["stoichiometry"].items()})
        cases["same stoichiometry under another ID"] = (duplicate, "OTHER_G3PAT160_lp has the same stoichiometry")
        backward = self.model.copy()
        reverse = cobra.Reaction("REV_" + row["id"], lower_bound=-1000.0, upper_bound=0.0)
        backward.add_reactions([reverse])
        reverse.add_metabolites({backward.metabolites.get_by_id(m): -c for m, c in row["stoichiometry"].items()})
        cases["reverse stoichiometry under another ID"] = (backward, "REV_G3PAT160_lp has the same stoichiometry")
        # After the step, a duplicate earlier in reaction order than the new reaction must still be found.
        earlier = self.model.copy()
        apply_mcnaugh_reactions(earlier)
        first = earlier.reactions[0]
        self.assertNotIn(first.id, {r["id"] for r in self.spec["reactions"] + self.spec["closures"]})
        first.subtract_metabolites(first.metabolites)
        first.add_metabolites({earlier.metabolites.get_by_id(m): c for m, c in row["stoichiometry"].items()})
        cases["applied model with an earlier duplicate"] = (earlier, f"{first.id} has the same stoichiometry")
        for name, (model, message) in cases.items():
            snapshot = network(model)
            with self.subTest(case=name), self.assertRaisesRegex(ValueError, message):
                apply_mcnaugh_reactions(model)
            self.assertEqual(mismatches(network(model)["reactions"], snapshot["reactions"]), [])

    def test_rewritten_added_reaction_fails_closed_without_edits(self):
        applied = self.model.copy()
        apply_mcnaugh_reactions(applied)
        rewrites = {
            "notes": lambda m, r: setattr(r, "notes", {**r.notes, "mass_balance_exception": "other"}),
            "name": lambda m, r: setattr(r, "name", "other"),
            "annotation": lambda m, r: setattr(r, "annotation", {**r.annotation, "sbo": "SBO:0000185"}),
            "group": lambda m, r: next(g for g in m.groups if g.name == "Glycerolipid metabolism").remove_members([r]),
        }
        for field, rewrite in rewrites.items():
            model = applied.copy()
            rewrite(model, model.reactions.G3PAT160_lp)
            snapshot, meta = network(model), metadata(model)
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "G3PAT160_lp exists with another definition"):
                apply_mcnaugh_reactions(model)
            self.assertEqual(mismatches(network(model)["reactions"], snapshot["reactions"]), [])
            self.assertEqual(mismatches(metadata(model), meta), [])

    def test_unexpected_scope_or_rows_fail_closed(self):
        for field, value in (("status", "draft"), ("change_id", "OTHER"), ("schema_version", 2)):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_mcnaugh_reactions(self.model, {**self.spec, field: value})
        spec = copy.deepcopy(self.spec)
        spec["reactions"][0]["stoichiometry"]["m0[C_lp]"] = 1.0
        with self.assertRaisesRegex(ValueError, "species or genes missing"):
            apply_mcnaugh_reactions(self.model, spec)
        spec = copy.deepcopy(self.spec)
        spec["reactions"][0]["gpr"] = "YALI1Z99999g"
        with self.assertRaisesRegex(ValueError, "species or genes missing"):
            apply_mcnaugh_reactions(self.model, spec)
        spec = copy.deepcopy(self.spec)
        spec["reactions"][0]["id"] = "R353"
        with self.assertRaisesRegex(ValueError, "exists with another definition"):
            apply_mcnaugh_reactions(self.model, spec)
        self.assertNotIn("G3PAT160_lp", self.model.reactions)

    def test_sbml_roundtrip_keeps_the_change(self):
        model = self.model.copy()
        apply_mcnaugh_reactions(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(path)
        self.assertEqual(mismatches(network(loaded)["reactions"], network(model)["reactions"]), [])
        touched = {row["id"] for row in self.spec["reactions"] + self.spec["closures"]}
        meta_loaded, meta_model = metadata(loaded), metadata(model)
        self.assertEqual(mismatches({k: meta_loaded[k] for k in touched}, {k: meta_model[k] for k in touched}), [])
        self.assertEqual(loaded.reactions.G3PAT160_lp.subsystem, "Glycerolipid metabolism")
        self.assertEqual(apply_mcnaugh_reactions(loaded)["status"], "already_correct")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).mcnaugh_reactions)
        self.assertTrue(parse_args(["--mcnaugh-reactions"]).mcnaugh_reactions)


if __name__ == "__main__":
    unittest.main()
