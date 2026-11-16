"""Reaction xref cleanup: exact removals, notes, fail-closed guards, idempotence, context, SBML roundtrip, flag; no LPs."""

import copy
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.bigg_crosswalk_annotations import load_decisions as load_crosswalk_decisions
from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.reaction_xref_cleanup import (
    KEYS, NOTES_KEY, apply_reaction_xref_cleanup, load_decisions)
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_"
            "r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_r2121_close_"
            "r1172_close_20261112.xml")


def table(model):
    """Annotation and notes of every reaction, metabolite and gene."""
    return {**{("reaction", r.id): (copy.deepcopy(r.annotation), dict(r.notes)) for r in model.reactions},
            **{("metabolite", m.id): (copy.deepcopy(m.annotation), dict(m.notes)) for m in model.metabolites},
            **{("gene", g.id): (copy.deepcopy(g.annotation), dict(g.notes)) for g in model.genes}}


def mismatches(got, expected):
    """Keys whose values differ; keeps a failure readable."""
    return sorted(k for k in set(got) | set(expected) if got.get(k) != expected.get(k))


def rebind(model, rid, key, value):
    reaction = model.reactions.get_by_id(rid)
    reaction.annotation = {**reaction.annotation, key: value}


class ReactionXrefCleanupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # A Path, not a str: cobra rejects str paths of 260+ characters on macOS.
        cls.model = read_sbml_model(MODEL.candidate_file(BASELINE))
        cls.manifest, cls.rows = load_decisions()
        _, cls.crosswalk_rows = load_crosswalk_decisions()

    def test_decision_set_matches_the_approved_scope(self):
        rows, manifest = self.rows, self.manifest
        self.assertEqual((len(rows), manifest["reactions"], len(manifest["notes"])), (374, 62, 62))
        self.assertEqual(sum(len(r["removed"].split(";")) for r in rows), 896)
        self.assertEqual({r["entity_id"] for r in rows}, set(manifest["notes"]))
        self.assertTrue(all(r["after"] == "" for r in rows))
        self.assertEqual(len(manifest["excluded_reactions"]), 51)
        self.assertFalse(set(manifest["excluded_reactions"]) & set(manifest["notes"]))
        self.assertEqual(manifest["base_model_sha256"],
                         "684f3752764fb8f25219008966953863141724f34a451dffb77c235aa3cb2f62")
        decided = {(r["entity_id"], r["key"]) for r in self.crosswalk_rows if r["entity"] == "reaction"}
        self.assertFalse({(r["entity_id"], r["key"]) for r in rows} & decided)
        # SBML escapes these in notes, so the read-back text would differ from the curated one.
        self.assertFalse(any(set(text) & set("<>&'\"") for text in manifest["notes"].values()))

    def test_removes_exactly_the_reviewed_values_and_adds_one_note(self):
        before = table(self.model)
        model = self.model.copy()
        # Compare the copy with itself: Model.copy() prints solver bounds differently ("0.0" vs "0").
        mathematics = semantics(model)
        result = apply_reaction_xref_cleanup(model)
        self.assertEqual((result["status"], result["rows"], result["reactions"], result["values_removed"]),
                         ("applied", 374, 62, 896))
        expected = copy.deepcopy(before)
        for row in self.rows:
            annotation, _ = expected[("reaction", row["entity_id"])]
            del annotation[row["key"]]
        for rid, text in self.manifest["notes"].items():
            expected[("reaction", rid)][1][NOTES_KEY] = text
        self.assertEqual(mismatches(table(model), expected), [])
        self.assertEqual(semantics(model), mathematics)
        self.assertEqual(mismatches(table(self.model), before), [])

    def test_second_run_is_already_correct(self):
        model = self.model.copy()
        apply_reaction_xref_cleanup(model)
        applied = table(model)
        self.assertEqual(apply_reaction_xref_cleanup(model)["status"], "already_correct")
        self.assertEqual(mismatches(table(model), applied), [])

    def test_stale_or_partly_applied_targets_fail_closed_without_edits(self):
        first, last = self.rows[0], self.rows[-1]
        stale = self.model.copy()
        rebind(stale, first["entity_id"], first["key"], "R00000")
        partial = self.model.copy()
        apply_reaction_xref_cleanup(partial)
        rebind(partial, last["entity_id"], last["key"], last["before"].split(";"))
        noted = self.model.copy()
        reaction = noted.reactions.get_by_id(first["entity_id"])
        reaction.notes = {**reaction.notes, NOTES_KEY: "other text"}
        for model, culprit in ((stale, first["decision_id"]), (partial, last["decision_id"]),
                               (noted, first["entity_id"])):
            snapshot = table(model)
            with self.subTest(culprit=culprit), \
                    self.assertRaisesRegex(ValueError, f"differ from the reviewed state .*'{culprit}'"):
                apply_reaction_xref_cleanup(model)
            self.assertEqual(mismatches(table(model), snapshot), [])

    def test_refuses_targets_whose_crosswalk_values_are_not_applied(self):
        # Without --bigg-crosswalk-annotations, R29 would keep the wrong MNXR while its note says the crosswalk MNXR
        # is kept.
        row = next(r for r in self.crosswalk_rows
                   if r["entity"] == "reaction" and r["entity_id"] == "R29" and r["key"] == "metanetx.reaction")
        model = self.model.copy()
        rebind(model, "R29", "metanetx.reaction", row["before"])
        snapshot = table(model)
        with self.assertRaisesRegex(ValueError, "BiGG crosswalk values on cleanup targets .*R29 metanetx.reaction"):
            apply_reaction_xref_cleanup(model)
        self.assertEqual(mismatches(table(model), snapshot), [])

    def test_unexpected_scope_or_rows_fail_closed(self):
        # Every case raises before any edit, so the shared model is checked rather than copied.
        snapshot = table(self.model)
        for field, value in (("status", "draft"), ("decision_set_id", "OTHER"), ("rows", 1), ("schema_version", 2),
                             ("keys", sorted(KEYS)[1:]), ("notes_key", "other"), ("notes", [])):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_reaction_xref_cleanup(self.model, {**self.manifest, field: value}, self.rows)
        first = self.rows[0]
        bad_rows = {
            "inconsistent values": [dict(first, key="ec-code"), dict(first, action="drop"), dict(first, removed=""),
                                    dict(first, removed="R99999"), dict(first, after="R99999")],
            "BiGG crosswalk decision": [dict(first, entity_id="R29", key="bigg.reaction")],
            # R1540 is a lipid-unlump exclusion; R1 is in neither the notes nor the exclusions.
            "excluded or has no note": [dict(first, entity_id="R1540"), dict(first, entity_id="R1")],
        }
        for message, variants in bad_rows.items():
            for bad in variants:
                with self.subTest(row=bad), self.assertRaisesRegex(ValueError, message):
                    apply_reaction_xref_cleanup(self.model, self.manifest, self.rows[1:] + [bad])
        with self.assertRaisesRegex(ValueError, "duplicate target"):
            apply_reaction_xref_cleanup(self.model, {**self.manifest, "rows": len(self.rows) + 1},
                                        self.rows + [first])
        missing = {**self.manifest, "notes": {**self.manifest["notes"], "R0": "note"}}
        with self.assertRaisesRegex(ValueError, "missing"):
            apply_reaction_xref_cleanup(self.model, missing, self.rows[1:] + [dict(first, entity_id="R0")])
        with self.assertRaisesRegex(ValueError, "name different reactions"):
            apply_reaction_xref_cleanup(self.model, missing, self.rows)
        self.assertEqual(mismatches(table(self.model), snapshot), [])

    def test_context_exit_restores_annotations_and_notes(self):
        model = self.model.copy()
        snapshot = table(model)
        with model:
            self.assertEqual(apply_reaction_xref_cleanup(model)["status"], "applied")
            self.assertNotEqual(mismatches(table(model), snapshot), [])
        self.assertEqual(mismatches(table(model), snapshot), [])

    def test_sbml_roundtrip_keeps_the_result(self):
        model = self.model.copy()
        apply_reaction_xref_cleanup(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(path)
        targets = {("reaction", rid) for rid in self.manifest["notes"]}
        self.assertEqual(mismatches({k: v for k, v in table(loaded).items() if k in targets},
                                    {k: v for k, v in table(model).items() if k in targets}), [])
        self.assertEqual(apply_reaction_xref_cleanup(loaded)["status"], "already_correct")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).reaction_xref_cleanup)
        self.assertTrue(parse_args(["--reaction-xref-cleanup"]).reaction_xref_cleanup)


if __name__ == "__main__":
    unittest.main()
