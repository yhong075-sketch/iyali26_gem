"""BiGG crosswalk annotations: exact targets, fail-closed guards, idempotence, SBML roundtrip and opt-in flag; no LPs."""

import copy
import tempfile
import unittest
from pathlib import Path

from cobra.io import read_sbml_model

from scripts.gem_annotate.bigg_crosswalk_annotations import apply_bigg_crosswalk_annotations, load_decisions
from scripts.gem_annotate.cli import parse_args
from scripts.gem_annotate.model_layout import MODEL
from scripts.gem_annotate.sbml import write_deterministic_sbml_model
from tests.test_coq9_curation import semantics

BASELINE = ("E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_"
            "r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_20261101.xml")


def annotation_table(model):
    return {**{("metabolite", m.id): copy.deepcopy(m.annotation) for m in model.metabolites},
            **{("reaction", r.id): copy.deepcopy(r.annotation) for r in model.reactions},
            **{("gene", g.id): copy.deepcopy(g.annotation) for g in model.genes}}


def mismatches(got, expected):
    """Keys whose values differ; keeps a failure readable (a full-model dict diff takes minutes)."""
    return sorted(k for k in set(got) | set(expected) if got.get(k) != expected.get(k))


def rebind(model, row, value):
    collection = model.metabolites if row["entity"] == "metabolite" else model.reactions
    obj = collection.get_by_id(row["entity_id"])
    obj.annotation = {**obj.annotation, row["key"]: value}


class BiggCrosswalkAnnotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # A Path, not a str: cobra rejects str paths of 260+ characters on macOS.
        cls.model = read_sbml_model(MODEL.candidate_file(BASELINE))
        cls.manifest, cls.rows = load_decisions()

    def test_decision_set_matches_the_approved_scope(self):
        self.assertEqual(len(self.rows), 896)
        self.assertEqual(self.manifest["actions"], {"write_his": 895, "drop_wrong": 1})
        self.assertEqual(len({(r["entity"], r["entity_id"]) for r in self.rows}), self.manifest["entities"])
        self.assertTrue(all(r["after"] == r["andrew_value"] for r in self.rows if r["action"] == "write_his"))

    def test_changes_exactly_the_reviewed_values(self):
        before = annotation_table(self.model)
        model = self.model.copy()
        # Compare the copy with itself: Model.copy() prints solver bounds differently ("0.0" vs "0").
        mathematics = semantics(model)
        result = apply_bigg_crosswalk_annotations(model)
        self.assertEqual((result["status"], result["rows"], result["entities"]), ("applied", 896, 595))
        expected = copy.deepcopy(before)
        for row in self.rows:
            values = row["after"].split(";")
            expected[(row["entity"], row["entity_id"])][row["key"]] = values[0] if len(values) == 1 else values
        self.assertEqual(mismatches(annotation_table(model), expected), [])
        self.assertEqual(mismatches(semantics(model), mathematics), [])
        self.assertEqual(mismatches(annotation_table(self.model), before), [])

    def test_second_run_is_already_correct(self):
        model = self.model.copy()
        apply_bigg_crosswalk_annotations(model)
        applied = annotation_table(model)
        self.assertEqual(apply_bigg_crosswalk_annotations(model)["status"], "already_correct")
        self.assertEqual(mismatches(annotation_table(model), applied), [])

    def test_stale_or_partly_applied_targets_fail_closed_without_edits(self):
        stale = self.model.copy()
        rebind(stale, self.rows[0], "MNXM0")
        partial = self.model.copy()
        apply_bigg_crosswalk_annotations(partial)
        rebind(partial, self.rows[-1], self.rows[-1]["before"].split(";"))
        for model, culprit in ((stale, self.rows[0]), (partial, self.rows[-1])):
            snapshot = annotation_table(model)
            with self.assertRaisesRegex(ValueError, f"differ from the reviewed state .*'{culprit['decision_id']}'"):
                apply_bigg_crosswalk_annotations(model)
            self.assertEqual(mismatches(annotation_table(model), snapshot), [])

    def test_unexpected_scope_or_rows_fail_closed(self):
        # Every case raises before any edit, so the shared model is checked rather than copied.
        snapshot = annotation_table(self.model)
        for field, value in (("status", "draft"), ("decision_set_id", "OTHER"), ("rows", 1), ("schema_version", 2)):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "scope/status"):
                apply_bigg_crosswalk_annotations(self.model, {**self.manifest, field: value}, self.rows)
        first = self.rows[0]
        bad_rows = {"unexpected entity/key/action": [dict(first, key="inchikey"), dict(first, action="write_ours"),
                                                     dict(first, after=""), dict(first, after=first["before"])],
                    "missing": [dict(first, entity_id="m0[C_cy]")]}
        for message, variants in bad_rows.items():
            for bad in variants:
                with self.subTest(row=bad), self.assertRaisesRegex(ValueError, message):
                    apply_bigg_crosswalk_annotations(self.model, self.manifest, self.rows[1:] + [bad])
        with self.assertRaisesRegex(ValueError, "duplicate target"):
            apply_bigg_crosswalk_annotations(self.model, {**self.manifest, "rows": len(self.rows) + 1},
                                             self.rows + [first])
        self.assertEqual(mismatches(annotation_table(self.model), snapshot), [])

    def test_sbml_roundtrip_keeps_the_values(self):
        model = self.model.copy()
        apply_bigg_crosswalk_annotations(model)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "model.xml"
            write_deterministic_sbml_model(model, path)
            loaded = read_sbml_model(path)
        targets = {(r["entity"], r["entity_id"]) for r in self.rows}
        self.assertEqual(mismatches({k: v for k, v in annotation_table(loaded).items() if k in targets},
                                    {k: v for k, v in annotation_table(model).items() if k in targets}), [])
        self.assertEqual(apply_bigg_crosswalk_annotations(loaded)["status"], "already_correct")

    def test_flag_is_opt_in(self):
        self.assertFalse(parse_args([]).bigg_crosswalk_annotations)
        self.assertTrue(parse_args(["--bigg-crosswalk-annotations"]).bigg_crosswalk_annotations)


if __name__ == "__main__":
    unittest.main()
