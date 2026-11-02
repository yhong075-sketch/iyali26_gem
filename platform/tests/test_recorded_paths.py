"""Absolute paths recorded in the former checkout's artifacts/ folder resolve like the relative form."""

import tempfile
import unittest
from pathlib import Path

from scripts.gem_annotate.config import load_project_paths

TASK = "recorded_path_test_fixture"
FORMER = Path("/former/checkout/iyali26_gem/artifacts") / TASK


class RecordedPathTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.paths = load_project_paths(self.temp.name)
        self.saved = self.paths.task_outputs / TASK / "run.json"
        self.saved.parent.mkdir(parents=True)
        self.saved.write_text("{}")

    def tearDown(self):
        self.temp.cleanup()

    def test_missing_former_artifacts_path_resolves_to_the_saved_file(self):
        self.assertEqual(self.paths.resolve_legacy_path(FORMER / "run.json"), self.saved.resolve())

    def test_missing_former_path_without_a_saved_file_is_left_as_given(self):
        for path in (FORMER / "new_output.json", Path("/former/checkout/other/artifacts") / TASK / "run.json"):
            self.assertEqual(self.paths.resolve_legacy_path(path), path)

    def test_existing_absolute_path_is_unchanged(self):
        self.assertEqual(self.paths.resolve_legacy_path(self.saved), self.saved.resolve())


if __name__ == "__main__":
    unittest.main()
