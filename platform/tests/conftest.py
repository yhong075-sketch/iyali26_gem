"""Keep test runs from writing into the model's report folders."""

import os
import tempfile
from pathlib import Path

from scripts.gem_annotate.model_layout import SCRATCH_DIR


def pytest_configure(config):
    # test_energy_candidate_behavior charges a solver budget and writes its outputs into the
    # task's report folder unless told otherwise; send both to a fresh temporary folder.
    scratch = Path(tempfile.mkdtemp(prefix="iyali26_tests_"))
    # Tests that build candidate files use SCRATCH_DIR: builders refuse outputs outside the repository.
    SCRATCH_DIR.mkdir(exist_ok=True)
    os.environ.setdefault("IYALI26_ENERGY_BUDGET", str(scratch / "budget.json"))
    os.environ.setdefault("IYALI26_ENERGY_TEST_OUTPUT", str(scratch / "energy_behavior"))
