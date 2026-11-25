"""Keep test runs from writing into the model's report folders."""

import dataclasses
import os
import tempfile
from pathlib import Path

import pytest

try:
    import cobra
except ImportError:  # The evidence CI job installs only pydantic, pyyaml and pytest.
    cobra = None

from scripts.gem_annotate.model_layout import SCRATCH_DIR


def pytest_configure(config):
    # test_energy_candidate_behavior charges a solver budget and writes its outputs into the
    # task's report folder unless told otherwise; send both to a fresh temporary folder.
    scratch = Path(tempfile.mkdtemp(prefix="iyali26_tests_"))
    # Tests that build candidate files use SCRATCH_DIR: builders refuse outputs outside the repository.
    SCRATCH_DIR.mkdir(exist_ok=True)
    os.environ.setdefault("IYALI26_ENERGY_BUDGET", str(scratch / "budget.json"))
    os.environ.setdefault("IYALI26_ENERGY_TEST_OUTPUT", str(scratch / "energy_behavior"))


@pytest.hookimpl(tryfirst=True)
def pytest_report_to_serializable(report):
    # Called by pytest-xdist only: a worker sends each report to the controller with execnet, which takes
    # plain values only. A sub-test label holding anything else (a Path, say) would fail its whole test
    # with a DumpError, so that value goes as its repr. Returning None leaves the serialising to pytest.
    context = getattr(report, "context", None)
    if context is not None:
        report.context = dataclasses.replace(
            context, kwargs={key: _sendable(value) for key, value in context.kwargs.items()})


def _sendable(value):
    import execnet  # pytest-xdist depends on it, so it is there whenever the hook above runs

    try:
        execnet.dumps(value)
    except execnet.DumpError:
        return repr(value)
    return value


@pytest.fixture(autouse=True, scope="module")
def restore_default_solver():
    # Some tools set cobra's process-wide default solver (tools/lipid/lipid_moiety_ledger.py sets
    # GLPK); restore it after each test file so later files read their models with the usual one.
    if cobra is None:
        yield
        return
    configuration = cobra.Configuration()
    solver = configuration.solver
    yield
    configuration.solver = solver
