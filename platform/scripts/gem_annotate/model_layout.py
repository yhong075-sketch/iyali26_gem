"""Locate a model folder's files from its model.toml, without other side effects."""

import os
import tomllib
from dataclasses import dataclass
from pathlib import Path

PLATFORM_ROOT = Path(__file__).resolve().parent.parent.parent
REPO_ROOT = PLATFORM_ROOT.parent
MODEL_DIR_ENV = "IYALI26_MODEL_DIR"
# Git-ignored folder for temporary outputs that must stay inside the repository.
SCRATCH_DIR = REPO_ROOT / ".scratch"


@dataclass(frozen=True, slots=True)
class ModelLayout:
    """Paths of one model folder, as declared in its model.toml."""

    root: Path
    start_model: Path
    canonical_model: Path
    curation: Path
    conditions: Path
    candidates: Path
    reports: Path

    def curation_file(self, name: str) -> Path:
        """Return the curation file with this unique name from any topic subfolder."""
        matches = [path for path in self.curation.rglob("*") if path.name == name and path.is_file()]
        if len(matches) > 1:
            raise FileExistsError(f"curation file name {name!r} is not unique: {matches}")
        return matches[0] if matches else self.curation / name

    def candidate_file(self, name: str) -> Path:
        """Return the candidate model file with this unique name from any subfolder."""
        matches = [path for path in self.candidates.rglob("*") if path.name == name and path.is_file()]
        if len(matches) > 1:
            raise FileExistsError(f"candidate file name {name!r} is not unique: {matches}")
        return matches[0] if matches else self.candidates / name


def load_model_layout(model_dir: str | Path | None = None) -> ModelLayout:
    """Read model.toml from the model folder (argument, environment, then repo default)."""
    root = Path(model_dir or os.environ.get(MODEL_DIR_ENV) or REPO_ROOT / "model").expanduser().resolve()
    with (root / "model.toml").open("rb") as handle:
        layout = tomllib.load(handle)["layout"]
    return ModelLayout(root=root, **{key: root / value for key, value in layout.items()})


MODEL = load_model_layout()
