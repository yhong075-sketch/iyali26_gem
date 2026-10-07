#!/usr/bin/env python3
"""Print a read-only state banner at session start: branch, local changes, research root, latest log entry."""
import os
import subprocess
from pathlib import Path

root = Path(os.environ.get("CLAUDE_PROJECT_DIR", "."))


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True).stdout.strip()


branch = git("rev-parse", "--abbrev-ref", "HEAD")
head = git("rev-parse", "--short", "HEAD")
changes = len([line for line in git("status", "--porcelain").splitlines() if line])
research = os.environ.get("IYALI26_RESEARCH_ROOT")
state = root / "model" / "STATE.md"
latest = next((line[3:].strip() for line in state.read_text().splitlines() if line.startswith("## ")), "none") if state.exists() else "missing"

print(
    f"[iYali26] branch={branch}@{head} | uncommitted changes={changes} | "
    f"IYALI26_RESEARCH_ROOT={'set' if research else 'UNSET (builds fail closed)'} | "
    f"latest STATE.md entry: {latest}"
)
