#!/bin/zsh
# Full pytest on HEAD b83e13a (extracted with git archive, no edits) and then on this worktree with the change.
# Usage: run_pytest_pair.zsh <before_tree> <worktree> <log_dir>
set -u
V="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/workflow-automation-review-6765da/platform/.venv/bin/python"
export IYALI26_RESEARCH_ROOT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
cd "$1/platform" && "$V" -m pytest -q -rfE -p no:cacheprovider > "$3/pytest_before_b83e13a.txt" 2>&1; echo "before exit=$?" >> "$3/pytest_exit.txt"
cd "$2/platform" && "$V" -m pytest -q -rfE -p no:cacheprovider > "$3/pytest_after.txt" 2>&1; echo "after exit=$?" >> "$3/pytest_exit.txt"
