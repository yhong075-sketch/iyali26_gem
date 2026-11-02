#!/bin/zsh
# Full pytest on origin/main (74146c9, extracted with git archive) and then on the combined branch, sequentially.
# Usage: run_pytest_pair.zsh <main_tree> <combined_repo> <log_dir>
set -u
V="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/workflow-automation-review-6765da/platform/.venv/bin/python"
export IYALI26_RESEARCH_ROOT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
cd "$1/platform" && "$V" -m pytest -q -rfE -p no:cacheprovider > "$3/pytest_before_main.txt" 2>&1; echo "before exit=$?" >> "$3/pytest_exit.txt"
cd "$2/platform" && "$V" -m pytest -q -rfE -p no:cacheprovider > "$3/pytest_after_combined.txt" 2>&1; echo "after exit=$?" >> "$3/pytest_exit.txt"
