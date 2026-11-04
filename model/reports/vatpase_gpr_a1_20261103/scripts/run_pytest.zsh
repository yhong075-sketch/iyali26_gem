#!/bin/zsh
# Full pytest in one tree. Usage: run_pytest.zsh <tree> <log_file>
# Before: origin/main b57272d3 extracted with git archive (no edits). After: this worktree with the change.
set -u
V="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/workflow-automation-review-6765da/platform/.venv/bin/python"
export IYALI26_RESEARCH_ROOT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
cd "$1/platform" && "$V" -m pytest -q -rfE -p no:cacheprovider > "$2" 2>&1; echo "$(basename $2) exit=$?" >> "$(dirname $2)/pytest_exit.txt"
