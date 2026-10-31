#!/bin/zsh
# Usage: zsh run_build.zsh SEED OUTPUT_MODEL LOG [EXTRA_FLAG ...]
# Runs the e1e7c868 argv (args_e1e7c868.txt) from platform/ with PYTHONHASHSEED=SEED, plus any extra flags.
set -u
WT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/trehalose-biomass-24e3c1"
HERE="${0:A:h}"
seed=$1; out=$2; log=$3; shift 3
args=("${(@f)$(cat "$HERE/args_e1e7c868.txt")}")
cd "$WT/platform"
PYTHONHASHSEED=$seed /opt/miniconda3/bin/python -m scripts.gem_annotate "${args[@]}" "$@" --output-model "$out" > "$log" 2>&1
echo "exit $?"
shasum -a 256 "$out"
