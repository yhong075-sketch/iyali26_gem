#!/bin/zsh
# Usage: zsh run_build.zsh SEED OUTPUT_MODEL LOG [EXTRA_FLAG ...]
# Runs the 684f3752 argv (args_684f3752.txt) from platform/ with PYTHONHASHSEED=SEED, plus any extra flags.
set -u
HERE="${0:A:h}"
WT="${HERE:h:h:h:h}"
seed=$1; out=$2; log=$3; shift 3
args=("${(@f)$(cat "$HERE/args_684f3752.txt")}")
cd "$WT/platform"
PYTHONHASHSEED=$seed /opt/miniconda3/bin/python -m scripts.gem_annotate "${args[@]}" "$@" --output-model "$out" > "$log" 2>&1
echo "exit $?"
shasum -a 256 "$out"
