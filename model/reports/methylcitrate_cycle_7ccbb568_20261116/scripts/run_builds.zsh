#!/bin/zsh
# Control (no new flag, seed 0), candidate seed 0 (formal file in model/candidates/), candidate seed 1 (artifacts).
set -u
HERE="${0:A:h}"; WT="${HERE:h:h:h:h}"; TASK="${HERE:h}"
A="${IYALI26_RESEARCH_ROOT:?}/artifacts/tasks/methylcitrate_cycle_7ccbb568_20261116"
zsh "$HERE/run_build.zsh" 0 "$A/control_seed0.xml" "$TASK/logs/build_control_seed0.log"
zsh "$HERE/run_build.zsh" 0 "$WT/model/candidates/E5_7ccbb568_methylcitrate_20261116.xml" \
    "$TASK/logs/build_candidate_seed0.log" --methylcitrate-cycle
zsh "$HERE/run_build.zsh" 1 "$A/candidate_seed1.xml" "$TASK/logs/build_candidate_seed1.log" --methylcitrate-cycle
