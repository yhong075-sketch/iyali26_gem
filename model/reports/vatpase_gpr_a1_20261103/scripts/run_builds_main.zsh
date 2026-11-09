#!/bin/zsh
# Merge check on origin/main 9b6eafd8: control with e584151f's argv (flag off) and the same argv plus --vatpase-a1-gpr,
# both seed 0, both to the scratch dir (no candidate is kept). Usage: run_builds_main.zsh <repo> <scratch_dir> <log_dir>
set -u
REPO=$1; S=$2; L=$3
PY=/opt/miniconda3/bin/python
build() {  # label extra out
  local -a A
  A=("${(@f)$($PY -c 'import json,sys;print("\n".join(json.load(open(sys.argv[1]))))' $REPO/model/reports/vatpase_gpr_a1_20261103/scripts/args_e584151f.json)}")
  [[ -n "$2" ]] && A+=("$2")
  print -r -- "PYTHONHASHSEED=0 $PY -m scripts.gem_annotate ${(q)A[@]} --output-model ${(q)3}" >> $L/build_commands_main.txt
  ( cd $REPO/platform && PYTHONHASHSEED=0 $PY -m scripts.gem_annotate "${A[@]}" --output-model "$3" ) > $L/build_main_$1.log 2>&1
  echo "$1 exit=$?" >> $L/build_results_main.txt
}
build control_e584151f "" $S/build/control_e584151f.xml
build vatpase_a1_on_e584151f --vatpase-a1-gpr $S/build/vatpase_a1_on_e584151f.xml
shasum -a 256 $S/build/control_e584151f.xml $S/build/vatpase_a1_on_e584151f.xml >> $L/build_results_main.txt
