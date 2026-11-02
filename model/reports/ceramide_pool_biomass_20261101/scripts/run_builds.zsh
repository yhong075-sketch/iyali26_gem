#!/bin/zsh
# Control and ceramide-pool builds. Usage: run_builds.zsh <repo> <scratch_dir> <log_dir> <which: scratch|candidate>
# args_c1f70ac4.json holds c1f70ac4's argv from its .build.json (without __main__ and --output-model).
# scratch: control (no new flag, seed 0) and the new build at seed 1, both to the scratch dir.
# candidate: the new build at seed 0 to model/candidates.
set -u
REPO=$1; S=$2; L=$3
PY=/opt/miniconda3/bin/python
NAME=E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_20261101.xml
build() {  # label seed extra out
  local -a A
  A=("${(@f)$($PY -c 'import json,sys;print("\n".join(json.load(open(sys.argv[1]))))' $S/args_c1f70ac4.json)}")
  [[ -n "$3" ]] && A+=("$3")
  print -r -- "PYTHONHASHSEED=$2 $PY -m scripts.gem_annotate ${(q)A[@]} --output-model ${(q)4}" >> $L/build_commands.txt
  ( cd $REPO/platform && PYTHONHASHSEED=$2 $PY -m scripts.gem_annotate "${A[@]}" --output-model "$4" ) > $L/build_$1.log 2>&1
  echo "$1 exit=$?" >> $L/build_results.txt
}
if [[ $4 == scratch ]]; then
  build control_c1f70ac4 0 "" $S/build/control_c1f70ac4.xml
  build ceramide_pool_seed1 1 --ceramide-pool-biomass $S/build/ceramide_pool_seed1.xml
  shasum -a 256 $S/build/control_c1f70ac4.xml $S/build/ceramide_pool_seed1.xml >> $L/build_results.txt
else
  build ceramide_pool_seed0 0 --ceramide-pool-biomass $REPO/model/candidates/$NAME
  shasum -a 256 $REPO/model/candidates/$NAME >> $L/build_results.txt
fi
