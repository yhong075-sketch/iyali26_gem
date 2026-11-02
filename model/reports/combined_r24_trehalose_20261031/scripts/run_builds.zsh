#!/bin/zsh
# Control and combined builds on the merged code. Usage: run_builds.zsh <repo> <scratch_dir> <log_dir>
# args_164542fb.json / args_f505571d.json hold each parent's argv from its .build.json (without __main__ and
# --output-model). The combined build is args_164542fb + --r24-homocitrate-aconitase.
set -u
REPO=$1; S=$2; L=$3
PY=/opt/miniconda3/bin/python
NAME=E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_20261031.xml
build() {  # label seed argsfile extra out
  local -a A
  A=("${(@f)$($PY -c 'import json,sys;print("\n".join(json.load(open(sys.argv[1]))))' $3)}")
  [[ -n "$4" ]] && A+=("$4")
  print -r -- "PYTHONHASHSEED=$2 $PY -m scripts.gem_annotate ${(q)A[@]} --output-model ${(q)5}" >> $L/build_commands.txt
  ( cd $REPO/platform && PYTHONHASHSEED=$2 $PY -m scripts.gem_annotate "${A[@]}" --output-model "$5" ) > $L/build_$1.log 2>&1
  echo "$1 exit=$?" >> $L/build_results.txt
}
build control_164542fb 0 $S/args_164542fb.json "" $S/build/control_164542fb.xml
build control_f505571d 0 $S/args_f505571d.json "" $S/build/control_f505571d.xml
build combined_seed0 0 $S/args_164542fb.json --r24-homocitrate-aconitase $REPO/model/candidates/$NAME
build combined_seed1 1 $S/args_164542fb.json --r24-homocitrate-aconitase $S/build/combined_seed1.xml
shasum -a 256 $S/build/control_164542fb.xml $S/build/control_f505571d.xml $REPO/model/candidates/$NAME $S/build/combined_seed1.xml >> $L/build_results.txt
