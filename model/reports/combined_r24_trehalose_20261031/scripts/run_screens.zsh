#!/bin/zsh
# Certified-infeasible knockout screens (rules A and B) on 164542fb and on the combined candidate, one after the other.
# ko_screen_certified.py is copied unchanged from ../../infeasible_ko_rule_20261030/scripts (benchmark D3).
# Usage: zsh run_screens.zsh SHORT_COPY_DIR COMBINED_SHA256
# Short copies: cobra on macOS rejects str SBML paths of 260+ characters; each copy is SHA-checked before its run.
# Run with the real HOME (Gurobi licence) and /opt/miniconda3/bin/python (scipy present).
set -eu
WT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/essential-gene-prediction-e76be3"
RR="${IYALI26_RESEARCH_ROOT:?IYALI26_RESEARCH_ROOT must be set}"
TASK="$WT/model/reports/combined_r24_trehalose_20261031"
short=${1:A}; combined_sha=$2
mkdir -p "$short" "$TASK/logs"
C="$WT/model/candidates"
P=E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split
cp "$C/${P}_treh_biomass_v2_20261030.xml" "$short/164542fb.xml"
cp "$C/${P}_treh_biomass_r24_aco_20261031.xml" "$short/combined.xml"
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
untargeted="$WT/model/expected/untargeted_genes_out_of_scope.json"
solver="$WT/model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py"
shasum -a 256 "$short/164542fb.xml" "$short/combined.xml" "$medium" "$profile" "$explist" "$untargeted" \
    "$WT/model/expected/mtdna_genes_out_of_scope.json" "$TASK/scripts/ko_screen_certified.py" "$solver" \
    > "$TASK/logs/sweep_inputs_sha256.txt"
[[ "$(shasum -a 256 "$short/164542fb.xml" | cut -d' ' -f1)" == 164542fb6056f01e257ecfaed1b456961e8a07e46e9e4e1942f79489d5014beb ]]
[[ "$(shasum -a 256 "$short/combined.xml" | cut -d' ' -f1)" == "$combined_sha" ]]
[[ "$(shasum -a 256 "$medium" | cut -d' ' -f1)" == 6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45 ]]
[[ "$(shasum -a 256 "$profile" | cut -d' ' -f1)" == 35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383 ]]
[[ "$(shasum -a 256 "$explist" | cut -d' ' -f1)" == 1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b ]]
[[ "$(shasum -a 256 "$solver" | cut -d' ' -f1)" == 69c42690c231284768dc7f512832c327bc2fc2bce508bf78fa1946a75e7e72f0 ]]
cd "$WT/platform"
set +e
for label in 164542fb combined; do
  PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts" /opt/miniconda3/bin/python \
      "$TASK/scripts/ko_screen_certified.py" "$short/$label.xml" "$medium" "$profile" "$explist" "$untargeted" \
      "$TASK/ko_screen_$label.tsv" "$TASK/ko_screen_${label}_summary.json" > "$TASK/logs/ko_screen_$label.log" 2>&1
  echo "screen $label exit $?" >> "$TASK/logs/ko_screen_exit.txt"
done
