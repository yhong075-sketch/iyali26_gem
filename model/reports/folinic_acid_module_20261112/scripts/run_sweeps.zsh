#!/bin/zsh
# Full checked single-gene sweeps (ko_screen_certified.py, unchanged) on the SHA-checked control and candidate copies.
# Usage: zsh run_sweeps.zsh   (real HOME for the Gurobi licence; /opt/miniconda3/bin/python)
set -eu
HERE="${0:A:h}"; WT="${HERE:h:h:h:h}"; TASK="${HERE:h}"
RR="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"; A="$RR/artifacts/tasks/folinic_acid_module_20261112"
[[ "$(shasum -a 256 "$A/control_seed0.xml" | cut -d' ' -f1)" == e584151fc79475fc2daa59e76c8c84ebfd06b0f3cb1fbd75fbb69cc9f31039d2 ]]
[[ "$(shasum -a 256 "$A/candidate_seed1.xml" | cut -d' ' -f1)" == d4832b29294ff6506517cd09b4119fa71307af612a67b5395e2d6141b5e259e4 ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
untargeted="$WT/model/expected/untargeted_genes_out_of_scope.json"
screen="$WT/model/reports/infeasible_ko_rule_20261030/scripts/ko_screen_certified.py"
solve="$WT/model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py"
shasum -a 256 "$A/control_seed0.xml" "$A/candidate_seed1.xml" "$medium" "$profile" "$explist" "$untargeted" "$screen" "$solve" \
    > "$TASK/logs/sweep_inputs_sha256.txt"
cd "$WT/platform"
export PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts"
for tag in control candidate; do
  [[ $tag == control ]] && m="$A/control_seed0.xml" || m="$A/candidate_seed1.xml"
  /opt/miniconda3/bin/python "$screen" "$m" "$medium" "$profile" "$explist" "$untargeted" \
      "$A/sweep_$tag.tsv" "$TASK/sweep_${tag}_summary.json" > "$TASK/logs/sweep_$tag.log" 2>&1 &
done
wait
echo done
