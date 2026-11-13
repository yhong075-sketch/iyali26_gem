#!/bin/zsh
# Full checked single-gene sweeps (ko_screen_certified.py, unchanged) on the SHA-checked control and candidate copies.
# Usage: zsh run_sweeps.zsh   (real HOME for the Gurobi licence; /opt/miniconda3/bin/python)
set -eu
HERE="${0:A:h}"; WT="${HERE:h:h:h:h}"; TASK="${HERE:h}"
RR="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"; A="$RR/artifacts/tasks/folinic_acid_module_20261112"
[[ "$(shasum -a 256 "$A/control684_seed0.xml" | cut -d' ' -f1)" == 684f3752764fb8f25219008966953863141724f34a451dffb77c235aa3cb2f62 ]]
[[ "$(shasum -a 256 "$A/candidate684_seed1.xml" | cut -d' ' -f1)" == 353cefdc16ef1dc510eae49c168fc999254c2652c3994287af45589621a90cd9 ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
untargeted="$WT/model/expected/untargeted_genes_out_of_scope.json"
screen="$WT/model/reports/infeasible_ko_rule_20261030/scripts/ko_screen_certified.py"
solve="$WT/model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py"
shasum -a 256 "$A/control684_seed0.xml" "$A/candidate684_seed1.xml" "$medium" "$profile" "$explist" "$untargeted" "$screen" "$solve" \
    > "$TASK/logs/sweep684_inputs_sha256.txt"
cd "$WT/platform"
export PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts"
for tag in control candidate; do
  [[ $tag == control ]] && m="$A/control684_seed0.xml" || m="$A/candidate684_seed1.xml"
  /opt/miniconda3/bin/python "$screen" "$m" "$medium" "$profile" "$explist" "$untargeted" \
      "$A/sweep684_$tag.tsv" "$TASK/sweep684_${tag}_summary.json" > "$TASK/logs/sweep684_$tag.log" 2>&1 &
done
wait
echo done
