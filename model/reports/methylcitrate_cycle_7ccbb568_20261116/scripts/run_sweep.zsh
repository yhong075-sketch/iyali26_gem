#!/bin/zsh
# Full certified single-gene sweep (ko_screen_certified.py, unchanged) on the SHA-checked candidate.
# Usage: zsh run_sweep.zsh CANDIDATE_SHA256   (real HOME for the Gurobi licence; /opt/miniconda3/bin/python)
set -eu
HERE="${0:A:h}"; WT="${HERE:h:h:h:h}"; TASK="${HERE:h}"
RR="${IYALI26_RESEARCH_ROOT:?}"; A="$RR/artifacts/tasks/methylcitrate_cycle_7ccbb568_20261116"
m="$WT/model/candidates/E5_7ccbb568_methylcitrate_20261116.xml"
[[ "$(shasum -a 256 "$m" | cut -d' ' -f1)" == "$1" ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
untargeted="$WT/model/expected/untargeted_genes_out_of_scope.json"
screen="$WT/model/reports/infeasible_ko_rule_20261030/scripts/ko_screen_certified.py"
solve="$WT/model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py"
shasum -a 256 "$m" "$medium" "$profile" "$explist" "$untargeted" "$WT/model/expected/mtdna_genes_out_of_scope.json" \
    "$screen" "$solve" > "$TASK/logs/sweep_inputs_sha256.txt"
[[ "$(shasum -a 256 "$medium" | cut -d' ' -f1)" == 6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45 ]]
[[ "$(shasum -a 256 "$profile" | cut -d' ' -f1)" == 35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383 ]]
[[ "$(shasum -a 256 "$explist" | cut -d' ' -f1)" == 1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b ]]
[[ "$(shasum -a 256 "$screen" | cut -d' ' -f1)" == 313d5d0125b0fdaf0bdcd5e353298c5881f543404407ce31c9bd2056ed32dbc7 ]]
[[ "$(shasum -a 256 "$solve" | cut -d' ' -f1)" == 69c42690c231284768dc7f512832c327bc2fc2bce508bf78fa1946a75e7e72f0 ]]
cd "$WT/platform"
export PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts"
/opt/miniconda3/bin/python "$screen" "$m" "$medium" "$profile" "$explist" "$untargeted" \
    "$A/sweep_candidate.tsv" "$TASK/sweep_candidate_summary.json" > "$TASK/logs/sweep_candidate.log" 2>&1
echo "sweep exit $?"
