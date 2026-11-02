#!/bin/zsh
# Run r24_candidate_kos.py on a SHA-checked short copy of the candidate (cobra rejects str paths of 260+ chars on macOS).
# Usage: zsh run_kos.zsh SHORT_COPY_DIR CANDIDATE_SHA256   (real HOME for the Gurobi licence; /opt/miniconda3/bin/python)
set -eu
HERE="${0:A:h}"
WT="${HERE:h:h:h:h}"
RR="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
TASK="${HERE:h}"
short=$1; sha=$2
mkdir -p "$short"
src="$WT/model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_r24_aco_20261030.xml"
cp "$src" "$short/candidate.xml"
[[ "$(shasum -a 256 "$short/candidate.xml" | cut -d' ' -f1)" == "$sha" ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
shasum -a 256 "$src" "$short/candidate.xml" "$medium" "$profile" "$explist" "$HERE/r24_candidate_kos.py" \
    "$WT/model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py" > "$TASK/logs/ko_inputs_sha256.txt"
cd "$WT/platform"
PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts" /opt/miniconda3/bin/python \
    "$HERE/r24_candidate_kos.py" "$short/candidate.xml" "$medium" "$profile" "$explist" "$TASK/candidate_kos.json" \
    > "$TASK/logs/candidate_kos.log" 2>&1
