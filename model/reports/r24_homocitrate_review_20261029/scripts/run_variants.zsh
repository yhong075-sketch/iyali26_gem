#!/bin/zsh
# Run r24_variants.py on a SHA-checked short copy of e1e7c868 (cobra rejects str paths of 260+ chars on macOS).
# Usage: zsh run_variants.zsh SHORT_COPY_DIR   (real HOME for the Gurobi licence; /opt/miniconda3/bin/python)
set -eu
WT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/nice-bouman-b161f4"
ACO="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/aconitase-aco2-split-1c492f"
RR="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
TASK="$WT/model/reports/r24_homocitrate_review_20261029"
short=$1
mkdir -p "$short"
src="$ACO/model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml"
cp "$src" "$short/e1e7c868.xml"
[[ "$(shasum -a 256 "$short/e1e7c868.xml" | cut -d' ' -f1)" == e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897 ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
shasum -a 256 "$src" "$short/e1e7c868.xml" "$medium" "$profile" "$explist" "$TASK/scripts/r24_variants.py" \
    "$WT/model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py" > "$TASK/logs/variant_inputs_sha256.txt"
cd "$WT/platform"
PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts" /opt/miniconda3/bin/python \
    "$TASK/scripts/r24_variants.py" "$short/e1e7c868.xml" "$medium" "$profile" "$explist" "$TASK/variants.json" \
    > "$TASK/logs/variants.log" 2>&1
