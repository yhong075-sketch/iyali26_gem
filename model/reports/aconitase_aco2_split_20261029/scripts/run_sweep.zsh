#!/bin/zsh
# Checked single-gene knockout sweep, 059316ef vs the aconitase candidate (Task A medium + PO1f, Gurobi 1e-9).
# Usage: zsh run_sweep.zsh SHORT_COPY_DIR
# Short copies: cobra on macOS rejects str SBML paths of 260+ characters; each copy is SHA-checked before the run.
# Run with the real HOME (Gurobi licence) and /opt/miniconda3/bin/python (scipy present).
set -eu
WT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/aconitase-aco2-split-1c492f"
RR="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
TASK="$WT/model/reports/aconitase_aco2_split_20261029"
short=$1
mkdir -p "$short"
base="$WT/model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_notes_v3_20261029.xml"
cand="$WT/model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml"
cp "$base" "$short/059316ef.xml"
cp "$cand" "$short/aco_split.xml"
shasum -a 256 "$base" "$cand" "$short/059316ef.xml" "$short/aco_split.xml" > "$TASK/logs/sweep_short_copies_sha256.txt"
[[ "$(shasum -a 256 "$short/059316ef.xml" | cut -d' ' -f1)" == 059316ef34417c88aadba61cfd5606134d27879a842b3835d2cc7df917a9248d ]]
[[ "$(shasum -a 256 "$short/aco_split.xml" | cut -d' ' -f1)" == "$(shasum -a 256 "$cand" | cut -d' ' -f1)" ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
shasum -a 256 "$medium" "$profile" "$explist" "$TASK/scripts/ko_compare_two.py" > "$TASK/logs/sweep_inputs_sha256.txt"
cd "$WT/platform"
PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts" /opt/miniconda3/bin/python \
    "$TASK/scripts/ko_compare_two.py" "$short/059316ef.xml" "$short/aco_split.xml" "$medium" "$profile" "$explist" \
    "$TASK/ko_compare.tsv" "$TASK/ko_compare_summary.json" > "$TASK/logs/ko_compare.log" 2>&1
echo "sweep exit $?"
