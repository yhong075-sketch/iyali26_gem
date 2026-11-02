#!/bin/zsh
# Checked single-gene knockout sweep, e1e7c868 vs the trehalose candidate (Task A medium + PO1f, Gurobi 1e-9).
# Usage: zsh run_sweep.zsh SHORT_COPY_DIR
# Short copies: cobra on macOS rejects str SBML paths of 260+ characters; each copy is SHA-checked before the run.
# Run with the real HOME (Gurobi licence) and /opt/miniconda3/bin/python (scipy present).
set -eu
WT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/trehalose-biomass-24e3c1"
RR="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem_research"
TASK="$WT/model/reports/trehalose_biomass_20261030"
short=$1
mkdir -p "$short"
base="$WT/model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml"
cand="$WT/model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_20261030.xml"
cp "$base" "$short/e1e7c868.xml"
cp "$cand" "$short/treh_biomass.xml"
shasum -a 256 "$base" "$cand" "$short/e1e7c868.xml" "$short/treh_biomass.xml" > "$TASK/logs/sweep_short_copies_sha256.txt"
[[ "$(shasum -a 256 "$short/e1e7c868.xml" | cut -d' ' -f1)" == e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897 ]]
[[ "$(shasum -a 256 "$short/treh_biomass.xml" | cut -d' ' -f1)" == "$(shasum -a 256 "$cand" | cut -d' ' -f1)" ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
shasum -a 256 "$medium" "$profile" "$explist" "$TASK/scripts/ko_compare_two.py" > "$TASK/logs/sweep_inputs_sha256.txt"
cd "$WT/platform"
PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts" /opt/miniconda3/bin/python \
    "$TASK/scripts/ko_compare_two.py" "$short/e1e7c868.xml" "$short/treh_biomass.xml" "$medium" "$profile" "$explist" \
    "$TASK/ko_compare.tsv" "$TASK/ko_compare_summary.json" > "$TASK/logs/ko_compare.log" 2>&1
echo "sweep exit $?"
