#!/bin/zsh
# Certified-infeasible knockout screen on e1e7c868 (Task A medium + PO1f, Gurobi 1e-9), rules A and B.
# Usage: zsh run_screen.zsh SHORT_COPY_DIR
# Short copy: cobra on macOS rejects str SBML paths of 260+ characters; the copy is SHA-checked before the run.
# Run with the real HOME (Gurobi licence) and /opt/miniconda3/bin/python (scipy present).
set -eu
WT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/cx34-infeasible-discussion-7b0c5b"
RR="${IYALI26_RESEARCH_ROOT:?IYALI26_RESEARCH_ROOT must be set}"
TASK="$WT/model/reports/infeasible_ko_rule_20261030"
short=${1:A}
mkdir -p "$short" "$TASK/logs"
model="$WT/model/candidates/E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_20261029.xml"
cp "$model" "$short/e1e7c868.xml"
shasum -a 256 "$model" "$short/e1e7c868.xml" > "$TASK/logs/sweep_short_copy_sha256.txt"
[[ "$(shasum -a 256 "$short/e1e7c868.xml" | cut -d' ' -f1)" == e1e7c8686c1fbbbbe1ab949000e873cd5f2d7089bc008a4ec129f7ba8dcdf897 ]]
medium="$WT/model/conditions/media/po1f_csm_leu_w29uptake_20261012.csv"
profile="$RR/state/strain_profiles/po1f_sd_leu.json"
explist="$RR/state/essentiality/repository/consensus_essential_genes.csv"
untargeted="$TASK/scripts/untargeted_genes_out_of_scope.json"
solver="$WT/model/reports/combined_tpp_heme_r2193_20261022/scripts/checked_solve.py"
shasum -a 256 "$medium" "$profile" "$explist" "$untargeted" "$WT/model/expected/mtdna_genes_out_of_scope.json" \
    "$TASK/scripts/ko_screen_certified.py" "$solver" > "$TASK/logs/sweep_inputs_sha256.txt"
# Inputs pinned in TASK.md section 1
[[ "$(shasum -a 256 "$medium" | cut -d' ' -f1)" == 6eca70fd2557f014ab8fa6ed3609ebf5fc181a0a76c131d6a377b14febb08f45 ]]
[[ "$(shasum -a 256 "$profile" | cut -d' ' -f1)" == 35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383 ]]
[[ "$(shasum -a 256 "$explist" | cut -d' ' -f1)" == 1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b ]]
[[ "$(shasum -a 256 "$solver" | cut -d' ' -f1)" == 69c42690c231284768dc7f512832c327bc2fc2bce508bf78fa1946a75e7e72f0 ]]
cd "$WT/platform"
set +e
PYTHONPATH="$PWD:../model/reports/combined_tpp_heme_r2193_20261022/scripts" /opt/miniconda3/bin/python \
    "$TASK/scripts/ko_screen_certified.py" "$short/e1e7c868.xml" "$medium" "$profile" "$explist" "$untargeted" \
    "$TASK/ko_screen.tsv" "$TASK/ko_screen_summary.json" > "$TASK/logs/ko_screen.log" 2>&1
code=$?
echo "sweep exit $code"
exit $code
