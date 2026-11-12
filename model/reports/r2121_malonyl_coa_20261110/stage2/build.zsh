#!/bin/zsh
# R2121 closure builds (TASK.md section 8). Usage: zsh build.zsh control|candidate|seed1
# Same argv as the e584151f .build.json, run from platform/ with miniconda Python 3.13.5 (as recorded there).
set -u
WT="/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.claude/worktrees/yali1-c06006g-gene-024dc8"
A="$IYALI26_RESEARCH_ROOT/artifacts/tasks/r2121_malonyl_coa_20261110/stage2"
NAME="E5_coq9_resp_pkg_r132_combined_dolichol_route_v2_um7_r2171_1to1_aars_thr_ser_split_r171_atp_synthase_gpr_r673_prolyl_split_aco_split_treh_biomass_r24_aco_ceramide_pool_bigg_annot_v2_mcnaugh_rxns_v3_r2121_close_20261112"
ARGS=(--research-root "$IYALI26_RESEARCH_ROOT" --offline --no-solve --coq9-curation qcycle --energy-candidate E5
  --coq9-biomass-dilution --coq9-respiratory-package --r132-mito-adh --restore-selection-cleared-gprs
  --restore-selection-overwritten-gprs --restore-selection-cpa-ura2 --remove-selection-tkl-fragment-gpr
  --remove-r153-peptidase-gpr --r2041-forward --r2190-ribulose-product --r1964-gpr --gpi-chain-links --gpi-coproducts
  --tpp-requirement --thiamine-phosphate-chemistry --thiamine-proton-lock --heme-iron-biomass --r2193-iron-gpr
  --lipoylation-module --acyl-acp-chemistry --cofactor-biomass-terms --r2135-acyl-coa-pool --r1726-cpt-open
  --dolichol-route --um7-nuclear-diffusion --r2171-pi-plc-1to1 --aars-thr-ser-split --r171-atp-synthase-gpr
  --r673-prolyl-split --aconitase-aco2-split --trehalose-biomass-removal --r24-homocitrate-aconitase
  --ceramide-pool-biomass --bigg-crosswalk-annotations --mcnaugh-reactions)
cd "$WT/platform"
case "$1" in
  control)   OUT="$A/builds/control_e584151f.xml"; SEED=0; EXTRA=() ;;
  candidate) OUT="../model/candidates/$NAME.xml";  SEED=0; EXTRA=(--r2121-malonyl-coa-close) ;;
  seed1)     OUT="$A/builds/seed1_r2121_close.xml"; SEED=1; EXTRA=(--r2121-malonyl-coa-close) ;;
  *) echo "unknown $1"; exit 2 ;;
esac
print -r -- "PYTHONHASHSEED=$SEED /opt/miniconda3/bin/python -m scripts.gem_annotate ${(q)ARGS[@]} ${(q)EXTRA[@]} --output-model ${(q)OUT}" >> "$A/logs/build_commands.txt"
PYTHONHASHSEED=$SEED /opt/miniconda3/bin/python -m scripts.gem_annotate "${ARGS[@]}" "${EXTRA[@]}" --output-model "$OUT" > "$A/logs/build_$1.log" 2>&1
rc=$?
echo "$1 exit=$rc" >> "$A/logs/build_results.txt"
[[ $rc -eq 0 ]] && shasum -a 256 "$OUT" >> "$A/logs/build_results.txt"
exit $rc
