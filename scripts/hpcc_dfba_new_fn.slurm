#!/bin/bash -l
#SBATCH --job-name=dfba_new_fn
#SBATCH --partition=batch
#SBATCH --account=iwheeldonlab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/rhome/yhong075/iyali26_gem_lipid_unlump_a014_runs/dfba_new_fn-%j.out
#SBATCH --error=/rhome/yhong075/iyali26_gem_lipid_unlump_a014_runs/dfba_new_fn-%j.err

set -uo pipefail

RUNS_ROOT=/rhome/yhong075/iyali26_gem_lipid_unlump_a014_runs
REPO="${REPO:-/rhome/yhong075/iyali26_gem_lipid_unlump_a014}"
PYTHON="${PYTHON:-/rhome/yhong075/iyali26_flow_phase2/.venv/bin/python}"
BASELINE_MODEL="${BASELINE_MODEL:-/rhome/yhong075/iyali26_gem_lipid_unlump_a014_inputs/model_bc2aac8fecd8f2f5.xml}"
CANDIDATE_MODEL="${CANDIDATE_MODEL:-$REPO/model_lipid_unlump_strict_sn_candidate_bc2aac8f_er_membrane_transport_metadata.xml}"
EXPECTED_BASELINE_SHA256="${EXPECTED_BASELINE_SHA256:-bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee}"
EXPECTED_CANDIDATE_SHA256="${EXPECTED_CANDIDATE_SHA256:-3b543b17b47b22a12a11745c3bc4ecf27e246bfa1d0f0d3865330742d1ac23f3}"
EXPECTED_RUNNER_SHA256="${EXPECTED_RUNNER_SHA256:-352f9f4980ecaf26cfce1eabb3b1c0f5661f4477ccd4c8fc0e9dfd63f0b2e962}"
EXPECTED_LOADER_SHA256="${EXPECTED_LOADER_SHA256:-f80f0195b993d44f6c61cb42b0987bc4491097de71080841e0e5e2662b4ce100}"
EXPECTED_FINGERPRINT_HELPER_SHA256="${EXPECTED_FINGERPRINT_HELPER_SHA256:-a3402f17eeb45ad14d4eb9e8374a55a4ce0fe41fe6c050d2a6756a60cb2042c3}"
EXPECTED_DFBA_MEDIUM_SHA256="${EXPECTED_DFBA_MEDIUM_SHA256:-7e189e312ddbae7cea711c5dfabbfbc6d0b10803f82a689e6d589d11b7ec5ec4}"
EXPECTED_EXPERIMENTAL_SHA256="${EXPECTED_EXPERIMENTAL_SHA256:-1e887f5ad4a95827a49b6c86894edaca410bdba3d264ff0d25193dedef3a659b}"
EXPECTED_GIT_COMMIT="${EXPECTED_GIT_COMMIT:-}"
EXPERIMENTAL_CSV="${EXPERIMENTAL_CSV:-/rhome/yhong075/iyali26_flow_phase2/data/essentiality/consensus_essential_genes.csv}"
DFBA_MEDIUM_CSV="${DFBA_MEDIUM_CSV:-$REPO/data/media/po1f_csm_leu_dfba.csv}"
DFBA_MODE="${DFBA_MODE:-single}"
DFBA_HOURS="${DFBA_HOURS:-24}"
DFBA_STEP_HOURS="${DFBA_STEP_HOURS:-0.1}"
DFBA_INITIAL_BIOMASS_GDW_L="${DFBA_INITIAL_BIOMASS_GDW_L:-0.05}"
DFBA_GROWTH_CUTOFF="${DFBA_GROWTH_CUTOFF:-0.01}"
DFBA_SOLVER="${DFBA_SOLVER:-glpk}"
DFBA_SHARD_COUNT="${DFBA_SHARD_COUNT:-32}"
EXPECTED_SHARD_COUNT="${EXPECTED_SHARD_COUNT:-32}"
RUN_ROOT="${RUN_ROOT:-}"
JOB_ID="${SLURM_JOB_ID:-manual}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
STATUS_PATH=""
PREPARE_STATUS=0

is_uint() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

prepare_output() {
  case "$DFBA_MODE" in
    single)
      OUTPUT_DIR="${OUTPUT_DIR:-$RUNS_ROOT/dfba_new_fn_${JOB_ID}}"
      mkdir -p "$OUTPUT_DIR" || PREPARE_STATUS=2
      STATUS_PATH="$OUTPUT_DIR/wrapper_status.json"
      ;;
    shard)
      if test -z "$RUN_ROOT" || ! is_uint "$TASK_ID" || ! is_uint "$DFBA_SHARD_COUNT"; then
        PREPARE_STATUS=2
        RUN_ROOT="${RUN_ROOT:-$RUNS_ROOT/dfba_new_fn_invalid_${JOB_ID}}"
        mkdir -p "$RUN_ROOT" || true
        STATUS_PATH="$RUN_ROOT/wrapper_status_${JOB_ID}.json"
        return
      fi
      OUTPUT_DIR="$RUN_ROOT/shards/$(printf '%03d' "$TASK_ID")"
      mkdir -p "$RUN_ROOT/shards" || PREPARE_STATUS=2
      if test "$PREPARE_STATUS" -eq 0 && ! mkdir "$OUTPUT_DIR"; then
        PREPARE_STATUS=2
      fi
      if test "$PREPARE_STATUS" -eq 0; then
        STATUS_PATH="$OUTPUT_DIR/wrapper_status.json"
      else
        STATUS_PATH="$RUN_ROOT/wrapper_status_shard_${TASK_ID}.json"
      fi
      ;;
    aggregate)
      if test -z "$RUN_ROOT"; then
        PREPARE_STATUS=2
        RUN_ROOT="$RUNS_ROOT/dfba_new_fn_invalid_${JOB_ID}"
      fi
      OUTPUT_DIR="${OUTPUT_DIR:-$RUN_ROOT}"
      test "$OUTPUT_DIR" = "$RUN_ROOT" || PREPARE_STATUS=2
      mkdir -p "$RUN_ROOT" || PREPARE_STATUS=2
      STATUS_PATH="$RUN_ROOT/wrapper_status.json"
      ;;
    *)
      PREPARE_STATUS=2
      RUN_ROOT="${RUN_ROOT:-$RUNS_ROOT/dfba_new_fn_invalid_${JOB_ID}}"
      mkdir -p "$RUN_ROOT" || true
      STATUS_PATH="$RUN_ROOT/wrapper_status_${JOB_ID}.json"
      ;;
  esac
}

check_contract() {
  local status=0
  test -n "$EXPECTED_GIT_COMMIT" || status=2
  test -n "$EXPECTED_EXPERIMENTAL_SHA256" || status=2
  test -n "$EXPERIMENTAL_CSV" || status=2
  test -n "$DFBA_MEDIUM_CSV" || status=2
  case "$DFBA_MODE" in single|shard|aggregate) ;; *) status=2 ;; esac
  test -z "$(git -C "$REPO" status --porcelain 2>/dev/null)" || status=2
  test "$(git -C "$REPO" rev-parse HEAD 2>/dev/null)" = "$EXPECTED_GIT_COMMIT" || status=2
  test -f "$REPO/scripts/dfba_new_fn_essentiality.py" || status=2
  test -f "$REPO/scripts/gem_annotate/validate_essential_genes.py" || status=2
  test -f "$REPO/scripts/lp_sn12_candidate.py" || status=2
  test -f "$BASELINE_MODEL" || status=2
  test -f "$CANDIDATE_MODEL" || status=2
  test -f "$EXPERIMENTAL_CSV" || status=2
  test -f "$DFBA_MEDIUM_CSV" || status=2
  if test -f "$REPO/scripts/dfba_new_fn_essentiality.py"; then
    test "$(sha256sum "$REPO/scripts/dfba_new_fn_essentiality.py" | awk '{print $1}')" = "$EXPECTED_RUNNER_SHA256" || status=2
  fi
  if test -f "$REPO/scripts/gem_annotate/validate_essential_genes.py"; then
    test "$(sha256sum "$REPO/scripts/gem_annotate/validate_essential_genes.py" | awk '{print $1}')" = "$EXPECTED_LOADER_SHA256" || status=2
  fi
  if test -f "$REPO/scripts/lp_sn12_candidate.py"; then
    test "$(sha256sum "$REPO/scripts/lp_sn12_candidate.py" | awk '{print $1}')" = "$EXPECTED_FINGERPRINT_HELPER_SHA256" || status=2
  fi
  if test -f "$BASELINE_MODEL"; then
    test "$(sha256sum "$BASELINE_MODEL" | awk '{print $1}')" = "$EXPECTED_BASELINE_SHA256" || status=2
  fi
  if test -f "$CANDIDATE_MODEL"; then
    test "$(sha256sum "$CANDIDATE_MODEL" | awk '{print $1}')" = "$EXPECTED_CANDIDATE_SHA256" || status=2
  fi
  if test -f "$DFBA_MEDIUM_CSV"; then
    test "$(sha256sum "$DFBA_MEDIUM_CSV" | awk '{print $1}')" = "$EXPECTED_DFBA_MEDIUM_SHA256" || status=2
  fi
  if test -f "$EXPERIMENTAL_CSV"; then
    test "$(sha256sum "$EXPERIMENTAL_CSV" | awk '{print $1}')" = "$EXPECTED_EXPERIMENTAL_SHA256" || status=2
  fi
  if test "$DFBA_MODE" = "shard" || test "$DFBA_MODE" = "aggregate"; then
    test "$DFBA_HOURS" = "24" || status=2
    test "$DFBA_STEP_HOURS" = "0.1" || status=2
    test "$DFBA_INITIAL_BIOMASS_GDW_L" = "0.05" || status=2
    test "$DFBA_GROWTH_CUTOFF" = "0.01" || status=2
    test "$DFBA_SOLVER" = "glpk" || status=2
    test "$DFBA_SHARD_COUNT" = "$EXPECTED_SHARD_COUNT" || status=2
    test -n "$RUN_ROOT" || status=2
  fi
  if test "$DFBA_MODE" = "shard"; then
    is_uint "$TASK_ID" || status=2
    is_uint "$DFBA_SHARD_COUNT" || status=2
    test "$TASK_ID" -lt "$DFBA_SHARD_COUNT" 2>/dev/null || status=2
    if test -n "${SLURM_ARRAY_TASK_COUNT:-}"; then
      test "$SLURM_ARRAY_TASK_COUNT" = "$DFBA_SHARD_COUNT" || status=2
    fi
  fi
  if test "$DFBA_MODE" = "aggregate"; then
    test -d "$RUN_ROOT/shards" || status=2
  fi
  return "$status"
}

write_wrapper_status() {
  local preflight_status="$1" runner_status="$2" postflight_status="$3"
  local summary_sha="" results_sha="" experimental_sha="" medium_sha=""
  local technical=false scientific_finding=false final_gate_deferred=false wrapper_passed=false
  if test -f "$OUTPUT_DIR/summary.json"; then
    summary_sha="$(sha256sum "$OUTPUT_DIR/summary.json" | awk '{print $1}')"
  fi
  if test -f "$OUTPUT_DIR/gene_results.tsv"; then
    results_sha="$(sha256sum "$OUTPUT_DIR/gene_results.tsv" | awk '{print $1}')"
  fi
  if test -f "$EXPERIMENTAL_CSV"; then
    experimental_sha="$(sha256sum "$EXPERIMENTAL_CSV" | awk '{print $1}')"
  fi
  if test -f "$DFBA_MEDIUM_CSV"; then
    medium_sha="$(sha256sum "$DFBA_MEDIUM_CSV" | awk '{print $1}')"
  fi
  if test "$DFBA_MODE" = "shard"; then
    final_gate_deferred=true
  fi
  if test "$DFBA_MODE" != "shard" && test "$runner_status" -eq 3; then
    scientific_finding=true
  fi
  if test "$preflight_status" -eq 0 && test "$postflight_status" -eq 0 && { test "$runner_status" -eq 0 || test "$scientific_finding" = true; }; then
    technical=true
  fi
  if test "$technical" = true && test "$runner_status" -eq 0; then
    wrapper_passed=true
  fi
  mkdir -p "$(dirname "$STATUS_PATH")"
  printf '{"job_id":"%s","array_job_id":"%s","array_task_id":"%s","mode":"%s","run_root":"%s","shard_count":"%s","git_commit":"%s","job_script_sha256":"%s","partition":"%s","account":"%s","cpus_per_task":"%s","memory_per_node":"%s","time_limit":"%s","experimental_sha256":"%s","dynamic_medium_sha256":"%s","preflight_exit_code":%d,"runner_exit_code":%d,"postflight_exit_code":%d,"technical_execution_completed":%s,"scientific_finding_recorded":%s,"final_scientific_gate_deferred_to_aggregate":%s,"wrapper_gate_passed":%s,"summary_sha256":"%s","gene_results_sha256":"%s","automatic_resubmit":false,"production_apply_performed":false}\n' \
    "$JOB_ID" "${SLURM_ARRAY_JOB_ID:-}" "$TASK_ID" "$DFBA_MODE" "$RUN_ROOT" "$DFBA_SHARD_COUNT" \
    "$(git -C "$REPO" rev-parse HEAD 2>/dev/null)" "$(sha256sum "$0" | awk '{print $1}')" \
    "${SLURM_JOB_PARTITION:-batch}" "${SLURM_JOB_ACCOUNT:-iwheeldonlab}" \
    "${SLURM_CPUS_PER_TASK:-1}" "${SLURM_MEM_PER_NODE:-4096}" "${SLURM_TIMELIMIT:-12:00:00}" \
    "$experimental_sha" "$medium_sha" "$preflight_status" "$runner_status" "$postflight_status" \
    "$technical" "$scientific_finding" "$final_gate_deferred" "$wrapper_passed" \
    "$summary_sha" "$results_sha" > "$STATUS_PATH"
}

prepare_output
if test "$PREPARE_STATUS" -ne 0; then
  write_wrapper_status "$PREPARE_STATUS" -1 -1
  exit "$PREPARE_STATUS"
fi

PREFLIGHT_STATUS=0
check_contract || PREFLIGHT_STATUS=$?
if test "$PREFLIGHT_STATUS" -ne 0; then
  write_wrapper_status "$PREFLIGHT_STATUS" -1 -1
  exit "$PREFLIGHT_STATUS"
fi

RUN_ARGS=(
  --baseline "$BASELINE_MODEL"
  --candidate "$CANDIDATE_MODEL"
  --experimental "$EXPERIMENTAL_CSV"
  --dynamic-medium "$DFBA_MEDIUM_CSV"
  --output-dir "$OUTPUT_DIR"
  --hours "$DFBA_HOURS"
  --step-hours "$DFBA_STEP_HOURS"
  --initial-biomass "$DFBA_INITIAL_BIOMASS_GDW_L"
  --growth-cutoff "$DFBA_GROWTH_CUTOFF"
  --solver "$DFBA_SOLVER"
)
if test "$DFBA_MODE" = "shard"; then
  RUN_ARGS+=(--shard-index "$TASK_ID" --shard-count "$DFBA_SHARD_COUNT")
elif test "$DFBA_MODE" = "aggregate"; then
  RUN_ARGS+=(--shard-count "$DFBA_SHARD_COUNT" --aggregate-shards "$RUN_ROOT/shards")
fi

"$PYTHON" "$REPO/scripts/dfba_new_fn_essentiality.py" "${RUN_ARGS[@]}"
RUNNER_STATUS=$?

POSTFLIGHT_STATUS=0
check_contract || POSTFLIGHT_STATUS=$?
write_wrapper_status 0 "$RUNNER_STATUS" "$POSTFLIGHT_STATUS"
if test "$POSTFLIGHT_STATUS" -ne 0; then
  exit "$POSTFLIGHT_STATUS"
fi
exit "$RUNNER_STATUS"
