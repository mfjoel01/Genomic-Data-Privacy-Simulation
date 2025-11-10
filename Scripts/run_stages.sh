#!/usr/bin/env bash
# Run pipeline stages x..y with Slurm dependencies.
# Defaults to 1 -> 5 (includes Stage 5 Method Comparison).

set -euo pipefail

# ------------- Root detection -------------
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd -P)"

find_root() {
  local d="$1"
  while [ "$d" != "/" ] && [ -n "$d" ]; do
    if [ -f "$d/R/config.R" ] && [ -d "$d/Scripts" ]; then
      printf '%s\n' "$(cd "$d" && pwd -P)"
      return 0
    fi
    d="$(cd "$d/.." && pwd -P)"
  done
  return 1
}

if [ -n "${PROJECT_ROOT:-}" ] && [ -f "${PROJECT_ROOT}/R/config.R" ] && [ -d "${PROJECT_ROOT}/Scripts" ]; then
  PROJECT_ROOT="$(cd "$PROJECT_ROOT" && pwd -P)"
else
  if ! PROJECT_ROOT="$(find_root "$SCRIPT_DIR")"; then
    echo "ERROR: Could not locate project root (need R/config.R and Scripts/) starting from: $SCRIPT_DIR" >&2
    exit 2
  fi
fi

cd "$PROJECT_ROOT"
echo "[run] Project root: $PROJECT_ROOT"
echo "[run] Invoked from: $SCRIPT_DIR"

# ------------- Defaults & CLI -------------
FROM=1
TO=5

BMI_MODE="${BMI_MODE:-}"          # optional; Stage 1
CONTAINER="${container:-}"        # optional; passed via --export as 'container=...'

# Stage 4 (Federated) env knobs
RUN_FED_SWEEP="${RUN_FED_SWEEP:-}"      # leave empty -> driver default (sweep ON)
FED_SITES_GRID="${FED_SITES_GRID:-}"    # e.g. "2,3,4,5,6,8,10"
FED_SITES="${FED_SITES:-}"              # single-K (still computed even if sweep runs)
FED_SPLIT_SEED="${FED_SPLIT_SEED:-}"    # deterministic split

DRY_RUN=0

usage() {
  echo "Usage:"
  echo "  run_stages.sh [FROM TO] [options]"
  echo "  run_stages.sh --from N --to M [options]"
  echo "  run_stages.sh --1to5        # 1..5"
  echo "  run_stages.sh --all         # 0..5"
  echo "  run_stages.sh --skip-prep   # bump FROM to >=1"
  echo
  echo "Options:"
  echo "  -f, --from N           First stage (0..5). Default: 1"
  echo "  -t, --to M             Last stage  (0..5). Default: 5"
  echo "  -b, --BMI_MODE VAL     Stage 1 BMI_MODE (PGS|SPARSE|NULL)"
  echo "  -c, --container IMG    Singularity image (exported as container=...)"
  echo "  -n, --dry-run          Print sbatch commands without submitting"
  echo
  echo "Stage 4 (Federated) options:"
  echo "      --no-fed-sweep     Set RUN_FED_SWEEP=0"
  echo "  -g, --fed-grid LIST    Comma list for K (e.g., 2,3,4,5,6,8,10)"
  echo "      --fed-sites K      Single K (FED_SITES)"
  echo "      --fed-seed N       FED_SPLIT_SEED"
}

is_stage_num() { case "$1" in 0|1|2|3|4|5) return 0;; *) return 1;; esac; }

# Positional FROM TO if first args are numbers
if [ $# -ge 1 ] && is_stage_num "${1#-}"; then
  FROM="$1"; shift
  if [ $# -ge 1 ] && is_stage_num "${1#-}"; then
    TO="$1"; shift
  fi
fi

# Flags
while [ $# -gt 0 ]; do
  case "$1" in
    -f|--from)        FROM="$2"; shift 2 ;;
    -t|--to)          TO="$2"; shift 2 ;;
    -b|--BMI_MODE)    BMI_MODE="$2"; shift 2 ;;
    -c|--container)   CONTAINER="$2"; shift 2 ;;
    -n|--dry-run|--dryrun) DRY_RUN=1; shift ;;
    --1to5|--15)      FROM=1; TO=5; shift ;;
    --1to4|--14)      FROM=1; TO=4; shift ;;
    --all|--0to5)     FROM=0; TO=5; shift ;;
    --skip-prep)      [ "$FROM" -lt 1 ] && FROM=1; shift ;;
    --no-fed-sweep)   RUN_FED_SWEEP="0"; shift ;;
    -g|--fed-grid)    FED_SITES_GRID="$2"; shift 2 ;;
    --fed-sites)      FED_SITES="$2"; shift 2 ;;
    --fed-seed)       FED_SPLIT_SEED="$2"; shift 2 ;;
    -h|--help)        usage; exit 0 ;;
    *) echo "Unknown argument: $1"; usage; exit 2 ;;
  esac
done

# ------------- Validation -------------
if ! is_stage_num "$FROM" || ! is_stage_num "$TO"; then
  echo "ERROR: FROM/TO must be integers in [0,5]. Got FROM='$FROM' TO='$TO'." >&2
  exit 2
fi
if [ "$FROM" -gt "$TO" ]; then
  echo "ERROR: FROM must be <= TO (FROM=$FROM, TO=$TO)." >&2
  exit 2
fi
if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: 'sbatch' not found in PATH." >&2
  exit 2
fi
if [ -n "$CONTAINER" ] && printf '%s' "$CONTAINER" | grep -qE '[[:space:]]'; then
  echo "ERROR: --container path contains spaces; cannot safely pass via --export." >&2
  exit 2
fi

# ------------- Stage map (indexed arrays) -------------
STAGE_NAME=()
STAGE_SCRIPT=()

STAGE_NAME[0]="Stage 0 – Prep"
STAGE_NAME[1]="Stage 1 – Phenotype"
STAGE_NAME[2]="Stage 2 – Baseline"
STAGE_NAME[3]="Stage 3 – Differential Privacy"
STAGE_NAME[4]="Stage 4 – Federated Learning"
STAGE_NAME[5]="Stage 5 – Method Comparison"

STAGE_SCRIPT[0]="Scripts/0_prep_driver.sh"
STAGE_SCRIPT[1]="Scripts/1_pheno_driver.sh"
STAGE_SCRIPT[2]="Scripts/2_baseline_driver.sh"
STAGE_SCRIPT[3]="Scripts/3_dp_driver.sh"
STAGE_SCRIPT[4]="Scripts/4_federated_driver.sh"
STAGE_SCRIPT[5]="Scripts/5_method_compare_driver.sh"

s="$FROM"
while [ "$s" -le "$TO" ]; do
  f="${STAGE_SCRIPT[$s]}"
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing script for stage $s: $f" >&2
    exit 2
  fi
  s=$((s+1))
done

# ------------- Helper to submit a stage -------------
submit_stage() {
  # args: stage dep_jobid
  local stage="$1"
  local dep_jid="${2:-}"
  local script="${STAGE_SCRIPT[$stage]}"

  local export_vars="ALL"
  [ -n "$BMI_MODE" ]       && [ "$stage" -eq 1 ] && export_vars+=",BMI_MODE=${BMI_MODE}"
  if [ "$stage" -eq 4 ]; then
    [ -n "$RUN_FED_SWEEP" ]  && export_vars+=",RUN_FED_SWEEP=${RUN_FED_SWEEP}"
    [ -n "$FED_SITES_GRID" ] && export_vars+=",FED_SITES_GRID=${FED_SITES_GRID}"
    [ -n "$FED_SITES" ]      && export_vars+=",FED_SITES=${FED_SITES}"
    [ -n "$FED_SPLIT_SEED" ] && export_vars+=",FED_SPLIT_SEED=${FED_SPLIT_SEED}"
  fi
  [ -n "$CONTAINER" ] && export_vars+=",container=${CONTAINER}"

  local dep_arg=()
  [ -n "$dep_jid" ] && dep_arg=( "--dependency=afterok:${dep_jid}" )

  if [ "$DRY_RUN" -eq 1 ]; then
    echo "[DRY-RUN] sbatch ${dep_arg[*]:-} --export=${export_vars} ${script}"
    # Fake job id for chaining demo
    echo "DRY${stage}"
    return 0
  fi

  local out jid
  if ! out="$(sbatch --parsable ${dep_arg:+${dep_arg[*]}} --export="${export_vars}" "${script}")"; then
    echo "ERROR: sbatch submission failed for stage $stage (${script})." >&2
    return 1
  fi
  jid="${out%%;*}"
  if ! printf '%s' "$jid" | grep -qE '^[0-9]+$'; then
    echo "ERROR: Could not parse job ID from sbatch output: $out" >&2
    return 1
  fi
  printf '%s\n' "$jid"
}

# ------------- Submit chained jobs -------------
declare -A JIDS
prev=""

s="$FROM"
while [ "$s" -le "$TO" ]; do
  echo "[run] Submitting ${STAGE_NAME[$s]}  ->  ${STAGE_SCRIPT[$s]}"
  if [ "$s" -eq 1 ] && [ -n "$BMI_MODE" ]; then
    echo "      BMI_MODE=${BMI_MODE}"
  fi
  if [ "$s" -eq 4 ]; then
    [ -n "$RUN_FED_SWEEP" ]  && echo "      RUN_FED_SWEEP=${RUN_FED_SWEEP}"
    [ -n "$FED_SITES_GRID" ] && echo "      FED_SITES_GRID=${FED_SITES_GRID}"
    [ -n "$FED_SITES" ]      && echo "      FED_SITES=${FED_SITES}"
    [ -n "$FED_SPLIT_SEED" ] && echo "      FED_SPLIT_SEED=${FED_SPLIT_SEED}"
  fi
  [ -n "$CONTAINER" ] && echo "      container override: ${CONTAINER}"

  jid="$(submit_stage "$s" "$prev")" || exit 1
  JIDS[$s]="$jid"
  prev="$jid"
  s=$((s+1))
done

# ------------- Summary -------------
summary="Jobs:"
s="$FROM"
while [ "$s" -le "$TO" ]; do
  summary+=" S${s}=${JIDS[$s]}"
  s=$((s+1))
done
echo "$summary"
echo "[run] Submitted ${FROM}..${TO} with afterok chaining. Tail job: ${JIDS[$TO]}"
