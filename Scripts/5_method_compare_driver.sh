#!/bin/bash
#SBATCH --job-name=stage5_methodcmp
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=60
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

container="${container:-/storage/singularity/mixtures.sif}"

export PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"

# Single-threaded math in case of internal parallelism
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1

# Method comparison tuning (override via --export or env)
export MC_REPS="${MC_REPS:-10}"              # default 10 runs per method
export PGS_BETA_SUBSAMPLE="${PGS_BETA_SUBSAMPLE:-200}"
export FED_SITES="${FED_SITES:-5}"
export FED_SPLIT_SEED="${FED_SPLIT_SEED:-1}"

echo "[Stage 5] host=$(hostname) user=$(whoami) date=$(date)"
echo "[Stage 5] PROJECT_ROOT=$PROJECT_ROOT"
echo "[Stage 5] Container=$container"
echo "[Stage 5] MC_REPS=$MC_REPS PGS_BETA_SUBSAMPLE=$PGS_BETA_SUBSAMPLE FED_SITES=$FED_SITES FED_SPLIT_SEED=$FED_SPLIT_SEED"

rscript_path="5_MethodComparison/R Code/5_driver.R"
if [[ ! -f "$rscript_path" ]]; then
  echo "FATAL: '$rscript_path' not found under '$PROJECT_ROOT'"; ls -la; exit 2
fi

singularity exec "$container" Rscript "$rscript_path"
