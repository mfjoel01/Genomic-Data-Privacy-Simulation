#!/bin/bash
#SBATCH --job-name=stage3_dp
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=60
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

container="${container:-/storage/singularity/mixtures.sif}"

export PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"

# --- IMPORTANT: use single-threaded BLAS inside each forked worker ---
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1

echo "[Stage 3] host=$(hostname) user=$(whoami) date=$(date)"
echo "[Stage 3] PROJECT_ROOT=$PROJECT_ROOT"
echo "[Stage 3] Container=$container"

# Future/parallelly fork settings (for R's plan(multicore))
export R_FUTURE_GLOBALS_MAXSIZE=$((16 * 1024 * 1024 * 1024))  # 16 GiB
export R_PARALLELLY_FORK_ENABLE=true

# NOTE the space in directory name — keep quotes
rscript_path="3_DifferentialPrivacy/R Code/3_driver.R"
if [[ ! -f "$rscript_path" ]]; then
  echo "FATAL: '$rscript_path' not found under '$PROJECT_ROOT'"; ls -la; exit 2
fi

singularity exec "$container" Rscript "$rscript_path"
