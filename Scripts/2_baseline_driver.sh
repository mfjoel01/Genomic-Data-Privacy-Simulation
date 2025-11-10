#!/bin/bash
#SBATCH --job-name=stage2_baseline
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

container="${container:-/storage/singularity/mixtures.sif}"

export PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "[Stage 2] host=$(hostname) user=$(whoami) date=$(date)"
echo "[Stage 2] PROJECT_ROOT=$PROJECT_ROOT"
echo "[Stage 2] Container=$container"

rscript_path="2_Baseline/R Code/2_driver.R"
if [[ ! -f "$rscript_path" ]]; then
  echo "FATAL: '$rscript_path' not found under '$PROJECT_ROOT'"; ls -la; exit 2
fi

singularity exec "$container" Rscript "$rscript_path"
