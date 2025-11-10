#!/bin/bash
#SBATCH --job-name=stage0_prep
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=60
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

container="${container:-/storage/singularity/mixtures.sif}"

# Use the directory where you invoked sbatch (Alderaan sets this)
export PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"

# Threads for math libs inside container
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "[Stage 0] host=$(hostname) user=$(whoami) date=$(date)"
echo "[Stage 0] PROJECT_ROOT=$PROJECT_ROOT"
echo "[Stage 0] Container=$container"

rscript_path="0_Prep/R Code/0_driver.R"
if [[ ! -f "$rscript_path" ]]; then
  echo "FATAL: '$rscript_path' not found under '$PROJECT_ROOT'"; ls -la; exit 2
fi

singularity exec "$container" Rscript "$rscript_path"
