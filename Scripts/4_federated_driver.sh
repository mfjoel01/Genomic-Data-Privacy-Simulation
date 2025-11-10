#!/bin/bash
#SBATCH --job-name=stage4_federated
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=60
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
# (optional) #SBATCH --time=24:00:00

set -euo pipefail

# -------- Container --------
container="${container:-/storage/singularity/mixtures.sif}"

# -------- Project root (Alderaan sets SLURM_SUBMIT_DIR) --------
export PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$PROJECT_ROOT"

# -------- IMPORTANT: single-threaded BLAS inside each R worker --------
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1

# -------- Parallel worker cap for Stage 4 R code --------
# 4b_fed_gwas.R honors FED_CORES; default = SLURM_CPUS_PER_TASK
export FED_CORES="${FED_CORES:-${SLURM_CPUS_PER_TASK:-8}}"

# -------- Federated settings --------
export FED_SITES="${FED_SITES:-5}"                          # single-K run (also used by sweep)
export FED_SPLIT_SEED="${FED_SPLIT_SEED:-1}"                # deterministic split across steps
export RUN_FED_SWEEP="${RUN_FED_SWEEP:-1}"                  # sweep runs by default
export FED_SITES_GRID="${FED_SITES_GRID:-2,3,4,5,6,8,10}"   # grid for sweep

echo "[Stage 4] host=$(hostname) user=$(whoami) date=$(date)"
echo "[Stage 4] PROJECT_ROOT=$PROJECT_ROOT"
echo "[Stage 4] Container=$container"
echo "[Stage 4] FED_SITES=$FED_SITES FED_SPLIT_SEED=$FED_SPLIT_SEED"
echo "[Stage 4] RUN_FED_SWEEP=$RUN_FED_SWEEP FED_SITES_GRID=$FED_SITES_GRID"
echo "[Stage 4] FED_CORES=$FED_CORES  (BLAS threads pinned: OMP/MKL/OPENBLAS/BLIS=1)"

rscript_path="4_FederatedLearning/R Code/4_driver.R"
if [[ ! -f "$rscript_path" ]]; then
  echo "FATAL: '$rscript_path' not found under '$PROJECT_ROOT'"; ls -la; exit 2
fi

singularity exec "$container" Rscript "$rscript_path"
