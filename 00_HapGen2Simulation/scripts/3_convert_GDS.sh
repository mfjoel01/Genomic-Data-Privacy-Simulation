#!/bin/bash
#SBATCH --job-name=gds_convert
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

set -euo pipefail

# container path can be overridden at submit time: sbatch --export=ALL,container=/path/to.img 3_convert_GDS.sh
container="${container:-/storage/singularity/mixtures.sif}"

singularity exec "$container" Rscript hapgen2/convert_GDS.R
# (example with array arg, if you ever need it)
# singularity exec "$container" Rscript hapgen2/convert_GDS.R "$SLURM_ARRAY_TASK_ID"
