#!/bin/bash
#SBATCH --job-name=hapgen2_10k
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
# no --time (uses partition default)

set -euo pipefail

# --- your setup ---
container="${container:-/storage/singularity/mixtures.sif}"
REF_DIR="ref"
HAP="${REF_DIR}/chr22_EUR_ref.hap"          # or .hap.gz present next to it
LEG="${REF_DIR}/chr22_EUR_ref.legend"       # or .legend.gz present next to it
MAP="input/1000G/genetic_map_chr22_combined_b37.txt"
DUMMY_FILE="${REF_DIR}/dummy_pos.txt"
OUT_DIR="hapgen2"
mkdir -p "${OUT_DIR}"

# ensure plain-text inputs exist (use local copies if only gz present)
[[ -s "$HAP" ]] || zcat "${HAP}.gz" > "$HAP"
[[ -s "$LEG" ]] || zcat "${LEG}.gz" > "$LEG"

# pick a valid dummy position: prefer your file, else 1st SNP in legend
DUMMY_POS=$( [[ -s "$DUMMY_FILE" ]] && cat "$DUMMY_FILE" || awk 'NR==2{print $2; exit}' "$LEG" )

# simulate 10,000 controls; add neutral dummy -dl to satisfy strict builds
singularity exec "$container" hapgen2 \
  -m "$MAP" -l "$LEG" -h "$HAP" \
  -o "${OUT_DIR}/chr22_EUR_sim_10k.gz" \
  -n 10000 0 \
  -dl "$DUMMY_POS" 0 1 1

echo "Done. Outputs in: ${OUT_DIR}"