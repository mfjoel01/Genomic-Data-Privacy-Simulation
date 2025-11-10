#!/bin/bash
#SBATCH --job-name=gen2vcf_chr22
#SBATCH --partition=math-alderaan
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
set -euo pipefail

# --- inputs (edit if paths differ) ---
in_gen="hapgen2/chr22_EUR_sim_10k.controls.gen.gz"
in_sample="hapgen2/chr22_EUR_sim_10k.controls.sample"
chr="22"

# --- environment / container ---
container="${container:-/storage/singularity/mixtures.sif}"

# --- outputs ---
OUT_DIR="vcf"
mkdir -p "${OUT_DIR}"
prefix="${OUT_DIR}/chr${chr}_EUR_sim_10k.controls"
vcf_gz="${prefix}.vcf.gz"

echo "[INFO] Inputs:"
echo "  GEN    : $in_gen"
echo "  SAMPLE : $in_sample"
echo "  CHR    : $chr"
echo "  OUT    : $vcf_gz"

# Helper: run a command inside the container
cexec() { singularity exec "$container" bash -lc "$*"; }
have()  { cexec "command -v \"$1\" >/dev/null 2>&1"; }

# ---------- Path A: qctool (preferred, preserves rsIDs, fills CHROM) ----------
if have qctool && have bgzip && have tabix; then
  echo "[INFO] Using qctool (rsIDs kept; -assume-chromosome ${chr})"
  cexec "
    qctool -g \"$in_gen\" -s \"$in_sample\" \
      -assume-chromosome \"$chr\" \
      -og '${prefix}.vcf' -ofiletype vcf
    bgzip -f '${prefix}.vcf'
    tabix -p vcf '$vcf_gz'
  "
  echo "[OK] Wrote $vcf_gz and index (qctool path)"
  exit 0
fi

# ---------- Path B: bcftools (requires column2 = CHROM:POS_REF_ALT) ----------
if have bcftools && have bgzip && have tabix; then
  echo "[INFO] Using bcftools convert --gensample2vcf"

  fixed_gen="${prefix}.fixed.for_bcftools.gen.gz"
  map_tsv="${prefix}.pos2rs.map.tsv.gz"   # CHROM  POS  rsID

  echo "[INFO] Rewriting .gen so: COL1=rsID, COL2=CHR:POS_REF_ALT"
  # Original (Oxford GEN): ID1 RSID POS REF ALT ...
  # We set:      ID1:=RSID     ID2:=CHR:POS_REF_ALT   (bcftools requires ID2 in that form)
  cexec "
    zcat \"$in_gen\" \
    | awk -v CHR=\"$chr\" 'BEGIN{OFS=\" \"}
        { id1=\$2; id2=CHR\":\"\$3\"_\"\$4\"_\"\$5; \$1=id1; \$2=id2; print }
      ' \
    | bgzip -c > \"$fixed_gen\"
  "

  echo "[INFO] Building CHROM:POS -> rsID map for ID replacement"
  cexec "
    zcat \"$fixed_gen\" \
    | awk 'BEGIN{OFS=\"\t\"}
        {
          # COL1=rsID, COL2=CHR:POS_REF_ALT
          split(\$2, a, \":\"); CHR=a[1];
          split(a[2], b, \"_\"); POS=b[1];
          print CHR, POS, \$1
        }
      ' \
    | bgzip -c > \"$map_tsv\"
    tabix -s 1 -b 2 -e 2 \"$map_tsv\"
  "

  echo "[INFO] Converting to VCF"
  tmp_vcf="${prefix}.tmp.idpos.vcf.gz"
  cexec "
    bcftools convert --gensample2vcf \"$fixed_gen\",\"$in_sample\" \
      -Oz -o \"$tmp_vcf\" --threads ${SLURM_CPUS_PER_TASK:-1}
    tabix -p vcf \"$tmp_vcf\"
  "

  echo "[INFO] Replacing VCF ID with rsID via bcftools annotate"
  cexec "
    bcftools annotate -a \"$map_tsv\" -c CHROM,POS,ID \"$tmp_vcf\" \
      -Oz -o \"$vcf_gz\"
    tabix -f -p vcf \"$vcf_gz\"
    rm -f \"$tmp_vcf\" \"$tmp_vcf.tbi\"
  "

  echo "[OK] Wrote $vcf_gz and index (bcftools path with rsIDs)"
  exit 0
fi

echo '[ERROR] Missing required tools inside the container (need qctool or bcftools, plus bgzip+tabix).' >&2
exit 2
