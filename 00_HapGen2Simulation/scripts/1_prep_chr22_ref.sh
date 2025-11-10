#!/bin/bash
set -euo pipefail

# ---------- config ----------
container="${container:-/storage/singularity/mixtures.sif}"
VCF_IN="input/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
PANEL="input/1000G/integrated_call_samples_v3.20130502.ALL.panel"

REF_DIR="ref"
EUR_VCF="${REF_DIR}/chr22_EUR_only.vcf.gz"
HAP_PREFIX="${REF_DIR}/chr22_EUR_ref"   # -> .hap/.legend(.gz) + .samples
MAP_FILE="input/1000G/genetic_map_chr22_combined_b37.txt"
# ----------------------------

echo "[prep] container: $container"
mkdir -p "$REF_DIR"

# sanity-check map presence early
if [[ ! -s "${MAP_FILE}" ]]; then
  echo "[prep][ERROR] Map not found: ${MAP_FILE}"
  echo "              Put genetic_map_chr22_combined_b37.txt in input/1000G/"
  exit 1
fi
echo "[prep] found map: ${MAP_FILE}"

# 1) make EUR sample list (robust to super_pop vs super_population)
echo "[prep] extracting EUR sample IDs from ${PANEL} ..."
awk '
  NR==1 {
    for (i=1;i<=NF;i++){
      h=tolower($i); gsub(/-/, "_", h);
      if (h=="sample") s=i;
      if (h=="super_pop" || h=="super_population") p=i;
    }
    if (!s || !p) { print "[prep][AWK] could not locate sample/super_pop column"; exit 2 }
    next
  }
  $p=="EUR" { print $s }
' "$PANEL" > "${REF_DIR}/eur_samples_ids.txt"

# 2) subset VCF to EUR, biallelic SNPs
echo "[prep] subsetting VCF to EUR bi-allelic SNPs..."
singularity exec "$container" bash -lc "
  bcftools view -S ${REF_DIR}/eur_samples_ids.txt -m2 -M2 -v snps -Oz -o ${EUR_VCF} ${VCF_IN} && \
  bcftools index -t ${EUR_VCF}
"

# 3) convert EUR VCF -> hap/legend/sample for HAPGEN2
echo "[prep] converting EUR VCF -> hap/legend/sample..."
singularity exec "$container" bash -lc "
  bcftools convert --haplegendsample ${HAP_PREFIX} --vcf-ids ${EUR_VCF}
"

# bcftools emits .hap.gz and .legend.gz and .samples (plural). Decompress & alias.
if [[ -s "${HAP_PREFIX}.hap.gz" ]]; then
  gunzip -f "${HAP_PREFIX}.hap.gz"
fi
if [[ -s "${HAP_PREFIX}.legend.gz" ]]; then
  gunzip -f "${HAP_PREFIX}.legend.gz"
fi
# optional: ensure a .sample alias exists alongside .samples
if [[ -s "${HAP_PREFIX}.samples" && ! -e "${HAP_PREFIX}.sample" ]]; then
  ln -sf "$(basename ${HAP_PREFIX}).samples" "${HAP_PREFIX}.sample"
fi

# 4) pick a polymorphic dummy disease locus (for fallback use)
#    We scan hap rows + legend positions and pick the first SNP with 0 < alt_count < max.
LEG="${HAP_PREFIX}.legend"
HAP="${HAP_PREFIX}.hap"
echo "[prep] selecting a polymorphic dummy locus..."
paste <(tail -n +2 "${LEG}" | awk '{print $2}') "${HAP}" | \
awk '{
  pos=$1; sum=0;
  for(i=2;i<=NF;i++) sum += $i;
  if (sum>0 && sum < (NF-1)) { print pos; exit }
}' > "${REF_DIR}/dummy_pos.txt"

if [[ ! -s "${REF_DIR}/dummy_pos.txt" ]]; then
  echo "[prep][WARN] Could not find polymorphic site for dummy locus; using first position from legend."
  awk 'NR==2{print $2; exit}' "${LEG}" > "${REF_DIR}/dummy_pos.txt"
fi
echo "[prep] dummy locus: $(cat ${REF_DIR}/dummy_pos.txt)"

echo "[prep] done. created:"
ls -lh "${HAP_PREFIX}.hap" "${HAP_PREFIX}.legend" "${MAP_FILE}" "${REF_DIR}/dummy_pos.txt"
