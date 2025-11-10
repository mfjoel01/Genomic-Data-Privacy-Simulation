#!/bin/bash
set -euo pipefail

mkdir -p input/1000G
mkdir -p input/maps
mkdir -p ref
mkdir -p hapgen2
mkdir -p plink
mkdir -p scripts

echo "[ok] created: input/{1000G,maps} ref hapgen2 plink scripts"
echo "[reminder] put these here:"
echo "  input/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
echo "  input/1000G/integrated_call_samples_v3.20130502.ALL.panel"
echo "  input/1000G/genetic_map_chr22_combined_b37.txt   # <— your map file"
