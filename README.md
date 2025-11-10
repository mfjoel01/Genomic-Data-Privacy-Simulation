# Genomic Data Privacy Simulation

End‑to‑end simulation for **secure GWAS & PGS benchmarking** on chr22.  
Compares **Baseline**, **Differential Privacy**, and **Federated Learning** pipelines.  
**No data is hosted in this repository** — scripts fetch/use public resources.

---

## What’s here
- Reproducible R + Bash/Slurm code for HPC environments.
- Stages: `00_HapGen2Simulation/` → `0_Prep/` → `1_Phenotype/` → `2_Baseline/` → `3_DifferentialPrivacy/` → `4_FederatedLearning/` → `5_MethodComparison/`.
- Slurm drivers in `Scripts/`.

## What you need (not in repo)

> Put the files below under `00_HapGen2Simulation/input/1000G/` unless noted.

**1) 1000 Genomes (IGSR) Phase 3, GRCh37 (chr22)**
- VCF: `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz`  
  • HTTPS mirror (UCSC):  
  https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  
  • Index (`.tbi`):  
  https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi  
  • IGSR Phase 3 overview: https://www.internationalgenome.org/data-portal/data-collection/phase-3

- Sample panel (superpopulation labels), **required to subset EUR**:  
  `integrated_call_samples_v3.20130502.ALL.panel`  
  HTTPS mirror: https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel

- **Genetic map (GRCh37)** for HAPGEN2:  
  `genetic_map_chr22_combined_b37.txt`  
  Get from the IMPUTE2/1000G Phase 3 archive (contains `genetic_map_chr*_combined_b37.txt`):  
  https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html  
  *(Download the tarball(s) from that page and extract the chr22 file.)*

**2) PGS Catalog scoring files (harmonized, GRCh37)**  
Download the harmonized scoring files (`*_hmPOS_GRCh37.txt.gz`) from the score pages:
- **T2D**: PGS003443 → “Download Score” → Harmonized → `PGS003443_hmPOS_GRCh37.txt.gz`  
  https://www.pgscatalog.org/score/PGS003443/
- **BMI**: PGS004994 → “Download Score” → Harmonized → `PGS004994_hmPOS_GRCh37.txt.gz`  
  https://www.pgscatalog.org/score/PGS004994/  
PGS Catalog download notes (hmPOS GRCh37/38): https://www.pgscatalog.org/downloads/

> Optional tools (installed in your container/environment): **bcftools**, **HAPGEN2**, **R** + packages used in the scripts.

---

## Quick start (HPC / Slurm)

```bash
# 0) Clone
git clone https://github.com/mfjoel01/Genomic-Data-Privacy-Simulation.git
cd Genomic-Data-Privacy-Simulation

# 1) Create expected folders for HAPGEN2 stage
bash 00_HapGen2Simulation/scripts/0_make_tree.sh

# 2) Place required files:
#    - .../input/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#    - .../input/1000G/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
#    - .../input/1000G/integrated_call_samples_v3.20130502.ALL.panel
#    - .../input/1000G/genetic_map_chr22_combined_b37.txt
#    - Data/pgs/PGS003443_hmPOS_GRCh37.txt.gz
#    - Data/pgs/PGS004994_hmPOS_GRCh37.txt.gz

# 3) (Optional) point to your Singularity/Apptainer image when submitting
#    export container=/path/to/your_image.sif

# 4) Run the stages with Slurm drivers (examples)
sbatch Scripts/1_pheno_driver.sh
sbatch Scripts/2_baseline_driver.sh
sbatch Scripts/3_dp_driver.sh
sbatch Scripts/4_federated_driver.sh
sbatch Scripts/5_method_compare_driver.sh
# or chain stages with:
bash Scripts/run_stages.sh 1 5
````

**Outputs:** each stage writes plots and text summaries under its `output/` directory.

---

## Notes

* Built for **HPC clusters** using **Slurm**.
* Uses IGSR/1000G Phase 3 on **GRCh37** (hg19).
* PGS weights come from the **PGS Catalog** and are applied on chr22 for the simulated cohort.
