source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(data.table))

# Load stage-0 outputs
ann <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))

# PGS directory + files (default to PROJECT_ROOT/Data/pgs; allow env overrides)
PGS_DIR      <- Sys.getenv("PGS_DIR",      file.path(PROJECT_ROOT, "Data", "pgs"))
PGS_T2D_FILE <- Sys.getenv("PGS_T2D_FILE", file.path(PGS_DIR, "PGS003443_hmPOS_GRCh37.txt.gz"))
PGS_BMI_FILE <- Sys.getenv("PGS_BMI_FILE", file.path(PGS_DIR, "PGS004994_hmPOS_GRCh37.txt.gz"))
#PGS_BMI_FILE <- Sys.getenv("PGS_BMI_FILE", file.path(PGS_DIR, "PGS003980_hmPOS_GRCh37.txt.gz"))

# Sanity checks
if (!file.exists(PGS_T2D_FILE) || !file.exists(PGS_BMI_FILE)) {
  cat("Looking for PGS files under:", PGS_DIR, "\n")
  if (dir.exists(PGS_DIR)) print(list.files(PGS_DIR))
  stop("Missing PGS files. Expected:\n  ", PGS_T2D_FILE, "\n  ", PGS_BMI_FILE)
}

pgs_t2d_hm <- fread(PGS_T2D_FILE, showProgress = FALSE)
pgs_bmi_hm <- fread(PGS_BMI_FILE, showProgress = FALSE)

# Keep chr22 only (accept "22" or 22)
pgs_t2d_chr22_all <- pgs_t2d_hm[chr_name %in% c("22", 22)]
pgs_bmi_chr22_all <- pgs_bmi_hm[chr_name %in% c("22", 22)]

# Overlap with retained SNPs
rs_keep <- ann$column_rsids
pgs_t2d_chr22 <- pgs_t2d_chr22_all[rsID %in% rs_keep]
pgs_bmi_chr22 <- pgs_bmi_chr22_all[rsID %in% rs_keep]

# Union effect allele table for harmonisation later
pgs_any_chr22 <- unique(rbindlist(list(
  pgs_t2d_chr22[, .(rsID, effect_allele, effect_weight)],
  pgs_bmi_chr22[, .(rsID, effect_allele, effect_weight)]
)))

# Save for later stages
saveRDS(pgs_t2d_chr22, file.path(DIR_GLOBAL, "pgs_t2d_chr22.rds"))
saveRDS(pgs_bmi_chr22, file.path(DIR_GLOBAL, "pgs_bmi_chr22.rds"))
saveRDS(pgs_any_chr22, file.path(DIR_GLOBAL, "pgs_any_chr22.rds"))

# Simple counts report
d <- stage_dirs("1_Phenotype")
lines <- c(
  sprintf("PGS003443 chr22 all: %d", nrow(pgs_t2d_chr22_all)),
  sprintf("PGS003443 chr22 overlap: %d", nrow(pgs_t2d_chr22)),
  sprintf("PGS004994 chr22 all: %d", nrow(pgs_bmi_chr22_all)),
  sprintf("PGS004994 chr22 overlap: %d", nrow(pgs_bmi_chr22))
)
dir.create(d$out_t, showWarnings = FALSE, recursive = TRUE)
writeLines(lines, file.path(d$out_t, "pgs_overlap_counts.txt"))
