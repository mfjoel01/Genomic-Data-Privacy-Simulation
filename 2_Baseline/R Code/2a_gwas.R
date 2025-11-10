source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(parallel); library(parallelly)
  library(SNPRelate);  library(gdsfmt)
})

## -- Inputs -------------------------------------------------------------------
gwas_dat <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
pgs_any  <- readRDS(file.path(DIR_GLOBAL, "pgs_any_chr22.rds"))

geno_mat     <- gwas_dat$G
column_rsids <- ann$column_rsids
CHR          <- ann$CHR
BP           <- ann$BP

## -- Keep polymorphic SNPs ----------------------------------------------------
snp_var <- matrixStats::colVars(geno_mat, na.rm = TRUE, useNames = FALSE)
keep    <- which(snp_var > 0)
geno_mat     <- geno_mat[, keep, drop = FALSE]
column_rsids <- column_rsids[keep]
CHR          <- CHR[keep]
BP           <- BP[keep]
n_snps       <- length(keep)

## -- Harmonize effect allele (flip to effect allele where PGS provides one) ---
cache_file <- file.path(DIR_GLOBAL, "alleles_cache.rds")
alleles_vec <- if (file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL

if (is.null(alleles_vec) || length(alleles_vec) != length(ann$snp_keep_idx)) {
  paths <- readRDS(file.path(DIR_GLOBAL, "paths_gds.rds"))
  res <- local({
    gf <- SNPRelate::snpgdsOpen(paths$gds_file, readonly = TRUE)
    on.exit(SNPRelate::snpgdsClose(gf), add = TRUE)
    full <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.allele"))
    list(alleles = full[ann$snp_keep_idx])
  })
  alleles_vec <- res$alleles
  saveRDS(alleles_vec, cache_file)
  cat(sprintf("Built allele cache: %s (n=%d)\n", cache_file, length(alleles_vec)))
} else {
  cat(sprintf("Loaded allele cache: %s (n=%d)\n", cache_file, length(alleles_vec)))
}
alleles_vec <- alleles_vec[keep]

alts <- data.table::tstrsplit(alleles_vec, "/")[[2]]
effect_allele_vec <- alts
pgs_match <- match(column_rsids, pgs_any$rsID)
has_pgs   <- !is.na(pgs_match)
effect_allele_vec[has_pgs] <- pgs_any$effect_allele[pgs_match[has_pgs]]
flip_mask <- effect_allele_vec != alts
if (any(flip_mask)) geno_mat[, flip_mask] <- 2 - geno_mat[, flip_mask]

## -- Covariates / outcomes ----------------------------------------------------
age   <- gwas_dat$cov$age
sex   <- gwas_dat$cov$sex
y_t2d <- gwas_dat$yT2D
y_bmi <- gwas_dat$yBMI

## -- Parallel setup -----------------------------------------------------------
mc_cores <- as.integer(Sys.getenv("BASELINE_CORES",
                                  parallelly::availableCores(omit = 1)))
use_fork <- parallelly::supportsMulticore() && .Platform$OS.type == "unix" && mc_cores > 1L

## -- Per-SNP worker (uses data= to avoid scoping issues) ----------------------
worker_fun <- function(i) {
  g <- geno_mat[, i]
  if (!is.numeric(g) || var(g, na.rm = TRUE) == 0) {
    return(list(FAILED = TRUE, REASON = "constant", SNP = column_rsids[i], CHR = CHR[i], BP = BP[i]))
  }
  df <- data.frame(y_bmi = y_bmi, y_t2d = y_t2d, g = g, age = age, sex = sex)

  out <- tryCatch({
    lm_bmi   <- lm(y_bmi ~ g + age + sex, data = df)
    sb       <- summary(lm_bmi)$coefficients
    beta_bmi <- unname(sb["g", "Estimate"])
    p_bmi    <- unname(sb["g", "Pr(>|t|)"])

    glm_t2d  <- glm(y_t2d ~ g + age + sex, data = df, family = binomial())
    sg       <- summary(glm_t2d)$coefficients
    beta_t2d <- unname(sg["g", "Estimate"])
    p_t2d    <- unname(sg["g", "Pr(>|z|)"])

    list(FAILED = FALSE, SNP = column_rsids[i], CHR = CHR[i], BP = BP[i],
         BETA_BMI = beta_bmi, BETA_T2D = beta_t2d, P_BMI = p_bmi, P_T2D = p_t2d)
  }, error = function(e) list(FAILED = TRUE, REASON = "model_error", SNP = column_rsids[i], CHR = CHR[i], BP = BP[i]))
  out
}

gwas_list <- if (use_fork) {
  parallel::mclapply(seq_len(n_snps), worker_fun, mc.cores = mc_cores)
} else {
  lapply(seq_len(n_snps), worker_fun)
}

all_dt <- data.table::rbindlist(gwas_list, use.names = TRUE, fill = TRUE)
gwas_full <- all_dt[FAILED == FALSE]

## Ensure the key columns are present even if 0 rows (defensive)
needed <- c("SNP","CHR","BP","BETA_BMI","BETA_T2D","P_BMI","P_T2D")
for (k in setdiff(needed, names(gwas_full))) gwas_full[, (k) := as.numeric(NA_real_)]

## Drop helper columns and save
cols_to_drop <- intersect(c("FAILED","REASON"), names(gwas_full))
if (length(cols_to_drop)) gwas_full[, (cols_to_drop) := NULL]

data.table::fwrite(gwas_full, file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.tsv.gz"), sep = "\t")
saveRDS(gwas_full, file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds"))

## Summary text
d <- stage_dirs("2_Baseline")
writeLines(sprintf("Attempted: %d | Succeeded: %d | Failed: %d",
                   n_snps, nrow(gwas_full), sum(all_dt$FAILED %in% TRUE, na.rm = TRUE)),
           file.path(d$out_t, "gwas_run_summary.txt"))
cat("✓ Saved baseline GWAS in Data/global and summary text\n")
