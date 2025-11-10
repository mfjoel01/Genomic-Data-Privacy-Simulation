## 4b_fed_gwas.R — Federated GWAS (parallel, column-safe)
source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(matrixStats)
  library(future); library(future.apply); library(parallelly)
})
source(file.path("4_FederatedLearning", "R Code", "4a_fed_utils.R"))

d4 <- stage_dirs("4_FederatedLearning")
dir.create(d4$out_t, showWarnings = FALSE, recursive = TRUE)

# ---- Inputs ------------------------------------------------------------------
gwas_dat <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))

G_full       <- gwas_dat$G
column_rsids <- ann$column_rsids
CHR          <- as.numeric(ann$CHR)
BP           <- ann$BP

# ---- NEW (CRITICAL): Allele harmonization to Baseline convention -------------
# Flip genotype columns so that the counted allele equals the PGS effect allele
# where available; otherwise keep ALT. This mirrors 2a_gwas.R behavior.
pgs_any <- readRDS(file.path(DIR_GLOBAL, "pgs_any_chr22.rds"))
cache_file <- file.path(DIR_GLOBAL, "alleles_cache.rds")

# Reuse cached alleles from Stage 2; rebuild from GDS if absent/mismatched
alleles_vec <- if (file.exists(cache_file)) tryCatch(readRDS(cache_file), error = function(e) NULL) else NULL
if (is.null(alleles_vec) || length(alleles_vec) != length(ann$snp_keep_idx)) {
  paths <- readRDS(file.path(DIR_GLOBAL, "paths_gds.rds"))
  res <- local({
    gf <- SNPRelate::snpgdsOpen(paths$gds_file, readonly = TRUE)
    on.exit(SNPRelate::snpgdsClose(gf), add = TRUE)
    full <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gf, "snp.allele"))
    full[ann$snp_keep_idx]
  })
  alleles_vec <- res
  saveRDS(alleles_vec, cache_file)
}

# ALT allele for each column (Stage-0 kept set → matches G_full columns)
alts <- data.table::tstrsplit(alleles_vec, "/")[[2]]
effect_allele_vec <- alts
pgs_match <- match(column_rsids, pgs_any$rsID)
has_pgs   <- !is.na(pgs_match)
effect_allele_vec[has_pgs] <- pgs_any$effect_allele[pgs_match[has_pgs]]

flip_mask <- effect_allele_vec != alts
if (any(flip_mask, na.rm = TRUE)) G_full[, flip_mask] <- 2 - G_full[, flip_mask]

# Reflect harmonized genotypes in gwas_dat before site construction
gwas_dat$G <- G_full

# ---- Keep polymorphic SNPs (variance is invariant to flipping) ---------------
keep <- which(matrixStats::colVars(as.matrix(G_full), na.rm = TRUE, useNames = FALSE) > 0)
G            <- G_full[, keep, drop = FALSE]
column_rsids <- column_rsids[keep]
CHR          <- CHR[keep]
BP           <- BP[keep]
p <- ncol(G); n <- nrow(G)

# ---- Federated split ----------------------------------------------------------
K <- as.integer(Sys.getenv("FED_SITES", "5")); K <- max(2L, min(10L, K))
seed <- as.integer(Sys.getenv("FED_SPLIT_SEED", "1"))
sites <- make_sites_from_gwas(gwas_dat, K, seed = seed)

# IMPORTANT: align each site's G to the kept columns so per-SNP index i matches
sites <- lapply(sites, function(s) { s$G <- s$G[, keep, drop = FALSE]; s })

# ---- Parallel plan (safe cap; override with FED_CORES) ------------------------
cap <- as.integer(Sys.getenv("FED_CORES", "8"))
nworkers <- max(1L, min(cap, parallelly::availableCores(omit = 1)))
use_fork <- parallelly::supportsMulticore() && .Platform$OS.type == "unix" && nworkers > 1L
if (use_fork) {
  options(parallelly.fork.enable = TRUE)
  plan(multicore, workers = nworkers)
} else {
  plan(sequential)
}

# ---- Per-SNP worker -----------------------------------------------------------
worker_fun <- function(i) {
  b_log <- s_log <- numeric(0)
  b_lin <- s_lin <- numeric(0)

  # Gather site-local stats for SNP i
  for (k in seq_len(K)) {
    gk <- sites[[k]]$G[, i]
    outL <- site_glm_logistic(sites[[k]]$yT2D, gk, sites[[k]]$cov)
    if (isTRUE(outL["ok"] == 1)) { b_log <- c(b_log, outL["beta"]); s_log <- c(s_log, outL["se"]) }
    outC <- site_lm_linear(sites[[k]]$yBMI,  gk, sites[[k]]$cov)
    if (isTRUE(outC["ok"] == 1)) { b_lin <- c(b_lin, outC["beta"]); s_lin <- c(s_lin, outC["se"]) }
  }
  cmbL <- combine_meta(b_log, s_log)
  cmbC <- combine_meta(b_lin, s_lin)

  c(beta_log = cmbL$beta, se_log = cmbL$se, p_log = cmbL$p, Q_log = cmbL$Q, I2_log = cmbL$I2, k_log = cmbL$k,
    beta_lin = cmbC$beta, se_lin = cmbC$se, p_lin = cmbC$p, Q_lin = cmbC$Q, I2_lin = cmbC$I2, k_lin = cmbC$k)
}

# ---- Run (parallel over SNPs) -------------------------------------------------
res <- future_lapply(
  X = seq_len(p), FUN = worker_fun,
  future.seed = 20251011, future.globals = FALSE
)
plan(sequential)  # restore

R <- do.call(rbind, res)

# ---- Collect ----------------------------------------------------------------------------
beta_log <- as.numeric(R[, "beta_log"]); se_log <- as.numeric(R[, "se_log"])
p_log    <- as.numeric(R[, "p_log"]);    Q_log <- as.numeric(R[, "Q_log"])
I2_log   <- as.numeric(R[, "I2_log"]);   k_log <- as.integer(R[, "k_log"])

beta_lin <- as.numeric(R[, "beta_lin"]); se_lin <- as.numeric(R[, "se_lin"])
p_lin    <- as.numeric(R[, "p_lin"]);    Q_lin <- as.numeric(R[, "Q_lin"])
I2_lin   <- as.numeric(R[, "I2_lin"]);   k_lin <- as.integer(R[, "k_lin"])

fed_gwas <- data.table(
  SNP = column_rsids, CHR = CHR, BP = BP,
  BETA_T2D_FED = beta_log, SE_T2D_FED = se_log, P_T2D_FED = p_log, Q_T2D = Q_log, I2_T2D = I2_log, K_T2D = k_log,
  BETA_BMI_FED = beta_lin, SE_BMI_FED = se_lin, P_BMI_FED = p_lin, Q_BMI = Q_lin, I2_BMI = I2_lin, K_BMI = k_lin
)

# ---- Persist & diagnostics ----------------------------------------------------
fout <- file.path(DIR_GLOBAL, "fed_gwas_chr22_full")
data.table::fwrite(fed_gwas, paste0(fout, ".tsv.gz"), sep = "\t")
saveRDS(fed_gwas, paste0(fout, ".rds"))
cat("✓ Saved FED GWAS:", paste0(fout, ".{rds,tsv.gz}"), "\n")

p_thr <- 5e-8
lambda_gc <- function(pvals) {
  pv <- pmin(pmax(as.numeric(pvals), .Machine$double.eps), 1 - 1e-15)
  chisq <- qchisq(1 - pv, df = 1)
  median(chisq[is.finite(chisq)], na.rm=TRUE) / qchisq(0.5, df = 1)
}
lambda_T2D <- lambda_gc(fed_gwas$P_T2D_FED)
lambda_BMI <- lambda_gc(fed_gwas$P_BMI_FED)
diag_txt <- sprintf(
  "Sites: %d | n=%d | p=%d\nlambdaGC (BMI)=%.3f | lambdaGC (T2D)=%.3f | BMI hits (P<5e-8)=%d | T2D hits (P<5e-8)=%d",
  K, n, p, lambda_BMI, lambda_T2D,
  sum(fed_gwas$P_BMI_FED < p_thr, na.rm=TRUE),
  sum(fed_gwas$P_T2D_FED < p_thr, na.rm=TRUE)
)
writeLines(diag_txt, file.path(d4$out_t, "gwas_run_summary.txt"))
writeLines(sprintf("Genome-wide p-threshold (fixed): %.3e | -log10: %.2f", p_thr, -log10(p_thr)),
           file.path(d4$out_t, "bonferroni_threshold.txt"))
