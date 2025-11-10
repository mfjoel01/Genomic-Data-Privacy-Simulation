# File: 4_FederatedLearning/R Code/4c_fed_pgs.R
## Federated PGS construction + evaluation (AUC for T2D, R^2 for BMI)
## PATCH: Guard NA indices before subsetting genotype columns.

source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(zoo)
})
source(file.path("4_FederatedLearning", "R Code", "4a_fed_utils.R"))

d4 <- stage_dirs("4_FederatedLearning")
dir.create(d4$out_t, showWarnings = FALSE, recursive = TRUE)

# Inputs
fed_gwas <- readRDS(file.path(DIR_GLOBAL, "fed_gwas_chr22_full.rds"))
ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
pgs_dat  <- readRDS(file.path(DIR_GLOBAL, "pgs_dat_chr22.rds"))

RS    <- fed_gwas$SNP
p_thr <- 5e-8

# Select SNPs for PGS (GWAS genome-wide significant)
sel_t2d <- which(is.finite(fed_gwas$P_T2D_FED) & fed_gwas$P_T2D_FED < p_thr)
sel_bmi <- which(is.finite(fed_gwas$P_BMI_FED) & fed_gwas$P_BMI_FED < p_thr)

# Map to G columns
rs2col  <- match(RS, ann$column_rsids)
idx_t2d <- rs2col[sel_t2d]; w_t2d <- fed_gwas$BETA_T2D_FED[sel_t2d]
idx_bmi <- rs2col[sel_bmi]; w_bmi <- fed_gwas$BETA_BMI_FED[sel_bmi]

# PATCH: drop NAs in column indices to avoid 'subscript contains NAs'
hit <- which(!is.na(idx_t2d)); idx_t2d <- idx_t2d[hit]; w_t2d <- w_t2d[hit]
hit <- which(!is.na(idx_bmi)); idx_bmi <- idx_bmi[hit]; w_bmi <- w_bmi[hit]

# Federated split (use same deterministic split as GWAS)
K    <- as.integer(Sys.getenv("FED_SITES", "5"))
K    <- max(2L, min(10L, K))
seed <- as.integer(Sys.getenv("FED_SPLIT_SEED", "1"))
sites <- make_sites_from_gwas(pgs_dat, K, seed=seed)

# Compute local PGS per site
scores_t2d <- vector("list", K)
scores_bmi <- vector("list", K)
for (k in seq_len(K)) {
  Gk <- sites[[k]]$G
  scores_t2d[[k]] <- if (length(idx_t2d)) as.vector(Gk[, idx_t2d, drop=FALSE] %*% w_t2d) else rep(0, nrow(Gk))
  scores_bmi[[k]] <- if (length(idx_bmi)) as.vector(Gk[, idx_bmi, drop=FALSE] %*% w_bmi) else rep(0, nrow(Gk))
}

# ---- Federated ROC for T2D ----
probs <- seq(0, 1, length.out = 61)
weights <- sapply(sites, function(s) length(s$yT2D))
thr <- global_thresholds_from_sites(scores_t2d, probs = probs, weights = weights)
roc <- aggregate_roc_counts(scores_t2d, lapply(sites, `[[`, "yT2D"), thresholds = thr)
auc <- auc_trap(roc$FPR, roc$TPR); auc <- max(0, min(1, auc))

# ---- Federated OLS for BMI ~ PGS + age + sex ----
stats_list <- lapply(seq_len(K), function(k) {
  X <- cbind(1, scores_bmi[[k]], sites[[k]]$cov$age, sites[[k]]$cov$sex)
  y <- sites[[k]]$yBMI
  ols_cross_stats(X, y)
})
sumstats <- sum_ols_stats(stats_list)
r2res <- r2_from_ols_sums(sumstats)

# ---- Per-site diagnostics (sizes + AUCs for T2D) ----------------------------
roc_one <- function(scores, y) {
  y <- as.integer(y)
  pos <- sum(y == 1); neg <- sum(y == 0)
  if (pos == 0 || neg == 0 || !any(is.finite(scores))) return(NA_real_)
  thr_local <- sort(unique(as.numeric(quantile(scores, probs = seq(0,1,length.out=61), na.rm=TRUE))))
  if (length(thr_local) < 2L) thr_local <- unique(c(min(scores, na.rm = TRUE), max(scores, na.rm = TRUE)))
  TPR <- FPR <- numeric(length(thr_local))
  for (i in seq_along(thr_local)) {
    t <- thr_local[i]
    TPR[i] <- sum(y==1 & scores>=t)/pos
    FPR[i] <- sum(y==0 & scores>=t)/neg
  }
  o <- order(FPR); F <- FPR[o]; T <- TPR[o]
  if (length(F) < 2L) return(NA_real_)
  sum(diff(F) * zoo::rollmean(T, 2))
}
site_sizes <- vapply(sites, function(s) length(s$yT2D), integer(1))
site_aucs  <- mapply(roc_one, scores_t2d, lapply(sites, `[[`, "yT2D"))
site_aucs  <- as.numeric(site_aucs)

# Persist summary
saveRDS(list(
  thresholds = thr, roc = roc, auc_t2d = auc,
  pgs_bmi_beta = r2res$beta, R2_bmi = r2res$R2, R2_adj_bmi = r2res$R2_adj,
  site_sizes = site_sizes, site_auc_t2d = site_aucs
), file.path(DIR_GLOBAL, "fed_pgs_summary_chr22.rds"))

# Write metrics text (top‑level only)
lines <- c(
  sprintf("Sites: %d", K),
  sprintf("Selected SNPs (T2D): %d", length(idx_t2d)),
  sprintf("Selected SNPs (BMI): %d", length(idx_bmi)),
  sprintf("AUC (T2D): %.3f", auc),
  sprintf("Adj. R^2 (BMI): %.3f", r2res$R2_adj)
)
writeLines(lines, file.path(d4$out_t, "pgs_metrics.txt"))
