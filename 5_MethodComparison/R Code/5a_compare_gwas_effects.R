# 5a_compare_gwas_effects.R
# Heatmap-like (smooth density) y=x comparisons of SNP betas:
#   (i)  BMI: Baseline β vs DP β (DP‑OLS; rescaled to per‑allele)
#   (ii) T2D: Baseline β vs Federated β (meta-analysis)

source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(matrixStats)
})

# ---------- helpers ----------
safe_png <- function(path, w=900, h=800) png_open(path, w, h)

harmonize_and_get_G <- function() {
  # Return list(G_keep, rs_keep, age_z, sex01, keep_idx, sd_g)
  ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
  gwas_dat <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
  pgs_any  <- readRDS(file.path(DIR_GLOBAL, "pgs_any_chr22.rds"))
  cache_file <- file.path(DIR_GLOBAL, "alleles_cache.rds")
  stopifnot(file.exists(cache_file))
  alleles_vec <- readRDS(cache_file)

  G  <- gwas_dat$G
  rs <- ann$column_rsids

  # keep: polymorphic
  keep <- which(matrixStats::colVars(G, na.rm = TRUE, useNames = FALSE) > 0)
  G <- G[, keep, drop = FALSE]
  rs <- rs[keep]
  alleles_vec <- alleles_vec[keep]

  # flip to "effect allele" preference (same as Stage 2 baseline GWAS)
  alts <- vapply(strsplit(alleles_vec, "/", fixed=TRUE), `[`, "", 2)
  eff  <- alts
  m <- match(rs, pgs_any$rsID); has <- !is.na(m)
  eff[has] <- pgs_any$effect_allele[m[has]]
  flip <- eff != alts
  if (any(flip, na.rm = TRUE)) G[, flip] <- 2 - G[, flip]

  list(
    G = G,
    rs = rs,
    age_z = as.numeric(scale(gwas_dat$cov$age)),
    sex01 = as.integer(gwas_dat$cov$sex),
    keep_idx = keep,
    sd_g = matrixStats::colSds(G, na.rm = TRUE)
  )
}

dp_bmi_betas_file <- file.path(DIR_GLOBAL, "dp_bmi_betas_chr22.rds")

compute_dp_bmi_betas <- function() {
  bud <- if (file.exists(file.path(DIR_GLOBAL, "dp_budgets.rds"))) {
    readRDS(file.path(DIR_GLOBAL, "dp_budgets.rds"))
  } else {
    list(eps_gwas_bmi = 2.5, delta_total = 1e-6)
  }

  src <- harmonize_and_get_G()
  G <- src$G; rs <- src$rs; age_z <- src$age_z; sex01 <- src$sex01; sd_g <- src$sd_g
  n <- nrow(G); p <- ncol(G)

  # DP‑OLS uses standardized columns for stable clipping/noise; then rescale β to per‑allele units.
  X_base <- cbind(1, matrix(NA_real_, n, 1), age_z, sex01)  # placeholder for g_z
  eps_per <- bud$eps_gwas_bmi / p
  delta_per <- bud$delta_total / p

  # y (BMI) from GWAS set
  y <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))$yBMI

  source(file.path("3_DifferentialPrivacy", "R Code", "3a_dp_helpers.R"))
  set.seed(20251011)
  beta_unscaled <- numeric(p)
  for (i in seq_len(p)) {
    g <- as.numeric(scale(G[, i]))
    X_base[, 2] <- g
    fit <- dp_ols_one(X_base, y, eps = eps_per, delta = delta_per, lambda = 1.0, L = 3.0, YB = 80)
    b_z <- fit$beta[2]
    s <- sd_g[i]; if (!is.finite(s) || s <= 0) s <- 1
    beta_unscaled[i] <- b_z / s
  }

  out <- list(rs=rs, beta=beta_unscaled)
  saveRDS(out, dp_bmi_betas_file)
  out
}

get_dp_bmi_betas <- function() {
  if (file.exists(dp_bmi_betas_file)) {
    tryCatch(readRDS(dp_bmi_betas_file), error=function(e) compute_dp_bmi_betas())
  } else compute_dp_bmi_betas()
}

# ---------- (i) BMI: Baseline vs DP ----------
base <- readRDS(file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds"))
dpb  <- get_dp_bmi_betas()

# align
m <- match(dpb$rs, base$SNP)
x <- base$BETA_BMI[m]
y <- dpb$beta
ok <- which(is.finite(x) & is.finite(y))

f <- file.path(stage_dirs("5_MethodComparison")$out_g, "effect_heatmap_BMI_base_vs_DP.png")
safe_png(f, 1000, 900)
par(mar=c(5,5,4,2)+0.1)
smoothScatter(x[ok], y[ok],
              xlab=expression(beta[Baseline]~"(BMI)"),
              ylab=expression(beta[DP]~"(BMI)"),
              main="GWAS SNP effect agreement — BMI (Baseline vs DP)")
abline(0,1,lty=2,col="gray40")
dev.off()

# ---------- (ii) T2D: Baseline vs Federated ----------
fed <- readRDS(file.path(DIR_GLOBAL, "fed_gwas_chr22_full.rds"))
m2 <- match(fed$SNP, base$SNP)
x2 <- base$BETA_T2D[m2]
y2 <- fed$BETA_T2D_FED
ok2 <- which(is.finite(x2) & is.finite(y2))

f2 <- file.path(stage_dirs("5_MethodComparison")$out_g, "effect_heatmap_T2D_base_vs_FED.png")
safe_png(f2, 1000, 900)
par(mar=c(5,5,4,2)+0.1)
smoothScatter(x2[ok2], y2[ok2],
              xlab=expression(beta[Baseline]~"(T2D)"),
              ylab=expression(beta[Federated]~"(T2D)"),
              main="GWAS SNP effect agreement — T2D (Baseline vs Federated)")
abline(0,1,lty=2,col="gray40")
dev.off()

writeLines(c(basename(f), basename(f2)),
           file.path(stage_dirs("5_MethodComparison")$out_t, "5a_outputs.txt"))
