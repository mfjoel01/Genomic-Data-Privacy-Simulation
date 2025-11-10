#!/usr/bin/env Rscript
## 3e_epsilon_sweep.R
## Epsilon sweep for DP GWAS utility (chr22) using private select-then-estimate

source(file.path("R", "config.R"))
source(file.path("3_DifferentialPrivacy", "R Code", "3a_dp_helpers.R"))
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(zoo)
})

if (!exists("png_open")) {
  png_open <- function(filename, width, height) {
    png(filename, width = width, height = height, units = "px",
        type = getOption("bitmapType", "cairo"))
  }
}

d3 <- stage_dirs("3_DifferentialPrivacy")
dir.create(d3$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d3$out_t, showWarnings = FALSE, recursive = TRUE)

# ---- Inputs ----
ann        <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
gwas_dat   <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
base_gwas  <- readRDS(file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds"))

G  <- gwas_dat$G
RS <- ann$column_rsids
CHR<- as.numeric(ann$CHR)
BP <- ann$BP

keep <- which(matrixStats::colVars(as.matrix(G), na.rm = TRUE, useNames = FALSE) > 0)
G  <- G[, keep, drop = FALSE]
RS <- RS[keep]; CHR <- CHR[keep]; BP <- BP[keep]
m  <- ncol(G)

y_t2d <- as.integer(gwas_dat$yT2D)
y_bmi <- as.numeric(gwas_dat$yBMI)
age_z <- as.numeric(scale(gwas_dat$cov$age))
sex01 <- as.integer(gwas_dat$cov$sex)

p_thr <- as.numeric(Sys.getenv("P_THR_CHR22", "5e-8"))

# grid and knobs
eps_grid <- as.numeric(strsplit(Sys.getenv("EPS_GRID", "0.5,1,1.5,2,3,4"), ",")[[1]])
ACCOUNTING <- tolower(Sys.getenv("DP_ACCOUNTING", "zcdp"))
delta_total <- as.numeric(Sys.getenv("DP_DELTA_TOTAL","1e-6"))
K_default <- as.integer(Sys.getenv("DP_TOPK", "1000"))
K_t2d <- as.integer(Sys.getenv("DP_TOPK_T2D", as.character(K_default)))
K_bmi <- as.integer(Sys.getenv("DP_TOPK_BMI", as.character(K_default)))
sel_frac_default <- as.numeric(Sys.getenv("DP_SEL_FRAC", "0.2"))
sel_frac_t2d <- as.numeric(Sys.getenv("DP_SEL_FRAC_T2D", as.character(sel_frac_default)))
sel_frac_bmi <- as.numeric(Sys.getenv("DP_SEL_FRAC_BMI", as.character(sel_frac_default)))

# baseline hits for recall
base_t2d <- base_gwas[match(RS, base_gwas$SNP), ]
base_bmi <- base_gwas[match(RS, base_gwas$SNP), ]
true_hits_t2d <- which(is.finite(base_t2d$P_T2D) & base_t2d$P_T2D < p_thr)
true_hits_bmi <- which(is.finite(base_bmi$P_BMI) & base_bmi$P_BMI < p_thr)
n_true_t2d <- length(true_hits_t2d); n_true_bmi <- length(true_hits_bmi)

eval_one <- function(eps_trait) {
  # split trait epsilon into selection and estimation
  eps_sel_t2d <- max(1e-6, eps_trait * sel_frac_t2d)
  eps_sel_bmi <- max(1e-6, eps_trait * sel_frac_bmi)
  eps_est_t2d <- max(1e-6, eps_trait - eps_sel_t2d)
  eps_est_bmi <- max(1e-6, eps_trait - eps_sel_bmi)

  if (ACCOUNTING == "zcdp") {
    rho_t2d_tot <- ac_per_query_zcdp(eps_total = eps_est_t2d, delta_total = delta_total * 0.5, m = 1)$rho_total
    rho_bmi_tot <- ac_per_query_zcdp(eps_total = eps_est_bmi, delta_total = delta_total * 0.5, m = 1)$rho_total
    rho_per_t2d <- rho_t2d_tot / max(1L, K_t2d)
    rho_per_bmi <- rho_bmi_tot / max(1L, K_bmi)
    eps_per_t2d <- eps_per_bmi <- NA_real_; delta_per_t2d <- delta_per_bmi <- NA_real_
  } else {
    ac_t2d <- ac_per_query(eps_total = eps_est_t2d, m = max(1L, K_t2d), delta_total = delta_total * 0.5)
    ac_bmi <- ac_per_query(eps_total = eps_est_bmi, m = max(1L, K_bmi), delta_total = delta_total * 0.5)
    eps_per_t2d   <- ac_t2d$eps;  delta_per_t2d <- ac_t2d$delta
    eps_per_bmi   <- ac_bmi$eps;  delta_per_bmi <- ac_bmi$delta
    rho_per_t2d <- rho_per_bmi <- NA_real_
  }

  # select
  set.seed(20251103)
  idx_t2d <- dp_topk_by_product(G, y_t2d, K = K_t2d, eps_sel = eps_sel_t2d, YB_sel = 1)
  idx_bmi <- dp_topk_by_product(G, y_bmi, K = K_bmi, eps_sel = eps_sel_bmi, YB_sel = 50)

  # fit on selected
  p_t2d <- rep(NA_real_, m); p_bmi <- rep(NA_real_, m)
  for (i in idx_t2d) {
    g_z <- as.numeric(scale(G[, i])); X <- cbind(1, g_z, age_z, sex01)
    fit <- if (ACCOUNTING == "zcdp") dp_ols_one(X, as.numeric(y_t2d), rho = rho_per_t2d, lambda = 10, L = 3, YB = 1)
           else dp_ols_one(X, as.numeric(y_t2d), eps = eps_per_t2d, delta = delta_per_t2d, lambda = 10, L = 3, YB = 1)
    p_t2d[i] <- fit$p[2]
  }
  for (i in idx_bmi) {
    g_z <- as.numeric(scale(G[, i])); X <- cbind(1, g_z, age_z, sex01)
    fit <- if (ACCOUNTING == "zcdp") dp_ols_one(X, y_bmi, rho = rho_per_bmi, lambda = 10, L = 3, YB = 50)
           else dp_ols_one(X, y_bmi, eps = eps_per_bmi, delta = delta_per_bmi, lambda = 10, L = 3, YB = 50)
    p_bmi[i] <- fit$p[2]
  }

  # recall at p_thr, relative to baseline true hits
  hits_t2d <- which(is.finite(p_t2d) & p_t2d < p_thr)
  hits_bmi <- which(is.finite(p_bmi) & p_bmi < p_thr)
  rec_t2d  <- if (n_true_t2d) sum(hits_t2d %in% true_hits_t2d) / n_true_t2d else NA_real_
  rec_bmi  <- if (n_true_bmi) sum(hits_bmi %in% true_hits_bmi) / n_true_bmi else NA_real_

  # fidelity on overlapping finite p-values
  l10 <- function(p) -log10(pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - 1e-15))
  o_t <- which(is.finite(p_t2d) & is.finite(base_t2d$P_T2D))
  o_b <- which(is.finite(p_bmi) & is.finite(base_bmi$P_BMI))
  rho_logp_t2d <- suppressWarnings(cor(l10(p_t2d[o_t]), l10(base_t2d$P_T2D[o_t]), method = "spearman"))
  rho_logp_bmi <- suppressWarnings(cor(l10(p_bmi[o_b]), l10(base_bmi$P_BMI[o_b]), method = "spearman"))

  data.frame(
    eps_total = eps_trait,
    K_t2d = length(idx_t2d), K_bmi = length(idx_bmi),
    recall_t2d = rec_t2d, recall_bmi = rec_bmi,
    spearman_rho_t2d = rho_logp_t2d, spearman_rho_bmi = rho_logp_bmi
  )
}

rows <- lapply(eps_grid, eval_one)
sweep_dt <- data.table::rbindlist(rows, fill = TRUE)

# Persist and plot
f_csv <- file.path(d3$out_t, "epsilon_sweep_metrics.csv")
f_rds <- file.path(DIR_GLOBAL, "epsilon_sweep_metrics.rds")
data.table::fwrite(sweep_dt, f_csv)
saveRDS(sweep_dt, f_rds)

png_open(file.path(d3$out_g, "epsilon_sweep_recall.png"), 900, 700)
plot(sweep_dt$eps_total, sweep_dt$recall_t2d, type = "o", pch = 16, lwd = 2,
     xlab = expression(epsilon~"(per trait)"), ylab = "Recall @ P threshold",
     main = sprintf("Privacy-Utility (select then estimate), P<%s", format(p_thr, digits=2, scientific=TRUE)),
     ylim = c(0,1))
lines(sweep_dt$eps_total, sweep_dt$recall_bmi, type = "o", pch = 17, lwd = 2)
legend("bottomright", inset = 0.02, legend = c("T2D", "BMI"), pch = c(16,17), lwd = 2, bty = "n")
dev.off()

png_open(file.path(d3$out_g, "epsilon_sweep_fidelity.png"), 900, 700)
plot(sweep_dt$eps_total, sweep_dt$spearman_rho_t2d, type = "o", pch = 16, lwd = 2,
     xlab = expression(epsilon~"(per trait)"), ylab = expression(Spearman~rho~"(-log"[10]*" p)"),
     main = "Fidelity to baseline vs epsilon")
lines(sweep_dt$eps_total, sweep_dt$spearman_rho_bmi, type = "o", pch = 17, lwd = 2)
abline(h = 1, lty = 3)
legend("bottomright", inset = 0.02, legend = c("T2D", "BMI"), pch = c(16,17), lwd = 2, bty = "n")
dev.off()

cat("Epsilon sweep complete.\n")
