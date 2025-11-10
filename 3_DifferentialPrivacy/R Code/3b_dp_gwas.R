## 3b_dp_gwas.R
## DP GWAS on chr22 using private select-then-estimate
## - Private top-K selection per trait
## - Spend remaining budget across the K selected SNPs only
## - zCDP or strong composition, same env knobs as before plus DP_TOPK/DP_SEL_FRAC

source(file.path("R", "config.R"))
source(file.path("3_DifferentialPrivacy", "R Code", "3a_dp_helpers.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(future)
  library(future.apply)
  library(parallelly)
  library(matrixStats)
})

size_gb <- suppressWarnings(as.numeric(Sys.getenv("FUTURE_GLOBALS_MAXSIZE_GB", "16")))
options(future.globals.maxSize = if (is.finite(size_gb) && size_gb > 0) size_gb * 1024^3 else +Inf)
plan(sequential)

# ---- Inputs ----
gwas_dat <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))

G           <- gwas_dat$G
storage.mode(G) <- "integer"
RS          <- ann$column_rsids
CHR         <- as.numeric(ann$CHR)
BP          <- ann$BP

# keep polymorphic
keep <- which(matrixStats::colVars(as.matrix(G), na.rm = TRUE, useNames = FALSE) > 0)
G   <- G[, keep, drop = FALSE]
RS  <- RS[keep]; CHR <- CHR[keep]; BP <- BP[keep]
m   <- ncol(G)

# covariates and outcomes
age_z <- as.numeric(scale(gwas_dat$cov$age))
sex01 <- as.integer(gwas_dat$cov$sex)
y_t2d <- as.integer(gwas_dat$yT2D)
y_bmi <- as.numeric(gwas_dat$yBMI)
rm(gwas_dat); gc()

# ---- Budgets and accounting ----
ACCOUNTING  <- tolower(Sys.getenv("DP_ACCOUNTING", "zcdp"))  # "zcdp" or "strong"
delta_total <- as.numeric(Sys.getenv("DP_DELTA_TOTAL","1e-6"))

# split epsilon across tasks (as in your pipeline)
eps_gwas_t2d <- as.numeric(Sys.getenv("DP_EPS_GWAS_T2D", "2.5"))
eps_gwas_bmi <- as.numeric(Sys.getenv("DP_EPS_GWAS_BMI", "2.5"))

# new knobs
K_default <- as.integer(Sys.getenv("DP_TOPK", "1000"))
K_t2d <- as.integer(Sys.getenv("DP_TOPK_T2D", as.character(K_default)))
K_bmi <- as.integer(Sys.getenv("DP_TOPK_BMI", as.character(K_default)))

sel_frac_default <- as.numeric(Sys.getenv("DP_SEL_FRAC", "0.2"))
sel_frac_t2d <- as.numeric(Sys.getenv("DP_SEL_FRAC_T2D", as.character(sel_frac_default)))
sel_frac_bmi <- as.numeric(Sys.getenv("DP_SEL_FRAC_BMI", as.character(sel_frac_default)))

# selection eps per trait
eps_sel_t2d <- max(1e-6, eps_gwas_t2d * sel_frac_t2d)
eps_sel_bmi <- max(1e-6, eps_gwas_bmi * sel_frac_bmi)

# estimation budget left per trait
eps_est_t2d <- max(1e-6, eps_gwas_t2d - eps_sel_t2d)
eps_est_bmi <- max(1e-6, eps_gwas_bmi - eps_sel_bmi)

# convert estimation budgets to per-query (across K, not m)
if (ACCOUNTING == "zcdp") {
  # reserve half of delta for each trait, same spirit as before
  rho_t2d_tot <- ac_per_query_zcdp(eps_total = eps_est_t2d, delta_total = delta_total * 0.5, m = 1)$rho_total
  rho_bmi_tot <- ac_per_query_zcdp(eps_total = eps_est_bmi, delta_total = delta_total * 0.5, m = 1)$rho_total
  rho_per_t2d <- rho_t2d_tot / max(1L, K_t2d)
  rho_per_bmi <- rho_bmi_tot / max(1L, K_bmi)
  eps_per_t2d <- eps_per_bmi <- NA_real_
  delta_per_t2d <- delta_per_bmi <- NA_real_
} else {
  ac_t2d <- ac_per_query(eps_total = eps_est_t2d, m = max(1L, K_t2d), delta_total = delta_total * 0.5)
  ac_bmi <- ac_per_query(eps_total = eps_est_bmi, m = max(1L, K_bmi), delta_total = delta_total * 0.5)
  eps_per_t2d   <- ac_t2d$eps;  delta_per_t2d <- ac_t2d$delta
  eps_per_bmi   <- ac_bmi$eps;  delta_per_bmi <- ac_bmi$delta
  rho_per_t2d <- rho_per_bmi <- NA_real_
}

# ---- Private selection: top-K by bounded marginal |sum((g-1)*clip(y))| ----
set.seed(20251103)
idx_t2d <- dp_topk_by_product(G, y_t2d, K = K_t2d, eps_sel = eps_sel_t2d, YB_sel = 1)
idx_bmi <- dp_topk_by_product(G, y_bmi, K = K_bmi, eps_sel = eps_sel_bmi, YB_sel = 50)

# ---- DP p-values on selected SNPs only ----
p_t2d_dp <- rep(NA_real_, m)
p_bmi_dp <- rep(NA_real_, m)

# T2D: Analyze-Gauss on LPM (same design matrix as before), YB=1
if (length(idx_t2d)) {
  for (i in idx_t2d) {
    g_z <- as.numeric(scale(G[, i]))
    X   <- cbind(1, g_z, age_z, sex01)
    fit <- if (ACCOUNTING == "zcdp") {
      dp_ols_one(X, as.numeric(y_t2d), rho = rho_per_t2d, lambda = 10.0, L = 3.0, YB = 1.0)
    } else {
      dp_ols_one(X, as.numeric(y_t2d), eps = eps_per_t2d, delta = delta_per_t2d, lambda = 10.0, L = 3.0, YB = 1.0)
    }
    p_t2d_dp[i] <- fit$p[2]
  }
}

# BMI: DP-OLS with YB=50
if (length(idx_bmi)) {
  for (i in idx_bmi) {
    g_z <- as.numeric(scale(G[, i]))
    X   <- cbind(1, g_z, age_z, sex01)
    fit <- if (ACCOUNTING == "zcdp") {
      dp_ols_one(X, y_bmi, rho = rho_per_bmi, lambda = 10.0, L = 3.0, YB = 50)
    } else {
      dp_ols_one(X, y_bmi, eps = eps_per_bmi, delta = delta_per_bmi, lambda = 10.0, L = 3.0, YB = 50)
    }
    p_bmi_dp[i] <- fit$p[2]
  }
}

# ---- Persist DP GWAS table ----
dp_gwas <- data.table(SNP = RS, CHR = CHR, BP = BP,
                      P_T2D_DP = as.numeric(p_t2d_dp),
                      P_BMI_DP = as.numeric(p_bmi_dp))

fout <- file.path(DIR_GLOBAL, "dp_gwas_chr22_full")
data.table::fwrite(dp_gwas, paste0(fout, ".tsv.gz"), sep = "\t")
saveRDS(dp_gwas, paste0(fout, ".rds"))
cat("Saved DP GWAS with private top-K: ", paste0(fout, ".{rds,tsv.gz}"), "\n")

# Small run summary for reproducibility
d <- stage_dirs("3_DifferentialPrivacy")
lines <- c(
  sprintf("ACCOUNTING: %s", ACCOUNTING),
  sprintf("m (polymorphic): %d", m),
  sprintf("K_t2d: %d (eps_sel=%.4f; est=%s)", length(idx_t2d), eps_sel_t2d,
          if (ACCOUNTING == "zcdp") sprintf("rho_per=%.3e", rho_per_t2d) else sprintf("eps_per=%.4f, delta_per=%.1e", eps_per_t2d, delta_per_t2d)),
  sprintf("K_bmi: %d (eps_sel=%.4f; est=%s)", length(idx_bmi), eps_sel_bmi,
          if (ACCOUNTING == "zcdp") sprintf("rho_per=%.3e", rho_per_bmi) else sprintf("eps_per=%.4f, delta_per=%.1e", eps_per_bmi, delta_per_bmi))
)
writeLines(lines, file.path(d$out_t, "dp_gwas_selection_summary.txt"))
