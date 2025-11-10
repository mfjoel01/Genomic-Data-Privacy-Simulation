#!/usr/bin/env Rscript
# 4e_fed_k_sweep.R — Sweep #sites (K) and plot utility vs K (runs from 4_driver.R)
source(file.path("R", "config.R"))
suppressPackageStartupMessages({ library(data.table) })

# ---- helpers ----
clip01 <- function(p) pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - 1e-15)
lambda_gc <- function(pvals) {
  pv <- clip01(pvals)
  chisq <- qchisq(1 - pv, df = 1)
  median(chisq[is.finite(chisq)], na.rm = TRUE) / qchisq(0.5, df = 1)
}
L10 <- function(x) -log10(clip01(x))
parse_grid <- function(var, default = c(2,3,4,5,6,8,10)) {
  s <- trimws(Sys.getenv(var, ""))
  if (!nzchar(s)) return(default)
  as.integer(strsplit(s, ",")[[1]])
}

d4 <- stage_dirs("4_FederatedLearning")
dir.create(d4$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d4$out_t, showWarnings = FALSE, recursive = TRUE)

# Baseline (optional; skip concordance if not present)
base_path <- file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds")
have_base <- file.exists(base_path)
if (have_base) base <- readRDS(base_path)

p_thr   <- 5e-8
K_grid  <- parse_grid("FED_SITES_GRID")   # default: 2,3,4,5,6,8,10
seed    <- as.integer(Sys.getenv("FED_SPLIT_SEED", "1"))

rows <- list()
for (K in K_grid) {
  # Keep split deterministic across re-runs
  Sys.setenv(FED_SITES = as.character(K))
  Sys.setenv(FED_SPLIT_SEED = as.character(seed))

  # Re-run federated GWAS/PGS in this shared R session
  source(file.path("4_FederatedLearning", "R Code", "4b_fed_gwas.R"), local = TRUE)
  source(file.path("4_FederatedLearning", "R Code", "4c_fed_pgs.R"),  local = TRUE)

  fed_gwas <- readRDS(file.path(DIR_GLOBAL, "fed_gwas_chr22_full.rds"))
  fed_pgs  <- readRDS(file.path(DIR_GLOBAL, "fed_pgs_summary_chr22.rds"))

  # Save K-tagged copies for traceability
  saveRDS(fed_gwas, file.path(DIR_GLOBAL, sprintf("fed_gwas_chr22_full_K%d.rds", K)))
  saveRDS(fed_pgs,  file.path(DIR_GLOBAL, sprintf("fed_pgs_summary_chr22_K%d.rds", K)))

  # Metrics
  auc_t2d   <- as.numeric(fed_pgs$auc_t2d)
  r2_bmi    <- as.numeric(fed_pgs$R2_adj_bmi)
  hits_t2d  <- sum(is.finite(fed_gwas$P_T2D_FED) & fed_gwas$P_T2D_FED < p_thr, na.rm = TRUE)
  hits_bmi  <- sum(is.finite(fed_gwas$P_BMI_FED) & fed_gwas$P_BMI_FED < p_thr, na.rm = TRUE)
  lam_t2d   <- lambda_gc(fed_gwas$P_T2D_FED)
  lam_bmi   <- lambda_gc(fed_gwas$P_BMI_FED)

  rho_t2d <- rho_bmi <- NA_real_
  if (have_base) {
    m <- match(fed_gwas$SNP, base$SNP)
    rho_t2d <- suppressWarnings(cor(L10(fed_gwas$P_T2D_FED), L10(base$P_T2D[m]),
                                    method = "spearman", use = "complete.obs"))
    rho_bmi <- suppressWarnings(cor(L10(fed_gwas$P_BMI_FED),  L10(base$P_BMI[m]),
                                    method = "spearman", use = "complete.obs"))
  }

  rows[[length(rows)+1]] <- data.frame(
    K = K,
    auc_t2d = auc_t2d,
    r2_bmi  = r2_bmi,
    hits_t2d = hits_t2d,
    hits_bmi = hits_bmi,
    lambda_t2d = lam_t2d,
    lambda_bmi = lam_bmi,
    rho_logp_t2d = rho_t2d,
    rho_logp_bmi = rho_bmi
  )
}

dt <- data.table::rbindlist(rows)
setorder(dt, K)

# Persist table
data.table::fwrite(dt, file.path(d4$out_t, "fed_k_sweep_metrics.csv"))
saveRDS(dt, file.path(DIR_GLOBAL, "fed_k_sweep_metrics.rds"))

# ---- Plots ----
png_open(file.path(d4$out_g, "k_sweep_auc_r2.png"), 900, 700)
yr <- range(c(dt$auc_t2d, dt$r2_bmi), na.rm = TRUE)
plot(dt$K, dt$auc_t2d, type="o", pch=16, lwd=2, xlab="# federated sites (K)",
     ylab="Utility", ylim=yr, main="Utility vs number of sites")
lines(dt$K, dt$r2_bmi, type="o", pch=17, lwd=2)
legend("bottomright", bty="n", pch=c(16,17), lwd=2, legend=c("AUC (T2D)", "Adj. R² (BMI)"))
dev.off()

png_open(file.path(d4$out_g, "k_sweep_hits.png"), 900, 700)
yr <- range(c(dt$hits_t2d, dt$hits_bmi), na.rm = TRUE)
plot(dt$K, dt$hits_t2d, type="o", pch=16, lwd=2, xlab="# federated sites (K)",
     ylab="# significant loci (P<5e-8)", main="Genome-wide discoveries vs K", ylim = yr)
lines(dt$K, dt$hits_bmi, type="o", pch=17, lwd=2)
legend("topleft", bty="n", pch=c(16,17), lwd=2, legend=c("T2D", "BMI"))
dev.off()

if (have_base) {
  png_open(file.path(d4$out_g, "k_sweep_fidelity.png"), 900, 700)
  plot(dt$K, dt$rho_logp_t2d, type="o", pch=16, lwd=2, ylim=c(0,1),
       xlab="# federated sites (K)", ylab=expression(Spearman~rho~"(-log"[10]*" p)"),
       main="Fidelity to baseline vs K")
  lines(dt$K, dt$rho_logp_bmi, type="o", pch=17, lwd=2)
  legend("bottomright", bty="n", pch=c(16,17), lwd=2, legend=c("T2D", "BMI"))
  dev.off()
}

png_open(file.path(d4$out_g, "k_sweep_lambda.png"), 900, 700)
yr <- range(c(dt$lambda_t2d, dt$lambda_bmi), na.rm = TRUE)
plot(dt$K, dt$lambda_t2d, type="o", pch=16, lwd=2, ylim=yr,
     xlab="# federated sites (K)", ylab=expression(lambda[GC]),
     main=expression("Genomic inflation " * lambda[GC] * " vs K"))
lines(dt$K, dt$lambda_bmi, type="o", pch=17, lwd=2)
abline(h=1, lty=3)
legend("bottomright", bty="n", pch=c(16,17), lwd=2, legend=c("T2D", "BMI"))
dev.off()

cat("✓ K-sweep complete. Saved metrics and plots.\n")
