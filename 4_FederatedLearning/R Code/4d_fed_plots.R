## 4_FederatedLearning/R Code/4d_fed_plots.R
## Plots for federated GWAS and PGS (NA safe)
source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(qqman); library(zoo)
})

d4 <- stage_dirs("4_FederatedLearning")
dir.create(d4$out_g, showWarnings = FALSE, recursive = TRUE)

fed_gwas <- readRDS(file.path(DIR_GLOBAL, "fed_gwas_chr22_full.rds"))
pgs_t2d  <- readRDS(file.path(DIR_GLOBAL, "pgs_t2d_chr22.rds"))
pgs_bmi  <- readRDS(file.path(DIR_GLOBAL, "pgs_bmi_chr22.rds"))
pgs_sum  <- readRDS(file.path(DIR_GLOBAL, "fed_pgs_summary_chr22.rds"))

base_col <- "#69B9C6"
p_thr   <- 5e-8
log_thr <- -log10(p_thr)

clipP <- function(p) { p <- as.numeric(p); pmin(pmax(p, .Machine$double.eps), 1 - 1e-15) }

# Ensure CHR numeric for qqman
fed_gwas$CHR <- as.numeric(fed_gwas$CHR)

pgs_rsids_t2d <- unique(pgs_t2d$rsID)
pgs_rsids_bmi <- unique(pgs_bmi$rsID)

# ---------- Manhattan: T2D ----------
png_open(file.path(d4$out_g, "manhattan_T2D_fed.png"), 1200, 600)
dt_t2d <- fed_gwas[is.finite(fed_gwas$CHR) & is.finite(fed_gwas$BP) & is.finite(fed_gwas$P_T2D_FED),
                   .(SNP, CHR, BP, P = clipP(P_T2D_FED))]
if (nrow(dt_t2d)) {
  qqman::manhattan(
    dt_t2d,
    main = "Federated GWAS, T2D",
    genomewideline = log_thr,
    suggestiveline = FALSE,
    #highlight = intersect(pgs_rsids_t2d, dt_t2d$SNP),
    col = base_col,
    xlim = range(dt_t2d$BP, na.rm = TRUE),
    xaxs = "i"
  )
} else { plot.new(); title(main = "Federated GWAS, T2D, no finite p") }
dev.off()

# ---------- Manhattan: BMI ----------
png_open(file.path(d4$out_g, "manhattan_BMI_fed.png"), 1200, 600)
dt_bmi <- fed_gwas[is.finite(fed_gwas$CHR) & is.finite(fed_gwas$BP) & is.finite(fed_gwas$P_BMI_FED),
                   .(SNP, CHR, BP, P = clipP(P_BMI_FED))]
if (nrow(dt_bmi)) {
  qqman::manhattan(
    dt_bmi,
    main = "Federated GWAS, BMI",
    genomewideline = log_thr,
    suggestiveline = FALSE,
    #highlight = intersect(pgs_rsids_bmi, dt_bmi$SNP),
    col = base_col,
    xlim = range(dt_bmi$BP, na.rm = TRUE),
    xaxs = "i"
  )
} else { plot.new(); title(main = "Federated GWAS, BMI, no finite p") }
dev.off()

# ---------- QQ plots ----------
png_open(file.path(d4$out_g, "qqplots_fed.png"), 1200, 600)
par(mfrow = c(1, 2))
pt <- dt_t2d$P; pb <- dt_bmi$P
if (length(pt) > 1) qq(pt, main = "FED QQ, T2D") else { plot.new(); title("FED QQ, T2D, no finite p") }
if (length(pb) > 1) qq(pb, main = "FED QQ, BMI") else { plot.new(); title("FED QQ, BMI, no finite p") }
dev.off()

# ---------- ROC plot (T2D) ----------
png_open(file.path(d4$out_g, "roc_T2D_fed.png"), 700, 600)
plot(pgs_sum$roc$FPR, pgs_sum$roc$TPR, type = "l", lwd = 2,
     xlab = "False Positive Rate", ylab = "True Positive Rate",
     main = sprintf("ROC, T2D (Federated PGS, AUC = %.3f)", pgs_sum$auc_t2d))
abline(0, 1, lty = 2)
dev.off()

# ---------- Binned BMI vs PGS (simple line from aggregated beta) ----------
png_open(file.path(d4$out_g, "bmi_vs_pgs_binned_fed.png"), 700, 600)
beta <- pgs_sum$pgs_bmi_beta
plot(c(-3, 3), c(-3, 3), type = "n",
     xlab = "Polygenic Score, BMI, relative", ylab = "BMI, centered",
     main = sprintf("BMI vs Federated PGS (adj. R^2 = %.3f)", pgs_sum$R2_adj_bmi))
abline(a = beta[1], b = beta[2], lwd = 2, col = "red")
dev.off()

# ================= Additional diagnostics =================

# (A) Per-site sizes and AUC bars (if present)
if (!is.null(pgs_sum$site_sizes)) {
  png_open(file.path(d4$out_g, "site_size_bars.png"), 800, 600)
  barplot(pgs_sum$site_sizes, main = "Federated site sizes", xlab = "Site", ylab = "N")
  dev.off()
}
if (!is.null(pgs_sum$site_auc_t2d)) {
  png_open(file.path(d4$out_g, "site_auc_bars.png"), 800, 600)
  barplot(pgs_sum$site_auc_t2d, main = "Per-site ROC AUC, T2D", xlab = "Site", ylab = "AUC")
  abline(h = pgs_sum$auc_t2d, lty = 2)
  dev.off()
}

# (B) I^2 histograms (heterogeneity)
png_open(file.path(d4$out_g, "i2_histograms.png"), 1200, 600)
par(mfrow = c(1, 2))
hist(fed_gwas$I2_T2D, breaks = 30, main = "I^2, T2D", xlab = "I^2", xlim = c(0, 1))
hist(fed_gwas$I2_BMI, breaks = 30, main = "I^2, BMI",  xlab = "I^2", xlim = c(0, 1))
dev.off()

# (C) Concordance with centralized baseline (if available)
base_path <- file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds")
if (file.exists(base_path)) {
  base_gwas <- readRDS(base_path)
  m <- match(fed_gwas$SNP, base_gwas$SNP)
  clip01 <- function(p) pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - 1e-15)
  L10 <- function(x) -log10(clip01(x))
  png_open(file.path(d4$out_g, "fed_vs_base_logp_scatter.png"), 1200, 600)
  par(mfrow = c(1, 2))
  plot(L10(base_gwas$P_T2D[m]), L10(fed_gwas$P_T2D_FED), pch = 20, cex = .5,
       xlab = "Baseline -log10 p, T2D", ylab = "Federated -log10 p, T2D",
       main = "-log10 p concordance, T2D"); abline(0, 1, lty = 2)
  plot(L10(base_gwas$P_BMI[m]), L10(fed_gwas$P_BMI_FED), pch = 20, cex = .5,
       xlab = "Baseline -log10 p, BMI", ylab = "Federated -log10 p, BMI",
       main = "-log10 p concordance, BMI");  abline(0, 1, lty = 2)
  dev.off()

  png_open(file.path(d4$out_g, "fed_vs_base_beta_bmi_scatter.png"), 700, 600)
  plot(base_gwas$BETA_BMI[m], fed_gwas$BETA_BMI_FED, pch = 20, cex = .5,
       xlab = "Baseline beta, BMI", ylab = "Federated beta, BMI",
       main = "Effect size concordance, BMI"); abline(0, 1, lty = 2)
  dev.off()
}

# (D) Heterogeneity vs significance
clip01 <- function(p) pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - 1e-15)
L10 <- function(x) -log10(clip01(x))
png_open(file.path(d4$out_g, "het_vs_sig_scatter.png"), 1200, 600)
par(mfrow = c(1, 2))
plot(L10(fed_gwas$P_T2D_FED), fed_gwas$I2_T2D, pch = 20, cex = .5,
     xlab = "-log10 p, T2D", ylab = "I^2", main = "I^2 vs significance, T2D")
plot(L10(fed_gwas$P_BMI_FED), fed_gwas$I2_BMI, pch = 20, cex = .5,
     xlab = "-log10 p, BMI",  ylab = "I^2", main = "I^2 vs significance, BMI")
dev.off()

# Threshold line (text)
writeLines(sprintf("Genome-wide p-threshold (fixed): %.3e | -log10: %.2f", p_thr, log_thr),
           file.path(d4$out_t, "bonferroni_threshold.txt"))
