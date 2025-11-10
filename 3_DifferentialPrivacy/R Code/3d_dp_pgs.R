## DP PGS diagnostics, chr22
## Emits:
##   - dp_roc_T2D.png                    (DP only, Baseline style)
##   - dp_roc_T2D_overlay.png            (DP over Baseline on same plot)
## Also writes dp_summary_chr22.rds and a text metrics file.

source(file.path("R", "config.R"))
source(file.path("3_DifferentialPrivacy", "R Code", "3a_dp_helpers.R"))
suppressPackageStartupMessages({
  library(zoo)
})

if (!exists("png_open")) {
  png_open <- function(filename, width, height) {
    png(filename, width = width, height = height, units = "px",
        type = getOption("bitmapType", "cairo"))
  }
}

# ---------- Load inputs ----------
budgets <- readRDS(file.path(DIR_GLOBAL, "dp_budgets.rds"))
base    <- readRDS(file.path(DIR_GLOBAL, "baseline_pgs_chr22.rds"))
pgs_dat <- readRDS(file.path(DIR_GLOBAL, "pgs_dat_chr22.rds"))
d       <- stage_dirs("3_DifferentialPrivacy")

y_bin   <- as.integer(pgs_dat$yT2D)
scores  <- as.numeric(base$score_t2d)

# ---------- Helper: non private ROC and AUC ----------
make_plain_roc <- function(score, y01) {
  # y01 must be 0 or 1
  ok <- is.finite(score) & (y01 %in% c(0L,1L))
  score <- score[ok]; y01 <- y01[ok]
  if (length(score) < 2L || length(unique(y01)) < 2L) {
    return(list(FPR = numeric(0), TPR = numeric(0), AUC = NA_real_))
  }
  # Sort by decreasing score
  ord <- order(score, decreasing = TRUE)
  s <- score[ord]; y <- y01[ord]
  P <- sum(y == 1L); N <- sum(y == 0L)
  # Threshold at each unique score, compute cumulative TP and FP
  # Use run-length encoding over sorted scores to handle ties efficiently
  r <- rle(s)
  ends <- cumsum(r$lengths)
  y_cum <- cumsum(y)
  tp_at_end <- y_cum[ends]
  fp_at_end <- ends - tp_at_end

  TPR <- tp_at_end / P
  FPR <- fp_at_end / N

  # Add origin
  FPR <- c(0, FPR)
  TPR <- c(0, TPR)

  # Trapezoid AUC
  if (length(FPR) >= 2L) {
    dF <- diff(FPR)
    Tmid <- zoo::rollmean(TPR, 2)
    auc <- pmin(1, pmax(0, sum(dF * Tmid)))
  } else {
    auc <- NA_real_
  }
  list(FPR = FPR, TPR = TPR, AUC = auc)
}

# ---------- DP ROC from privatized histograms ----------
dp_roc <- make_dp_roc(scores, y_bin, eps_total = budgets$eps_plot_roc, n_pts = 60)
okroc  <- is.finite(dp_roc$FPR) & is.finite(dp_roc$TPR)
F_dp   <- dp_roc$FPR[okroc]
TPR_dp <- dp_roc$TPR[okroc]
if (length(F_dp) >= 2L) {
  ord   <- order(F_dp, TPR_dp)
  F_dp  <- F_dp[ord]
  TPR_dp<- TPR_dp[ord]
  dF    <- diff(F_dp)
  Tmid  <- zoo::rollmean(TPR_dp, 2)
  dp_auc <- pmin(1, pmax(0, sum(dF * Tmid)))
} else {
  dp_auc <- NA_real_
}

# ---------- Baseline ROC (non private) ----------
bl_roc <- make_plain_roc(scores, y_bin)
F_bl   <- bl_roc$FPR
TPR_bl <- bl_roc$TPR
bl_auc <- bl_roc$AUC

# ---------- Save DP only ROC, matching Baseline style ----------
png_open(file.path(d$out_g, "dp_roc_T2D.png"), 700, 600)
if (length(F_dp) >= 2L) {
  plot(F_dp, TPR_dp, type = "l", lwd = 2,
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = sprintf("ROC: T2D (AUROC = %s)",
                      ifelse(is.finite(dp_auc), sprintf("%.3f", dp_auc), "NA")))
  abline(0, 1, lty = 2)
} else {
  plot.new()
  title(main = "ROC: T2D (insufficient DP points)")
}
dev.off()

# ---------- Save overlay: DP over Baseline on same axes ----------
png_open(file.path(d$out_g, "dp_roc_T2D_overlay.png"), 700, 600)
if (length(F_dp) >= 2L && length(F_bl) >= 2L) {
  # Draw Baseline first to underlay
  plot(F_bl, TPR_bl, type = "l", lwd = 2,
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = sprintf("ROC: T2D, DP vs Baseline (DP=%.3f, BL=%.3f)",
                      ifelse(is.finite(dp_auc), dp_auc, NA_real_),
                      ifelse(is.finite(bl_auc), bl_auc, NA_real_)))
  lines(F_dp, TPR_dp, lwd = 2, col = 2)
  abline(0, 1, lty = 2)
  legend("bottomright",
         legend = c("Baseline", "DP"),
         lwd = 2, col = c(1, 2), bty = "n")
} else if (length(F_bl) >= 2L) {
  plot(F_bl, TPR_bl, type = "l", lwd = 2,
       xlab = "False Positive Rate", ylab = "True Positive Rate",
       main = sprintf("ROC: T2D Baseline only (BL=%.3f)", ifelse(is.finite(bl_auc), bl_auc, NA_real_)))
  abline(0, 1, lty = 2)
} else {
  plot.new(); title(main = "ROC: T2D (no curves available)")
}
dev.off()

# ---------- Binned BMI vs PGS (unchanged) ----------
bins_target <- 10L
YB <- 50
brks_raw <- quantile(base$score_bmi, probs = seq(0, 1, length.out = bins_target + 1), na.rm = TRUE)
edges <- unique(as.numeric(brks_raw))
if (length(edges) < 2L) {
  rng <- range(base$score_bmi, na.rm = TRUE)
  edges <- unique(pretty(rng, n = bins_target + 1))
}
bin_id <- cut(base$score_bmi, edges, include.lowest = TRUE, labels = FALSE, right = TRUE)
bins <- (length(edges) - 1L)

bin_n   <- tabulate(bin_id, nbins = bins)
bin_sum <- as.numeric(tapply(pgs_dat$yBMI, bin_id, function(v) sum(clip(v, -YB, YB))))
bin_sum[is.na(bin_sum)] <- 0

eps_cnt <- budgets$eps_plot_bins * 0.4
eps_sum <- budgets$eps_plot_bins * 0.6

bin_n_dp   <- pmax(1, round(bin_n + rLaplace(bins, 1/eps_cnt)))
bin_sum_dp <- bin_sum + rLaplace(bins, YB/eps_sum)
kappa <- 5
bin_mu_dp  <- bin_sum_dp / (bin_n_dp + kappa)

x_rank <- seq_len(bins)

png_open(file.path(d$out_g, "dp_binned_bmi_vs_pgs.png"), 700, 600)
if (any(is.finite(bin_mu_dp))) {
  yr <- range(bin_mu_dp[is.finite(bin_mu_dp)])
  if (!all(is.finite(yr))) yr <- c(0,1)
  plot(x_rank, bin_mu_dp, pch = 19, cex = 0.8,
       xlab = "Quantile bin (rank)", ylab = "DP mean BMI per bin",
       main = sprintf("DP BMI vs PGS (epsilon = %.2f, %d bins, YB = %d)", budgets$eps_plot_bins, bins, YB),
       ylim = yr, xlim = range(x_rank))
} else {
  plot.new(); title(main = "DP BMI vs PGS (no finite DP bin means)")
}
dev.off()

# ---------- Persist summary ----------
saveRDS(list(
  dp_auc = dp_auc,
  bl_auc = bl_auc,
  roc = list(
    DP = list(FPR = F_dp, TPR = TPR_dp),
    Baseline = list(FPR = F_bl, TPR = TPR_bl)
  )),
  file.path(DIR_GLOBAL, "dp_summary_chr22.rds")
)

# ---------- Text metrics ----------
lines_out <- c(
  sprintf("DP AUC (T2D): %s", ifelse(is.finite(dp_auc), sprintf("%.3f", dp_auc), "NA")),
  sprintf("Baseline AUC (T2D): %s", ifelse(is.finite(bl_auc), sprintf("%.3f", bl_auc), "NA"))
)
writeLines(lines_out, file.path(d$out_t, "dp_pgs_metrics.txt"))
