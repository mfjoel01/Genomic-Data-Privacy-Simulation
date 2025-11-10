source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(pROC))

base <- readRDS(file.path(DIR_GLOBAL, "baseline_pgs_chr22.rds"))
pgs_dat <- readRDS(file.path(DIR_GLOBAL, "pgs_dat_chr22.rds"))

auroc <- base$auroc_t2d
score_bmi <- base$score_bmi
r2_adj <- base$r2_adj_bmi
roc_obj <- base$roc_t2d

d <- stage_dirs("2_Baseline")

png_open(file.path(d$out_g, "roc_T2D_baseline.png"), 700, 600)
plot(roc_obj, main = sprintf("ROC – T2D (AUROC = %.3f)", auroc), lwd = 2)
abline(0, 1, lty = 2)
dev.off()

png_open(file.path(d$out_g, "bmi_vs_pgs_scatter.png"), 700, 600)
plot(score_bmi, pgs_dat$yBMI, pch = 20, cex = 0.6,
     xlab = "Polygenic Score (BMI – PGS004994)", ylab = "BMI",
     main = sprintf("BMI vs PGS (adj. R^2 = %.3f)", r2_adj))
abline(lm(pgs_dat$yBMI ~ score_bmi), col = "red", lwd = 2)
dev.off()
