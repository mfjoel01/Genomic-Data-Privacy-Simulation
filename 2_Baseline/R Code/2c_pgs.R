source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(pROC))

pgs_dat  <- readRDS(file.path(DIR_GLOBAL, "pgs_dat_chr22.rds"))
pgs_eff  <- readRDS(file.path(DIR_GLOBAL, "pgs_effects_chr22.rds"))

score_t2d <- as.vector(pgs_dat$G[, pgs_eff$t2d$causal_idx, drop = FALSE] %*% pgs_eff$t2d$effect_sz)
score_bmi <- as.vector(pgs_dat$G[, pgs_eff$bmi$causal_idx, drop = FALSE] %*% pgs_eff$bmi$effect_sz)

roc_obj <- pROC::roc(pgs_dat$yT2D, score_t2d, quiet = TRUE)
auroc   <- as.numeric(roc_obj$auc)

lm_out  <- lm(pgs_dat$yBMI ~ score_bmi + pgs_dat$cov$age + pgs_dat$cov$sex)
r2_adj  <- summary(lm_out)$adj.r.squared

baseline_pgs <- list(score_t2d = score_t2d, score_bmi = score_bmi,
                     auroc_t2d = auroc, r2_adj_bmi = r2_adj, roc_t2d = roc_obj)
saveRDS(baseline_pgs, file.path(DIR_GLOBAL, "baseline_pgs_chr22.rds"))

d <- stage_dirs("2_Baseline")
writeLines(c(sprintf("AUROC (T2D): %.3f", auroc),
             sprintf("Adj. R^2 (BMI): %.3f", r2_adj)),
           file.path(d$out_t, "pgs_metrics.txt"))
cat("✓ Saved baseline_pgs_chr22.rds and metrics\n")
