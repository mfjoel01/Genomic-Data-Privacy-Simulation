# 2_Baseline/R Code/2b_manhattan_qq.R
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qqman))

gwas_dt  <- readRDS(file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds"))
pgs_t2d  <- readRDS(file.path(DIR_GLOBAL, "pgs_t2d_chr22.rds"))
pgs_bmi  <- readRDS(file.path(DIR_GLOBAL, "pgs_bmi_chr22.rds"))

setDT(gwas_dt)

# ---- normalize column names (tolerant to variants) ----
normalize_cols <- function(dt) {
  pick_and_rename <- function(target, candidates) {
    ln <- tolower(names(dt))
    idx <- which(ln %in% tolower(candidates))[1]
    if (!is.na(idx)) data.table::setnames(dt, old = names(dt)[idx], new = target)
  }
  pick_and_rename("SNP",   c("snp","rsid","rsID","rs_id"))
  pick_and_rename("CHR",   c("chr","chrom","chromosome","CHR"))
  pick_and_rename("BP",    c("bp","pos","position","BP"))
  pick_and_rename("P_T2D", c("P_T2D","p_t2d","P.t2d","pt2d","p.t2d"))
  pick_and_rename("P_BMI", c("P_BMI","p_bmi","P.bmi","pbmi","p.bmi"))
  required <- c("SNP","CHR","BP","P_T2D","P_BMI")
  missing  <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(sprintf("Missing required columns after normalization: %s. Have: %s",
                 paste(missing, collapse = ", "),
                 paste(names(dt), collapse = ", ")))
  }
}
normalize_cols(gwas_dt)

# ---- constants ----
base_col <- "#2F8CA6"
p_thr    <- 5e-8
log_thr  <- -log10(p_thr)

# Ensure numeric CHR for qqman
gwas_dt[, CHR := as.numeric(CHR)]

# Highlight sets from PGS libraries
hi_t2d <- intersect(unique(pgs_t2d$rsID), gwas_dt$SNP)
hi_bmi <- intersect(unique(pgs_bmi$rsID), gwas_dt$SNP)

d <- stage_dirs("2_Baseline")

# ---------- diagnostics: lambdaGC and hit counts ----------
clip01 <- function(p) pmin(pmax(p, .Machine$double.eps), 1 - 1e-15)
chisq_BMI <- qchisq(1 - clip01(gwas_dt$P_BMI), df = 1)
chisq_T2D <- qchisq(1 - clip01(gwas_dt$P_T2D), df = 1)
lambda_BMI <- median(chisq_BMI, na.rm = TRUE) / qchisq(0.5, 1)
lambda_T2D <- median(chisq_T2D, na.rm = TRUE) / qchisq(0.5, 1)
diag_txt <- sprintf(
  "lambdaGC (BMI)=%.3f | lambdaGC (T2D)=%.3f | BMI hits (P<5e-8)=%d | T2D hits (P<5e-8)=%d",
  lambda_BMI, lambda_T2D,
  sum(gwas_dt$P_BMI < 5e-8, na.rm = TRUE),
  sum(gwas_dt$P_T2D < 5e-8, na.rm = TRUE)
)
writeLines(diag_txt, file.path(d$out_t, "qq_diagnostics.txt"))

# ---------- Manhattan plots (single color, no left gap) ----------
# T2D
png_open(file.path(d$out_g, "manhattan_T2D_baseline.png"), 1200, 600)
qqman::manhattan(
  gwas_dt[, .(SNP, CHR, BP, P = P_T2D)],
  main = "GWAS, T2D (baseline, PGS003443)",
  genomewideline = log_thr,
  suggestiveline = FALSE,
  #highlight = hi_t2d,
  col = base_col,
  xlim = range(gwas_dt$BP, na.rm = TRUE),
  xaxs = "i"
)
dev.off()

# BMI
png_open(file.path(d$out_g, "manhattan_BMI_baseline.png"), 1200, 600)
qqman::manhattan(
  gwas_dt[, .(SNP, CHR, BP, P = P_BMI)],
  main = "GWAS, BMI (baseline, PGS004994)",
  genomewideline = log_thr,
  suggestiveline = FALSE,
  #highlight = hi_bmi,
  col = base_col,
  xlim = range(gwas_dt$BP, na.rm = TRUE),
  xaxs = "i"
)
dev.off()

# ---------- QQ plots ----------
png_open(file.path(d$out_g, "qqplots_baseline.png"), 1200, 600)
par(mfrow = c(1, 2))
qq(gwas_dt$P_T2D, main = sprintf("QQ, T2D (lambda = %.2f)", lambda_T2D))
qq(gwas_dt$P_BMI, main = sprintf("QQ, BMI (lambda = %.2f)", lambda_BMI))
dev.off()

writeLines(sprintf("Genome-wide p-threshold (fixed): %.3e | -log10: %.2f", p_thr, log_thr),
           file.path(d$out_t, "bonferroni_threshold.txt"))
