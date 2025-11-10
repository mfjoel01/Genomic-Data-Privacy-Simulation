## Stage 3 plots, chr22 only

source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(qqman))
suppressPackageStartupMessages(library(data.table))

if (!exists("png_open")) {
  png_open <- function(filename, width, height) {
    png(filename, width = width, height = height, units = "px",
        type = getOption("bitmapType", "cairo"))
  }
}

dp_gwas <- readRDS(file.path(DIR_GLOBAL, "dp_gwas_chr22_full.rds"))
pgs_t2d <- readRDS(file.path(DIR_GLOBAL, "pgs_t2d_chr22.rds"))
pgs_bmi <- readRDS(file.path(DIR_GLOBAL, "pgs_bmi_chr22.rds"))

base_col <- "#cee7eb"
p_thr   <- 5e-8
log_thr <- -log10(p_thr)

setDT(dp_gwas)
dp_gwas[, CHR := as.numeric(CHR)]

clipP <- function(p) { p <- as.numeric(p); pmin(pmax(p, .Machine$double.eps), 1 - 1e-15) }

pick_highlights <- function(dt, pgs_rsids, max_n = 1000L) {
  hi <- intersect(unique(pgs_rsids), dt$SNP)
  if (!length(hi)) return(character(0))
  pos <- match(hi, dt$SNP)
  ord <- order(dt$P[pos], na.last = NA)
  hi[ord[seq_len(min(length(ord), max_n))]]
}

d <- stage_dirs("3_DifferentialPrivacy")

# T2D
dt_t2d <- dp_gwas[is.finite(P_T2D_DP)][, .(SNP, CHR, BP, P = clipP(P_T2D_DP))]
hi_t2d <- pick_highlights(dt_t2d, pgs_t2d$rsID, max_n = 1000L)

png_open(file.path(d$out_g, "manhattan_T2D_DP.png"), 1200, 600)
if (nrow(dt_t2d) > 0) {
  qqman::manhattan(
    dt_t2d,
    main = "DP GWAS, T2D (Analyze‑Gauss LPM)",
    genomewideline = log_thr,
    suggestiveline = FALSE,
    highlight = hi_t2d,
    col = base_col,
    xlim = range(dt_t2d$BP, na.rm = TRUE),
    xaxs = "i"
  )
} else { plot.new(); title(main = "DP GWAS, T2D, no finite p-values") }
dev.off()

# BMI
dt_bmi <- dp_gwas[is.finite(P_BMI_DP)][, .(SNP, CHR, BP, P = clipP(P_BMI_DP))]
hi_bmi <- pick_highlights(dt_bmi, pgs_bmi$rsID, max_n = 1000L)

png_open(file.path(d$out_g, "manhattan_BMI_DP.png"), 1200, 600)
if (nrow(dt_bmi) > 0) {
  qqman::manhattan(
    dt_bmi,
    main = "DP GWAS, BMI (Analyze‑Gauss OLS)",
    genomewideline = log_thr,
    suggestiveline = FALSE,
    highlight = hi_bmi,
    col = base_col,
    xlim = range(dt_bmi$BP, na.rm = TRUE),
    xaxs = "i"
  )
} else { plot.new(); title(main = "DP GWAS, BMI, no finite p-values") }
dev.off()

# QQ
png_open(file.path(d$out_g, "qqplots_DP.png"), 1200, 600)
par(mfrow = c(1, 2))
if (nrow(dt_t2d) > 0 && any(is.finite(dt_t2d$P))) {
  qq(dt_t2d$P, main = "DP QQ, T2D")
} else {
  plot.new(); title(main = "DP QQ, T2D, no finite p")
}
if (nrow(dt_bmi) > 0 && any(is.finite(dt_bmi$P))) {
  qq(dt_bmi$P, main = "DP QQ, BMI")
} else {
  plot.new(); title(main = "DP QQ, BMI, no finite p")
}
dev.off()

writeLines(sprintf("Genome‑wide p threshold (fixed): %.3e | -log10: %.2f", p_thr, log_thr),
           file.path(d$out_t, "dp_bonferroni_threshold.txt"))
