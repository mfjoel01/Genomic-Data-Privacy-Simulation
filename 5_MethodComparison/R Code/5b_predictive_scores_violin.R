# File: 5_MethodComparison/R Code/5b_predictive_scores_violin.R
# Build per‑person linear predictors (scores) using method‑specific GWAS betas
# and plot violin distributions (Baseline / DP / Federated) for BMI and T2D.
# PATCH: NA‑safe weights + densities to prevent density() crash on NA values.

source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(matrixStats)
})

d5 <- stage_dirs("5_MethodComparison")
safe_png <- function(path, w=900, h=800) png_open(path, w, h)

# -------- helpers --------
align_env <- local({
  ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
  gwas_dat <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
  G        <- gwas_dat$G
  keep     <- which(matrixStats::colVars(G, na.rm = TRUE, useNames = FALSE) > 0)
  list(G=G[, keep, drop=FALSE], rs = ann$column_rsids[keep], keep=keep)
})

# Map rs_w/beta_w onto rs_G. Treat non‑finite betas as 0 so scores aren't NA.
align_weights <- function(rs_G, rs_w, beta_w) {
  v <- numeric(length(rs_G))
  m <- match(rs_G, rs_w)
  if (!is.null(beta_w)) {
    bw <- as.numeric(beta_w)
    bw[!is.finite(bw)] <- 0
    hit <- which(!is.na(m))
    if (length(hit)) v[hit] <- bw[m[hit]]
  }
  v
}

# Score is G %*% w with non‑finite w treated as 0
scores_from <- function(w) {
  w <- as.numeric(w); w[!is.finite(w)] <- 0
  drop(align_env$G %*% w)
}

# NA‑safe violin drawer (filters non‑finite values; handles constants)
draw_violins <- function(vecs, title, outfile) {
  nm <- names(vecs); k <- length(vecs)
  dens <- lapply(vecs, function(v) {
    vv <- as.numeric(v); vv <- vv[is.finite(vv)]
    if (length(vv) > 1 && stats::sd(vv) > 0) density(vv) else list(x = 0, y = 0)
  })
  y_max <- max(vapply(dens, function(d) max(d$y, na.rm = TRUE), 0.0))
  x_rng <- range(unlist(lapply(dens, `[[`, "x")), na.rm = TRUE)

  safe_png(outfile, 1000, 750)
  par(mar=c(5,6,4,2)+0.1)
  plot(0,0,type="n", xlim=c(0.5, k+0.5), ylim=x_rng,
       xaxt="n", xlab="", ylab="Score", main=title)
  for (i in seq_len(k)) {
    d <- dens[[i]]
    if (!is.list(d) || any(!is.finite(d$y))) next
    y <- d$x
    w <- if (y_max > 0) d$y / y_max * 0.4 else d$y * 0
    polygon(c(i - w, rev(i + w)), c(y, rev(y)), border = NA, col = adjustcolor("gray60", 0.6))
    vi <- as.numeric(vecs[[i]]); vi <- vi[is.finite(vi)]
    if (length(vi)) {
      points(rep(i, length(vi)), vi, pch = 16, cex = 0.3, col = adjustcolor("black", 0.15))
      segments(i - 0.35, mean(vi), i + 0.35, mean(vi), lwd = 3)
    }
  }
  axis(1, at=seq_len(k), labels=nm, las=2)
  dev.off()
}

# -------- load betas / compute DP ones if needed --------
base <- readRDS(file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds"))
fed  <- readRDS(file.path(DIR_GLOBAL, "fed_gwas_chr22_full.rds"))
dpb  <- if (file.exists(file.path(DIR_GLOBAL, "dp_bmi_betas_chr22.rds"))) {
          readRDS(file.path(DIR_GLOBAL, "dp_bmi_betas_chr22.rds"))
        } else {
          # Generate via 5a if missing
          source(file.path("5_MethodComparison", "R Code", "5a_compare_gwas_effects.R"), local = TRUE)
          readRDS(file.path(DIR_GLOBAL, "dp_bmi_betas_chr22.rds"))
        }

# DP logistic per‑allele β via noisy aggregated counts (3x2)
dp_t2d_betas_file <- file.path(DIR_GLOBAL, "dp_t2d_betas_chr22.rds")
get_dp_t2d_betas <- function() {
  if (file.exists(dp_t2d_betas_file)) return(readRDS(dp_t2d_betas_file))

  bud <- if (file.exists(file.path(DIR_GLOBAL, "dp_budgets.rds"))) {
    readRDS(file.path(DIR_GLOBAL, "dp_budgets.rds"))
  } else list(eps_gwas_t2d = 2.5)

  gwas_dat <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
  ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
  G <- gwas_dat$G; y <- as.integer(gwas_dat$yT2D); rs <- ann$column_rsids
  keep <- which(matrixStats::colVars(G, na.rm = TRUE, useNames = FALSE) > 0)
  G <- G[, keep, drop=FALSE]; rs <- rs[keep]
  eps <- max(1e-6, bud$eps_gwas_t2d / ncol(G))

  # Laplace helper
  source(file.path("3_DifferentialPrivacy", "R Code", "3a_dp_helpers.R"))

  # Build noisy 3x2 tables and fit aggregated‑binomial logistic: glm(cbind(case,control) ~ g)
  set.seed(20251011)
  beta <- numeric(ncol(G))
  for (i in seq_len(ncol(G))) {
    g_fac <- factor(G[, i], levels=c(0,1,2))
    tab <- table(g_fac, y)
    base <- matrix(as.numeric(tab), 3, 2)

    # add Laplace noise per cell; clamp and stabilize with tiny pseudocount
    noisy <- base + matrix(rLaplace(6, 1/eps), 3, 2)
    noisy[!is.finite(noisy) | noisy < 0] <- 0
    noisy <- noisy + 1e-6

    # aggregated logistic with 3 rows (g=0,1,2); drop empty rows after noise
    df <- data.frame(g = 0:2, case = noisy[,2], ctrl = noisy[,1])
    df <- df[(df$case + df$ctrl) > 0, , drop=FALSE]
    if (nrow(df) < 2) { beta[i] <- NA_real_; next }

    fit <- try(suppressWarnings(glm(cbind(case, ctrl) ~ g, data=df, family=binomial())), silent=TRUE)
    beta[i] <- if (inherits(fit, "try-error") || !("g" %in% names(coef(fit)))) NA_real_ else unname(coef(fit)["g"])
  }
  out <- list(rs=rs, beta=beta)
  saveRDS(out, dp_t2d_betas_file)
  out
}
dpt <- get_dp_t2d_betas()

# -------- build weights aligned to G (keep subset) --------
rsG <- align_env$rs

w_base_bmi <- align_weights(rsG, base$SNP, base$BETA_BMI)
w_base_t2d <- align_weights(rsG, base$SNP, base$BETA_T2D)

w_fed_bmi  <- align_weights(rsG, fed$SNP,  fed$BETA_BMI_FED)
w_fed_t2d  <- align_weights(rsG, fed$SNP,  fed$BETA_T2D_FED)

w_dp_bmi   <- align_weights(rsG, dpb$rs,   dpb$beta)
w_dp_t2d   <- align_weights(rsG, dpt$rs,   dpt$beta)

# -------- scores --------
s_base_bmi <- scores_from(w_base_bmi)
s_dp_bmi   <- scores_from(w_dp_bmi)
s_fed_bmi  <- scores_from(w_fed_bmi)

s_base_t2d <- scores_from(w_base_t2d)
s_dp_t2d   <- scores_from(w_dp_t2d)
s_fed_t2d  <- scores_from(w_fed_t2d)

# -------- violins --------
f1 <- file.path(d5$out_g, "violin_scores_BMI.png")
draw_violins(list(Baseline=s_base_bmi, DP=s_dp_bmi, Federated=s_fed_bmi),
             "Predictive score distributions — BMI", f1)

f2 <- file.path(d5$out_g, "violin_scores_T2D.png")
draw_violins(list(Baseline=s_base_t2d, DP=s_dp_t2d, Federated=s_fed_t2d),
             "Predictive score distributions — T2D", f2)

writeLines(c(basename(f1), basename(f2)),
           file.path(d5$out_t, "5b_outputs.txt"))
