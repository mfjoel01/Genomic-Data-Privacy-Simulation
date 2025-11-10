# 5e_run_to_run_variability.R
# Run the three pipelines multiple times, save per-run betas, and
# plot run-to-run variability as pirate-style distributions of Δβ.
#
# Artifacts per run are cached under:
#   Data/global/run2run/{BASE,DP,FED}/run_##.rds
#
# Env knobs:
#   MC_REPS             default 10      number of runs per method
#   PGS_BETA_SUBSAMPLE  default 200     number of SNPs per trait to evaluate
#   RUN2RUN_SEED        default 20251030 base random seed
#   FED_SITES           default 5       sites for Federated
#   FED_SPLIT_SEED      default 1       base split seed for Federated
#
# Notes:
# - Baseline is deterministic, so its 10 runs are identical, but we cache them
#   anyway to keep the interface symmetric.
# - DP betas are computed on genotypes harmonized to the effect allele to match
#   Baseline and Federated orientation.

source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(data.table); library(matrixStats)
})

# ---------- dirs and helpers ----------
d5 <- stage_dirs("5_MethodComparison")
dir.create(d5$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d5$out_t, showWarnings = FALSE, recursive = TRUE)

RUN2RUN_DIR <- file.path(DIR_GLOBAL, "run2run")
dir.create(RUN2RUN_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RUN2RUN_DIR, "BASE"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RUN2RUN_DIR, "DP"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(RUN2RUN_DIR, "FED"),  showWarnings = FALSE, recursive = TRUE)

MC_REPS            <- as.integer(Sys.getenv("MC_REPS", "10"))
PGS_BETA_SUBSAMPLE <- as.integer(Sys.getenv("PGS_BETA_SUBSAMPLE", "200"))
RUN2RUN_SEED       <- as.integer(Sys.getenv("RUN2RUN_SEED", "20251030"))
FED_SITES          <- as.integer(Sys.getenv("FED_SITES", "5"))
FED_SPLIT_SEED     <- as.integer(Sys.getenv("FED_SPLIT_SEED", "1"))

safe_png <- function(path, w=1100, h=800) png_open(path, w, h)

# ---------- load shared data ----------
ann      <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
gwas_dat <- readRDS(file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
base_gwas <- readRDS(file.path(DIR_GLOBAL, "baseline_gwas_chr22_full.rds"))

# Budgets for DP
bud <- if (file.exists(file.path(DIR_GLOBAL, "dp_budgets.rds"))) {
  readRDS(file.path(DIR_GLOBAL, "dp_budgets.rds"))
} else list(eps_gwas_bmi=2.5, eps_gwas_t2d=2.5, delta_total=1e-6)

# ---------- genotype harmonization to effect allele ----------
# Aligns G to "effect allele" preference used by Baseline and Federated.
harmonize_and_get_env <- function() {
  pgs_any  <- readRDS(file.path(DIR_GLOBAL, "pgs_any_chr22.rds"))
  cache_file <- file.path(DIR_GLOBAL, "alleles_cache.rds")
  if (!file.exists(cache_file)) stop("alleles_cache.rds not found. Run Stage 2 first.")

  G  <- gwas_dat$G
  rs <- ann$column_rsids

  # keep polymorphic columns
  keep <- which(matrixStats::colVars(G, na.rm = TRUE, useNames = FALSE) > 0)
  G <- G[, keep, drop = FALSE]
  rs <- rs[keep]

  alleles_vec <- readRDS(cache_file)[keep]
  alts <- vapply(strsplit(alleles_vec, "/", fixed=TRUE), `[`, "", 2)
  eff  <- alts
  m <- match(rs, pgs_any$rsID); has <- !is.na(m)
  eff[has] <- pgs_any$effect_allele[m[has]]
  flip <- eff != alts
  if (any(flip, na.rm = TRUE)) G[, flip] <- 2 - G[, flip]

  list(
    G = G,
    rs = rs,
    keep = keep,
    sd_g = matrixStats::colSds(G, na.rm = TRUE),
    age_z = as.numeric(scale(gwas_dat$cov$age)),
    sex01 = as.integer(gwas_dat$cov$sex),
    y_bmi = as.numeric(gwas_dat$yBMI),
    y_t2d = as.integer(gwas_dat$yT2D)
  )
}
env <- harmonize_and_get_env()

# ---------- choose a stable SNP subset per trait and cache it ----------
subset_file <- file.path(RUN2RUN_DIR, "subset_indices.rds")
if (file.exists(subset_file)) {
  subs <- readRDS(subset_file)
  sub_bmi <- subs$sub_bmi
  sub_t2d <- subs$sub_t2d
} else {
  effects <- readRDS(file.path(DIR_GLOBAL, "pgs_effects_chr22.rds"))
  map_pos <- function(causal_idx, keep) {
    pos <- match(causal_idx, keep)
    pos <- pos[!is.na(pos)]
    unique(pos)
  }
  pool_bmi <- map_pos(effects$bmi$causal_idx, env$keep)
  pool_t2d <- map_pos(effects$t2d$causal_idx, env$keep)
  # fallbacks if pools are tiny
  if (!length(pool_bmi)) pool_bmi <- seq_len(ncol(env$G))
  if (!length(pool_t2d)) pool_t2d <- seq_len(ncol(env$G))

  set.seed(RUN2RUN_SEED + 101)
  pick_subset <- function(pool, K) {
    pool <- pool[is.finite(pool)]
    if (length(pool) <= K) sort(unique(pool)) else sort(sample(pool, K))
  }
  sub_bmi <- pick_subset(pool_bmi, PGS_BETA_SUBSAMPLE)
  sub_t2d <- pick_subset(pool_t2d, PGS_BETA_SUBSAMPLE)

  saveRDS(list(sub_bmi = sub_bmi, sub_t2d = sub_t2d), subset_file)
}

rs_bmi <- env$rs[sub_bmi]
rs_t2d <- env$rs[sub_t2d]

# ---------- reference baseline betas on the same subset ----------
m_bmi <- match(rs_bmi, base_gwas$SNP)
m_t2d <- match(rs_t2d, base_gwas$SNP)
beta_base_bmi <- as.numeric(base_gwas$BETA_BMI[m_bmi])
beta_base_t2d <- as.numeric(base_gwas$BETA_T2D[m_t2d])

# ---------- method runners ----------
cache_path <- function(method, run_id) {
  file.path(RUN2RUN_DIR, method, sprintf("run_%02d.rds", run_id))
}

# Baseline: deterministic. Still cache for symmetry.
run_baseline_once <- function(run_id) {
  list(
    rs_bmi = rs_bmi,
    rs_t2d = rs_t2d,
    idx_bmi = sub_bmi,
    idx_t2d = sub_t2d,
    beta_bmi = beta_base_bmi,
    beta_t2d = beta_base_t2d,
    seed = NA_integer_
  )
}

# DP: DP-OLS for BMI, DP aggregated logistic for T2D, on harmonized G.
run_dp_once <- function(run_id) {
  source(file.path("3_DifferentialPrivacy", "R Code", "3a_dp_helpers.R"))
  set.seed(RUN2RUN_SEED + 1000 + run_id)

  p_keep <- ncol(env$G)
  eps_bmi   <- bud$eps_gwas_bmi / p_keep
  eps_t2d   <- bud$eps_gwas_t2d / p_keep
  delta_ps  <- if (is.null(bud$delta_total)) 1e-6 else bud$delta_total
  delta_bmi <- delta_ps / p_keep

  # BMI DP betas
  beta_bmi <- numeric(length(sub_bmi))
  for (j in seq_along(sub_bmi)) {
    i <- sub_bmi[j]
    g_z <- as.numeric(scale(env$G[, i]))
    X <- cbind(1, g_z, env$age_z, env$sex01)
    fit <- dp_ols_one(X, env$y_bmi, eps=eps_bmi, delta=delta_bmi, lambda=1.0, L=3.0, YB=80)
    s <- env$sd_g[i]; if (!is.finite(s) || s <= 0) s <- 1
    beta_bmi[j] <- fit$beta[2] / s
  }

  # T2D DP betas via 3x2 noisy counts
  beta_t2d <- numeric(length(sub_t2d))
  for (j in seq_along(sub_t2d)) {
    i <- sub_t2d[j]
    g_fac <- factor(env$G[, i], levels = c(0,1,2))
    tab <- table(g_fac, env$y_t2d)
    baseM <- matrix(as.numeric(tab), 3, 2)
    noisy <- baseM + matrix(rLaplace(6, 1/eps_t2d), 3, 2)
    noisy[!is.finite(noisy) | noisy < 0] <- 0
    noisy <- noisy + 1e-6
    df <- data.frame(g = 0:2, case = noisy[,2], ctrl = noisy[,1])
    df <- df[(df$case + df$ctrl) > 0, , drop=FALSE]
    if (nrow(df) < 2) { beta_t2d[j] <- NA_real_; next }
    fit <- try(suppressWarnings(glm(cbind(case, ctrl) ~ g, data=df, family=binomial())), silent=TRUE)
    beta_t2d[j] <- if (inherits(fit, "try-error") || !("g" %in% names(coef(fit)))) NA_real_ else unname(coef(fit)["g"])
  }

  list(
    rs_bmi = rs_bmi,
    rs_t2d = rs_t2d,
    idx_bmi = sub_bmi,
    idx_t2d = sub_t2d,
    beta_bmi = beta_bmi,
    beta_t2d = beta_t2d,
    seed = RUN2RUN_SEED + 1000 + run_id
  )
}

# Federated: meta analysis across K sites built from harmonized G.
run_fed_once <- function(run_id) {
  source(file.path("4_FederatedLearning", "R Code", "4a_fed_utils.R"))
  # Build the same structure as make_sites_from_gwas would, but with our harmonized G
  make_sites_local <- function(G, yBMI, yT2D, cov, K, seed) {
    n <- nrow(G)
    idxs <- split_indices(n, K, seed = seed)
    lapply(idxs, function(ii) list(
      idx = ii,
      G   = G[ii, , drop = FALSE],
      yT2D = yT2D[ii],
      yBMI = yBMI[ii],
      cov  = cov[ii, , drop = FALSE]
    ))
  }

  seed_here <- FED_SPLIT_SEED + run_id
  sites <- make_sites_local(env$G, env$y_bmi, env$y_t2d, gwas_dat$cov, K = FED_SITES, seed = seed_here)

  # BMI linear meta
  beta_bmi <- numeric(length(sub_bmi))
  for (j in seq_along(sub_bmi)) {
    i <- sub_bmi[j]
    b <- s <- numeric(0)
    for (k in seq_len(FED_SITES)) {
      res <- site_lm_linear(sites[[k]]$yBMI, sites[[k]]$G[, i], sites[[k]]$cov)
      if (res["ok"] == 1) { b <- c(b, res["beta"]); s <- c(s, res["se"]) }
    }
    cmb <- combine_meta(b, s)
    beta_bmi[j] <- cmb$beta
  }

  # T2D logistic meta
  beta_t2d <- numeric(length(sub_t2d))
  for (j in seq_along(sub_t2d)) {
    i <- sub_t2d[j]
    b <- s <- numeric(0)
    for (k in seq_len(FED_SITES)) {
      res <- site_glm_logistic(sites[[k]]$yT2D, sites[[k]]$G[, i], sites[[k]]$cov)
      if (res["ok"] == 1) { b <- c(b, res["beta"]); s <- c(s, res["se"]) }
    }
    cmb <- combine_meta(b, s)
    beta_t2d[j] <- cmb$beta
  }

  list(
    rs_bmi = rs_bmi,
    rs_t2d = rs_t2d,
    idx_bmi = sub_bmi,
    idx_t2d = sub_t2d,
    beta_bmi = beta_bmi,
    beta_t2d = beta_t2d,
    seed = seed_here
  )
}

# ---------- run or resume cached runs ----------
run_or_load <- function(method, run_id) {
  f <- cache_path(method, run_id)
  if (file.exists(f)) return(readRDS(f))
  obj <- switch(method,
    "BASE" = run_baseline_once(run_id),
    "DP"   = run_dp_once(run_id),
    "FED"  = run_fed_once(run_id),
    stop("Unknown method: ", method)
  )
  saveRDS(obj, f)
  obj
}

methods <- c("BASE","DP","FED")
for (m in methods) {
  for (r in seq_len(MC_REPS)) {
    cat(sprintf("[run2run] method=%s run=%d/%d\n", m, r, MC_REPS))
    invisible(run_or_load(m, r))
  }
}

# ---------- assemble Δβ and plot as pirate ----------
load_all <- function(method) {
  runs <- lapply(seq_len(MC_REPS), function(r) readRDS(cache_path(method, r)))
  list(
    BMI = do.call(cbind, lapply(runs, `[[`, "beta_bmi")),
    T2D = do.call(cbind, lapply(runs, `[[`, "beta_t2d"))
  )
}

all_BASE <- load_all("BASE")
all_DP   <- load_all("DP")
all_FED  <- load_all("FED")

# Δβ = method beta minus baseline beta (per run, per SNP)
# Flatten across runs and SNPs for pirate distributions
as_deltas <- function(mat_method, beta_base_vec, label) {
  # mat_method: rows=SNPs, cols=runs
  if (is.null(dim(mat_method))) mat_method <- matrix(mat_method, ncol = MC_REPS)
  del <- sweep(mat_method, 1, beta_base_vec, FUN = "-")
  data.frame(method = label, delta = as.numeric(del))
}

df_BMI <- rbind(
  as_deltas(all_BASE$BMI, beta_base_bmi,  "Baseline (BMI)"),
  as_deltas(all_DP$BMI,   beta_base_bmi,  "DP (BMI)"),
  as_deltas(all_FED$BMI,  beta_base_bmi,  "Federated (BMI)")
)
df_T2D <- rbind(
  as_deltas(all_BASE$T2D, beta_base_t2d,  "Baseline (T2D)"),
  as_deltas(all_DP$T2D,   beta_base_t2d,  "DP (T2D)"),
  as_deltas(all_FED$T2D,  beta_base_t2d,  "Federated (T2D)")
)

plot_pirate <- function(dat, title, outfile) {
  methods <- unique(dat$method)
  safe_png(outfile, 1100, 800)
  par(mar=c(8,6,4,2)+0.1)
  y <- dat$delta[is.finite(dat$delta)]
  yr <- if (length(y)) range(y, na.rm = TRUE) else c(-1, 1)
  plot(0,0,type="n", xlim=c(0.5, length(methods)+0.5), ylim=yr,
       xaxt="n", xlab="", ylab=expression(Delta*beta), main=title)
  for (i in seq_along(methods)) {
    v <- dat$delta[dat$method == methods[i]]
    v <- v[is.finite(v)]
    if (!length(v)) next
    xj <- jitter(rep(i, length(v)), amount=0.25)
    points(xj, v, pch=16, cex=0.25, col=adjustcolor("black", 0.2))
    segments(i-0.35, mean(v, na.rm=TRUE), i+0.35, mean(v, na.rm=TRUE), lwd=3)
  }
  axis(1, at=seq_along(methods), labels=methods, las=2)
  dev.off()
}

out_png <- file.path(d5$out_g, "run2run_pirate.png")
plot_pirate(rbind(df_BMI, df_T2D), "Run-to-run variability: Δβ across methods and traits", out_png)

# ---------- write a small text summary ----------
writeLines(c(
  sprintf("Runs per method: %d", MC_REPS),
  sprintf("SNPs per trait:  %d", PGS_BETA_SUBSAMPLE),
  sprintf("Sites (FED):     %d", FED_SITES),
  sprintf("Seed base:       %d", RUN2RUN_SEED),
  "Cache dir: Data/global/run2run",
  paste("Output plot:", basename(out_png))
), file.path(d5$out_t, "5e_outputs.txt"))

cat("✓ Run-to-run variability complete. Plot and cache saved.\n")
