## 4a_fed_utils.R — helpers for federated GWAS/PGS (no raw data leaves sites)
suppressPackageStartupMessages({
  library(data.table)
})

# Deterministic split of n rows into K sites (roughly equal sizes)
split_indices <- function(n, K, seed=1) {
  stopifnot(K >= 2)
  set.seed(seed)
  idx <- sample.int(n, n)  # shuffle
  sizes <- rep(floor(n/K), K); sizes[seq_len(n %% K)] <- sizes[seq_len(n %% K)] + 1L
  starts <- c(0L, cumsum(sizes))
  lapply(seq_len(K), function(k) idx[(starts[k]+1):starts[k+1]])
}

# Turn full gwas_dat into a list of site-local objects
make_sites_from_gwas <- function(gwas_dat, K, seed=1) {
  G <- gwas_dat$G
  yT2D <- gwas_dat$yT2D
  yBMI <- gwas_dat$yBMI
  cov  <- gwas_dat$cov
  n <- nrow(G)
  idxs <- split_indices(n, K, seed)
  sites <- vector("list", K)
  for (k in seq_len(K)) {
    ii <- idxs[[k]]
    sites[[k]] <- list(
      idx = ii,
      G   = G[ii, , drop = FALSE],
      yT2D = yT2D[ii],
      yBMI = yBMI[ii],
      cov  = cov[ii, , drop = FALSE]
    )
  }
  sites
}

# ----- Site-level regressions (local, hardened) -----

site_glm_logistic <- function(y, g, cov) {
  # y: binary (0/1); g: numeric vector; cov: data.frame with columns age, sex
  y <- as.integer(y); g <- as.numeric(g)
  # quick exits for degenerate sites
  if (!any(is.finite(g)) || var(g, na.rm=TRUE) == 0 ||
      sum(y == 1, na.rm=TRUE) == 0 || sum(y == 0, na.rm=TRUE) == 0) {
    return(c(beta=NA_real_, se=NA_real_, n=length(y), ok=0))
  }
  df <- data.frame(y=y, g=g, age=cov$age, sex=cov$sex)

  # primary: standard GLM (suppress separation/rank warnings)
  fit <- try(suppressWarnings(glm(y ~ g + age + sex, data=df, family=binomial())), silent=TRUE)
  if (!inherits(fit, "try-error")) {
    co <- try(summary(fit)$coefficients, silent=TRUE)
    if (!inherits(co, "try-error") && is.finite(unname(co["g","Std. Error"]))) {
      return(c(beta = unname(co["g","Estimate"]),
               se   = unname(co["g","Std. Error"]),
               n    = nrow(df), ok=1))
    }
  }

  # fallback: penalized logistic regression (handles separation)
  if (requireNamespace("logistf", quietly = TRUE)) {
    fit2 <- try(logistf::logistf(y ~ g + age + sex, data=df), silent=TRUE)
    if (!inherits(fit2, "try-error")) {
      betag <- suppressWarnings(as.numeric(coef(fit2)["g"]))
      seg   <- suppressWarnings(as.numeric(fit2$se["g"]))
      return(c(beta = betag, se = seg, n = nrow(df), ok=1))
    }
  }

  c(beta=NA_real_, se=NA_real_, n=nrow(df), ok=0)
}

site_lm_linear <- function(y, g, cov) {
  g <- as.numeric(g)
  if (!any(is.finite(g)) || var(g, na.rm=TRUE) == 0) {
    return(c(beta=NA_real_, se=NA_real_, n=length(y), ok=0))
  }
  df <- data.frame(y=y, g=g, age=cov$age, sex=cov$sex)
  out <- tryCatch({
    m <- lm(y ~ g + age + sex, data=df)
    co <- summary(m)$coefficients
    c(beta = unname(co["g","Estimate"]), se = unname(co["g","Std. Error"]), n = nrow(df), ok=1)
  }, error=function(e) c(beta=NA_real_, se=NA_real_, n=length(y), ok=0))
  out
}

# ----- Meta-analysis combiner -----
combine_meta <- function(beta, se) {
  ok <- which(is.finite(beta) & is.finite(se) & se > 0)
  k <- length(ok)
  if (k == 0) return(list(beta=NA_real_, se=NA_real_, z=NA_real_, p=NA_real_, Q=NA_real_, I2=NA_real_, p_het=NA_real_, k=0))
  b <- beta[ok]; s <- se[ok]; w <- 1/(s^2)
  bhat <- sum(w*b) / sum(w)
  sehat <- sqrt(1/sum(w))
  z <- bhat / sehat
  p <- 2*pnorm(-abs(z))
  Q <- sum(w * (b - bhat)^2)
  df <- max(1, k-1)
  I2 <- max(0, min(1, (Q - df)/Q))
  p_het <- pchisq(Q, df=df, lower.tail = FALSE)
  list(beta=bhat, se=sehat, z=z, p=p, Q=Q, I2=I2, p_het=p_het, k=k)
}

lambda_gc <- function(pvals) {
  pv <- pvals[is.finite(pvals)]
  if (!length(pv)) return(NA_real_)
  chisq <- qchisq(pmin(pmax(1 - pv, .Machine$double.eps), 1 - 1e-15), df=1)
  median(chisq, na.rm=TRUE) / qchisq(0.5, df=1)
}

# ----- Aggregated OLS for PGS performance -----

ols_cross_stats <- function(X, y) {
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  list(XtX = XtX, Xty = Xty, yty = sum(y*y), sumy = sum(y), n = length(y))
}

sum_ols_stats <- function(stats_list) {
  p <- nrow(stats_list[[1]]$XtX)
  XtX <- matrix(0, p, p)
  Xty <- matrix(0, p, 1)
  yty <- 0; sumy <- 0; n <- 0
  for (s in stats_list) {
    XtX <- XtX + s$XtX
    Xty <- Xty + s$Xty
    yty <- yty + s$yty
    sumy <- sumy + s$sumy
    n <- n + s$n
  }
  list(XtX=XtX, Xty=Xty, yty=yty, sumy=sumy, n=n)
}

r2_from_ols_sums <- function(sumstats) {
  # Fit beta = solve(X'X, X'y); compute SSE and SST
  beta <- tryCatch(solve(sumstats$XtX, sumstats$Xty), error=function(e) rep(NA_real_, nrow(sumstats$XtX)))
  beta <- as.numeric(beta)
  ybar <- sumstats$sumy / sumstats$n
  # SSE = y'y - 2 beta^T X'y + beta^T X'X beta
  sse <- sumstats$yty - 2*sum(beta * sumstats$Xty) + as.numeric(t(beta) %*% sumstats$XtX %*% beta)
  sst <- sumstats$yty - sumstats$n * ybar^2
  r2  <- if (sst > 0) 1 - sse/sst else NA_real_
  p   <- nrow(sumstats$XtX)
  r2_adj <- if (sst > 0) 1 - (sse/(sumstats$n - p)) / (sst/(sumstats$n - 1)) else NA_real_
  list(beta = beta, SSE = sse, SST = sst, R2 = r2, R2_adj = r2_adj)
}

# ----- Federated ROC from aggregated counts -----

global_thresholds_from_sites <- function(scores_list, probs = seq(0,1,length.out=61), weights = NULL) {
  # Compute site-local quantiles, then weighted average across sites
  K <- length(scores_list)
  Q <- sapply(scores_list, function(v) as.numeric(quantile(v, probs=probs, na.rm=TRUE)))
  if (is.null(weights)) weights <- rep(1/K, K) else weights <- weights/sum(weights)
  # Weighted average across columns (sites)
  drop(Q %*% matrix(weights, nrow=K))
}

aggregate_roc_counts <- function(scores_list, y_list, thresholds) {
  K <- length(scores_list)
  n_thr <- length(thresholds)
  tp <- rep(0, n_thr); fp <- rep(0, n_thr); pos <- 0; neg <- 0
  for (k in seq_len(K)) {
    s <- scores_list[[k]]; y <- y_list[[k]]
    pos <- pos + sum(y == 1); neg <- neg + sum(y == 0)
    for (i in seq_len(n_thr)) {
      t <- thresholds[i]
      tp[i] <- tp[i] + sum(y == 1 & s >= t)
      fp[i] <- fp[i] + sum(y == 0 & s >= t)
    }
  }
  TPR <- if (pos > 0) tp/pos else rep(0, n_thr)
  FPR <- if (neg > 0) fp/neg else rep(0, n_thr)
  ord <- order(FPR, TPR) # ensure sorted for trapezoid
  list(FPR = pmin(pmax(cummax(FPR[ord]), 0), 1),
       TPR = pmin(pmax(cummax(TPR[ord]), 0), 1))
}

## Robust trapezoid AUC
auc_trap <- function(FPR, TPR) {
  # Sort by FPR (then TPR), clamp to [0,1], and make sure we have ≥2 points.
  o <- order(FPR, TPR, na.last = NA)
  F <- pmin(pmax(FPR[o], 0), 1)
  T <- pmin(pmax(TPR[o], 0), 1)
  n <- length(F)
  if (n < 2L) return(NA_real_)
  # Standard trapezoid rule: midpoints length is (n-1), same as diff(F).
  dF   <- diff(F)
  Tmid <- 0.5 * (T[-1] + T[-n])
  auc  <- sum(dF * Tmid)
  # Be conservative with tiny numeric drift
  pmin(pmax(auc, 0), 1)
}
