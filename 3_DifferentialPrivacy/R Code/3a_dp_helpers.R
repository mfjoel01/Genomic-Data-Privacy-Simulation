## 3a_dp_helpers.R
## DP helpers for Stage 3
## - Analyze Gauss style DP-OLS with noise-aware SEs and DP-safe SSE
## - Composition helpers (strong and zCDP)
## - One-shot DP top-K selection (Gumbel) on bounded marginal scores
## - DP ROC for PGS figures (unchanged interface)

# ----------------- Basic RNG and clipping -----------------
rLaplace <- function(n, scale) {
  u <- runif(n, -0.5, 0.5)
  -scale * sign(u) * log(1 - 2 * abs(u))
}

clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

clip_rows_l2 <- function(X, L) {
  for (i in seq_len(nrow(X))) {
    nrm <- sqrt(sum(X[i, ]^2))
    s <- min(1, L / (nrm + 1e-12))
    X[i, ] <- X[i, ] * s
  }
  X
}

sym_gauss <- function(p, sigma) {
  M  <- matrix(0, p, p)
  ut <- upper.tri(M, diag = TRUE)
  z  <- rnorm(sum(ut), 0, sigma)
  M[ut] <- z
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

ensure_spd <- function(M, jitter_start = 1e-6, max_iter = 8L, floor = 1e-8) {
  M <- 0.5 * (M + t(M))
  p <- nrow(M)
  jit <- jitter_start
  for (i in seq_len(max_iter)) {
    ok <- try(chol(M + diag(jit, p)), silent = TRUE)
    if (!inherits(ok, "try-error")) return(M + diag(jit, p))
    jit <- jit * 10
  }
  ev <- eigen(M, symmetric = TRUE)
  ev$values[ev$values < floor] <- floor
  ev$vectors %*% (ev$values * t(ev$vectors))
}

# ----------------- Composition helpers -----------------
# Strong composition (advanced composition) per-query allocation
ac_per_query <- function(eps_total, m, delta_total, delta_prime_frac = 0.5) {
  stopifnot(eps_total > 0, m >= 1, delta_total > 0)
  delta_prime <- max(1e-12, delta_total * delta_prime_frac)
  B <- sqrt(2 * m * log(1 / delta_prime))
  eps_per <- ( -B + sqrt(B * B + 4 * m * eps_total) ) / (2 * m)
  delta_per <- (delta_total * (1 - delta_prime_frac)) / m
  list(eps = max(eps_per, 1e-6), delta = max(delta_per, 1e-12))
}

# zCDP allocator: eps_total, delta_total -> rho_total, then per-query rho
# Uses eps = rho + 2*sqrt(rho*log(1/delta)) to back out rho_total
ac_per_query_zcdp <- function(eps_total, delta_total, m) {
  stopifnot(eps_total > 0, delta_total > 0, m >= 1)
  L <- log(1 / delta_total)
  t <- -sqrt(L) + sqrt(L + eps_total)         # t = sqrt(rho_total)
  rho_total <- max(t^2, 1e-16)
  rho_per   <- rho_total / m
  list(rho_per = rho_per, rho_total = rho_total)
}

# ----------------- DP OLS (Analyze Gauss style) -----------------
# Returns beta, se, t, p, sse, sigma2
dp_ols_one <- function(X, y,
                       eps = NULL, delta = NULL,   # classic Gaussian mechanism
                       rho = NULL,                 # zCDP parameter per query
                       lambda = 1.0, L = 3.0, YB = 80,
                       wS = 0.2, wt = 0.7, wy = 0.1,
                       debug = NULL, dbg_id = NA_character_) {
  n <- nrow(X); p <- ncol(X)
  stopifnot(n > p, L > 0, YB > 0)

  # Preprocess with row L2 clipping and y clipping
  Xc <- clip_rows_l2(X, L)
  yc <- clip(y, -YB, YB)

  # Sufficient statistics
  S  <- crossprod(Xc)
  t_ <- crossprod(Xc, yc)
  yy <- sum(yc * yc)

  # Sensitivities for Gaussian mechanism
  sens_S  <- L^2
  sens_t  <- L * YB
  sens_yy <- YB^2

  if (!is.null(rho)) {
    sig_S  <- sens_S  / sqrt(2 * rho * wS)
    sig_t  <- sens_t  / sqrt(2 * rho * wt)
    sig_yy <- sens_yy / sqrt(2 * rho * wy)
  } else {
    stopifnot(!is.null(eps), !is.null(delta), eps > 0, delta > 0)
    w <- c(wS, wt, wy); w <- w / sum(w); wS <- w[1]; wt <- w[2]; wy <- w[3]
    k <- sqrt(2 * log(1.25 / delta))
    sig_S  <- sens_S  * k / (eps * wS)
    sig_t  <- sens_t  * k / (eps * wt)
    sig_yy <- sens_yy * k / (eps * wy)
  }

  # Privatize stats
  S_t  <- S  + sym_gauss(p, sig_S)
  t_t  <- t_ + matrix(rnorm(p, 0, sig_t), p, 1)
  yy_t <- yy + rnorm(1, 0, sig_yy)

  # Ridge solve with SPD guard
  K    <- ensure_spd(S_t + diag(lambda, p))
  beta <- tryCatch(solve(K, t_t), error = function(e) rep(NA_real_, p))

  # DP-safe SSE: SSE = y'y - 2 beta' X'y + beta' (X'X) beta, all DP quantities
  SSE <- as.numeric(yy_t - 2 * drop(crossprod(beta, t_t)) + drop(t(beta) %*% S_t %*% beta))
  dof <- max(1, n - p)
  sse_fallback <- FALSE
  if (!is.finite(SSE) || SSE <= 0) {
    # No non-DP fallback. Use a conservative DP fallback from yy_t only.
    sse_fallback <- TRUE
    SSE <- max(1e-6, abs(yy_t))
  }
  sigma2 <- SSE / dof

  # Approximate covariance and Wald tests
  Vb <- tryCatch(solve(K), error = function(e) diag(1e6, p))
  se <- sqrt(pmax(diag(Vb), 1e-12) * sigma2)
  se[!is.finite(se) | se <= 0] <- 1e6
  tvals <- as.numeric(beta / se)
  pvals <- 2 * pt(-abs(tvals), df = dof)

  if (!is.null(debug)) {
    # Append one record per call
    try({
      eigS <- try(eigen(S_t, only.values = TRUE)$values, silent = TRUE)
      dbg <- list(
        id = dbg_id, n = n, p = p,
        eps = eps, delta = delta, rho = rho,
        sig_S = sig_S, sig_t = sig_t, sig_yy = sig_yy,
        eigSmin = if (inherits(eigS, "try-error")) NA_real_ else min(eigS),
        condK = tryCatch(kappa(K), error = function(e) NA_real_),
        sse = SSE, sse_fallback = sse_fallback,
        beta = as.numeric(beta), se = se, t = tvals, p = pvals[2]
      )
      debug[[length(debug) + 1L]] <<- dbg
    }, silent = TRUE)
  }

  list(beta = as.numeric(beta), se = se, t = tvals, p = pvals, sse = SSE, sigma2 = sigma2)
}

# ----------------- One-shot DP top-K selection (Gumbel) -----------------
# Utility: bounded marginal score u_i = |sum_j ( (g_ij - 1) * clip(y_j, -YB, YB) )|
# Sensitivity per record is YB, so Gumbel scale b = YB / eps_sel
# Returns integer indices into columns of G (1-based)
.sample_gumbel <- function(n, scale) { -scale * log(-log(runif(n))) }

dp_topk_by_product <- function(G, y, K, eps_sel, YB_sel) {
  K <- as.integer(K); if (!is.finite(K) || K < 1) return(integer(0))
  eps_sel <- max(1e-8, as.numeric(eps_sel))
  y_c <- clip(as.numeric(y), -YB_sel, YB_sel)
  # sum((g-1)*y_c) = crossprod(y_c, G) - sum(y_c)
  s <- as.numeric(crossprod(y_c, G)) - sum(y_c)
  u <- abs(s)
  b <- YB_sel / eps_sel
  noisy <- u + .sample_gumbel(length(u), b)
  ord <- order(noisy, decreasing = TRUE)
  ord[seq_len(min(K, length(ord)))]
}

# ----------------- DP ROC via one-shot histograms (for PGS figs) -----------------
make_dp_roc <- function(scores, y, eps_total, n_pts = 60, eps_den_frac = 0.1) {
  scores <- as.numeric(scores); y <- as.integer(y)
  ok <- is.finite(scores) & is.finite(y)
  scores <- scores[ok]; y <- y[ok]
  if (!length(scores)) return(list(FPR = c(0,1), TPR = c(0,1)))

  edges <- quantile(scores, probs = seq(0, 1, length.out = n_pts + 1), na.rm = TRUE)
  edges <- unique(as.numeric(edges))
  if (length(edges) < 2L) edges <- unique(range(scores, na.rm = TRUE))
  bins  <- length(edges) - 1L

  pos <- scores[y == 1]; neg <- scores[y == 0]
  bp  <- cut(pos, edges, include.lowest = TRUE, labels = FALSE, right = TRUE)
  bn  <- cut(neg, edges, include.lowest = TRUE, labels = FALSE, right = TRUE)

  cnt_pos <- tabulate(bp, nbins = bins)
  cnt_neg <- tabulate(bn, nbins = bins)

  eps_den <- max(1e-12, eps_total * eps_den_frac)
  eps_vec <- max(1e-12, eps_total - eps_den)

  n_pos_dp <- pmax(1, round(sum(cnt_pos) + rLaplace(1, 1/eps_den)))
  n_neg_dp <- pmax(1, round(sum(cnt_neg) + rLaplace(1, 1/eps_den)))

  cnt_pos_dp <- pmax(0, cnt_pos + rLaplace(bins, 1/eps_vec))
  cnt_neg_dp <- pmax(0, cnt_neg + rLaplace(bins, 1/eps_vec))

  TPR <- cumsum(rev(cnt_pos_dp)) / as.numeric(n_pos_dp)
  FPR <- cumsum(rev(cnt_neg_dp)) / as.numeric(n_neg_dp)
  list(FPR = rev(FPR), TPR = rev(TPR))
}

# ----------------- Optional DP 2x2 chi-square (kept for reference) -----------------
dp_chisq_p <- function(g, y01, eps_cell) {
  g_fac <- factor(g, levels = c(0,1))
  y_fac <- factor(y01, levels = c(0,1))
  tab   <- table(g_fac, y_fac)
  noisy <- matrix(as.numeric(tab), 2, 2) + matrix(rLaplace(4, 1/eps_cell), 2, 2)
  noisy[!is.finite(noisy) | noisy < 0] <- 0
  noisy <- noisy + 1e-6
  rs <- rowSums(noisy); cs <- colSums(noisy); tot <- sum(noisy)
  if (!is.finite(tot) || tot <= 0) return(1.0)
  E <- outer(rs, cs, function(r, c) r * c / tot); E[E <= 0] <- 1e-12
  stat <- sum((noisy - E)^2 / E)
  pchisq(stat, df = 1, lower.tail = FALSE)
}
