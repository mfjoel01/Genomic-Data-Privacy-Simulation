source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(SNPRelate); library(gdsfmt); library(data.table); library(matrixStats); library(pROC)
})

paths <- readRDS(file.path(DIR_GLOBAL, "paths_gds.rds"))
ann   <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
covs  <- readRDS(file.path(DIR_GLOBAL, "covariates_chr22.rds"))

pgs_t2d_chr22 <- readRDS(file.path(DIR_GLOBAL, "pgs_t2d_chr22.rds"))
pgs_bmi_chr22 <- readRDS(file.path(DIR_GLOBAL, "pgs_bmi_chr22.rds"))

# ---------- Open GDS and pull genotype matrix for retained SNPs ----------
g_res <- local({
  genofile <- snpgdsOpen(paths$gds_file, readonly = TRUE)
  on.exit(snpgdsClose(genofile), add = TRUE)

  req <- c("sample.id","snp.id","snp.allele","snp.rs.id","snp.chromosome","snp.position")
  miss <- req[vapply(req, function(nm) is.null(index.gdsn(genofile, nm, silent=TRUE)), logical(1))]
  if (length(miss)) stop("GDS missing nodes: ", paste(miss, collapse=", "))

  snp_keep_ids <- ann$snp_keep_ids    # IDs to read
  snp_keep_idx <- ann$snp_keep_idx    # indices for annotation

  eur_samples <- covs$eur_samples

  geno_mat <- snpgdsGetGeno(
    genofile,
    sample.id   = eur_samples,
    snp.id      = snp_keep_ids,
    snpfirstdim = FALSE,
    with.id     = FALSE
  )

  all_rsids    <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
  column_rsids <- all_rsids[snp_keep_idx]

  list(geno_mat = geno_mat, column_rsids = column_rsids)
})

G            <- g_res$geno_mat
column_rsids <- g_res$column_rsids
n            <- nrow(G)
p            <- ncol(G)

#cat(sprintf("Genotype matrix: %d samples X %d SNPs\n", n, p))

# ---------- Collapse effect-size tables by rsID (used by PGS mode) ----------
w_t2d_dt <- data.table(pgs_t2d_chr22[, .(rsID, effect_weight)])[ , .(weight = sum(effect_weight)), by = rsID]
w_bmi_dt <- data.table(pgs_bmi_chr22[, .(rsID, effect_weight)])[ , .(weight = sum(effect_weight)), by = rsID]

# ---------- T2D genetic component (always from PGS overlaps) ----------
causal_idx_t2d <- match(w_t2d_dt$rsID, column_rsids); valid_t2d <- which(!is.na(causal_idx_t2d))
causal_idx_t2d <- causal_idx_t2d[valid_t2d]
effect_sz_t2d  <- w_t2d_dt$weight[valid_t2d]
stopifnot(length(causal_idx_t2d) > 0)

g_raw_t2d <- as.vector(G[, causal_idx_t2d, drop = FALSE] %*% effect_sz_t2d)
g_std_t2d <- as.vector(scale(g_raw_t2d))

# =========================
# BMI genetic component mode
#   BMI_MODE = "PGS"    -> use PGS004994 overlaps (original behavior)
#   BMI_MODE = "SPARSE" -> use a small, random, disjoint causal set on chr22
#   BMI_MODE = "NULL"   -> no genetic signal for BMI (calibration)
# =========================
BMI_MODE <- toupper(Sys.getenv("BMI_MODE", "PGS"))

if (BMI_MODE == "PGS") {
  causal_idx_bmi <- match(w_bmi_dt$rsID, column_rsids); valid_bmi <- which(!is.na(causal_idx_bmi))
  causal_idx_bmi <- causal_idx_bmi[valid_bmi]
  effect_sz_bmi  <- w_bmi_dt$weight[valid_bmi]
  stopifnot(length(causal_idx_bmi) > 0)
  g_raw_bmi <- as.vector(G[, causal_idx_bmi, drop = FALSE] %*% effect_sz_bmi)
} else if (BMI_MODE == "SPARSE") {
  set.seed(77)
  pgs_rsids_bmi <- unique(pgs_bmi_chr22$rsID)
  not_pgs <- which(!(column_rsids %in% pgs_rsids_bmi))
  pool <- if (length(not_pgs) > 500) not_pgs else seq_len(p)

  K <- 300L          # number of causal SNPs
  eff_sd <- 0.1      # tune for desired h2
  causal_idx_bmi <- sort(sample(pool, K))
  effect_sz_bmi  <- rnorm(K, mean = 0, sd = eff_sd)
  g_raw_bmi <- as.vector(G[, causal_idx_bmi, drop = FALSE] %*% effect_sz_bmi)
} else if (BMI_MODE == "NULL") {
  causal_idx_bmi <- integer(0)
  effect_sz_bmi  <- numeric(0)
  g_raw_bmi <- rep(0, n)
} else {
  stop("Unknown BMI_MODE: ", BMI_MODE)
}

# Safe scaling (avoid NA when vector is constant)
g_std_bmi <- {
  v <- as.vector(scale(g_raw_bmi))
  if (anyNA(v)) rep(0, n) else v
}

# ---------- Phenotype parameters ----------
h2_bmi <- as.numeric(Sys.getenv("H2_BMI_OVERRIDE", NA))
if (is.na(h2_bmi)) h2_bmi <- if (BMI_MODE == "NULL") 0 else 0.002
h2_t2d <- 0.1
prev   <- 0.119

age_coef_t2d <- 0.05; sex_coef_t2d <- 0.35
age_coef_bmi <- 0.02; sex_coef_bmi <- 2.0

var_e_bmi  <- if (h2_bmi > 0) (1 - h2_bmi) / h2_bmi else 1e6
var_e_liab <- (1 - h2_t2d) / h2_t2d

covariates <- covs$cov

# ---------- Build GWAS phenotypes ----------
set.seed(42)
bmi_gwas  <- 25 + g_std_bmi + age_coef_bmi * covariates$age + sex_coef_bmi * covariates$sex + rnorm(n, 0, sqrt(var_e_bmi))
liab_gwas <- g_std_t2d + age_coef_t2d * covariates$age + sex_coef_t2d * covariates$sex + rnorm(n, 0, sqrt(var_e_liab))
t2d_gwas  <- as.integer(liab_gwas > quantile(liab_gwas, 1 - prev))

# ---------- Build PGS evaluation phenotypes (independent noise) ----------
set.seed(43)
bmi_pgs  <- 25 + g_std_bmi + age_coef_bmi * covariates$age + sex_coef_bmi * covariates$sex + rnorm(n, 0, sqrt(var_e_bmi))
liab_pgs <- g_std_t2d + age_coef_t2d * covariates$age + sex_coef_t2d * covariates$sex + rnorm(n, 0, sqrt(var_e_liab))
t2d_pgs  <- as.integer(liab_pgs > quantile(liab_pgs, 1 - prev))

# ---------- Persist datasets ----------
gwas_dat <- list(G = G, yT2D = t2d_gwas, yBMI = bmi_gwas, cov = covariates)
pgs_dat  <- list(G = G, yT2D = t2d_pgs,  yBMI = bmi_pgs,  cov = covariates)

saveRDS(gwas_dat, file.path(DIR_GLOBAL, "gwas_dat_chr22_full.rds"))
saveRDS(pgs_dat,  file.path(DIR_GLOBAL, "pgs_dat_chr22.rds"))

# Save effect indices (reflect actual BMI mode)
saveRDS(list(
  t2d = list(causal_idx = causal_idx_t2d, effect_sz = effect_sz_t2d),
  bmi = list(causal_idx = causal_idx_bmi, effect_sz = effect_sz_bmi)
), file.path(DIR_GLOBAL, "pgs_effects_chr22.rds"))

cat(sprintf("✓ Saved: gwas_dat_chr22_full.rds, pgs_dat_chr22.rds, pgs_effects_chr22.rds\n"))
cat(sprintf("✓ BMI_MODE=%s | h2_bmi=%.3f | causal_bmi=%d | causal_t2d=%d\n",
            BMI_MODE, h2_bmi, length(causal_idx_bmi), length(causal_idx_t2d)))
