# 0c_filter.R (patched)
source(file.path("R", "config.R"))
suppressPackageStartupMessages({
  library(SNPRelate); library(gdsfmt); library(matrixStats); library(data.table)
})

p <- readRDS(file.path(DIR_GLOBAL, "paths_gds.rds"))
stopifnot(file.exists(p$gds_file))

local({
  g <- SNPRelate::snpgdsOpen(p$gds_file, readonly = TRUE)
  on.exit(SNPRelate::snpgdsClose(g), add = TRUE)

  cat("[0c] Opened:", p$gds_file, "\n")
  # Verify required nodes (no root fiddling)
  need <- c("sample.id","snp.id","snp.chromosome","snp.position","snp.rs.id")
  missing <- need[vapply(need, function(nm) is.null(gdsfmt::index.gdsn(g, nm, silent = TRUE)), logical(1))]
  if (length(missing)) stop("Missing nodes in GDS: ", paste(missing, collapse = ", "))

  # Read IDs & counts
  snp_ids   <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.id"))
  n_gds_raw <- length(snp_ids)

  # Filters: SelectSNP returns IDs; map to indices
  keep_ids <- SNPRelate::snpgdsSelectSNP(g, maf = 0.01, missing.rate = 0.05)
  keep_idx <- match(keep_ids, snp_ids); keep_idx <- keep_idx[!is.na(keep_idx)]
  n_maf_ok <- length(keep_idx); n_ld_keep <- n_maf_ok  # (no LD pruning here)

  # Samples & covariates (toy)
  samples <- as.character(gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "sample.id")))
  set.seed(42)
  covariates <- data.frame(
    age = round(pmin(pmax(rnorm(length(samples), 40, 15), 18), 90)),
    sex = rbinom(length(samples), size = 1L, prob = 0.5),
    row.names = samples
  )
  saveRDS(list(cov = covariates, eur_samples = samples),
          file.path(DIR_GLOBAL, "covariates_chr22.rds"))
  writeLines(samples, file.path(DIR_GLOBAL, "eur_samples_ids.txt"))

  # SNP annotations (subset by indices)
  CHR <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.chromosome"))[keep_idx]
  BP  <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.position"))  [keep_idx]
  all_rs <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(g, "snp.rs.id"))
  col_rs <- all_rs[keep_idx]

  saveRDS(
    list(
      column_rsids = col_rs,
      CHR = CHR, BP = BP,
      snp_keep_ids = keep_ids,
      snp_keep_idx = keep_idx,
      # back-compat for downstream code that previously used 'snp_keep'
      snp_keep     = keep_idx,
      n_gds_raw = n_gds_raw,
      n_maf_ok = n_maf_ok,
      n_ld_keep = n_ld_keep
    ),
    file.path(DIR_GLOBAL, "snp_annotations_chr22.rds")
  )

  cat("✓ Stage 0 artifacts saved in ", DIR_GLOBAL, "\n", sep = "")
})
