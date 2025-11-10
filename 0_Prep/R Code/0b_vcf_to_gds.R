## ===== [0b] VCF → GDS (idempotent) =====
source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(SNPRelate))

vcf_file <- file.path(DIR_VCF, "chr22_EUR_sim_10k.controls.vcf.gz")
idx_file <- paste0(vcf_file, ".tbi")

dir.create(DIR_GDS, showWarnings = FALSE, recursive = TRUE)
gds_file <- file.path(DIR_GDS, "chr22_EUR_sim_10k.controls.gds")

cat("[0b] VCF → GDS (skip if exists)\n")

need_rebuild <- TRUE
if (file.exists(gds_file) && file.size(gds_file) > 0L) {
  cat("• Found existing GDS:", gds_file, "— validating...\n")
  need_rebuild <- FALSE
  tryCatch({
    SNPRelate::snpgdsSummary(gds_file)  # quick sanity check
  }, error = function(e) {
    cat("! Existing GDS failed validation — will rebuild. Reason:", conditionMessage(e), "\n")
    need_rebuild <<- TRUE
  })
  
  ## Optional: rebuild if source VCF is newer than the GDS
  # if (!need_rebuild && file.exists(vcf_file) && file.mtime(vcf_file) > file.mtime(gds_file)) {
  #   cat("↻ VCF is newer than GDS — rebuilding…\n")
  #   need_rebuild <- TRUE
  # }
}

if (need_rebuild) {
  ## Only require VCF/TBI if we truly need to rebuild
  stopifnot(file.exists(vcf_file), file.exists(idx_file))
  if (file.exists(gds_file)) suppressWarnings(try(unlink(gds_file), silent = TRUE))
  cat("→ Converting VCF to GDS…\n")
  SNPRelate::snpgdsVCF2GDS(
    vcf.fn = vcf_file,
    out.fn = gds_file,
    method = "biallelic.only",
    ignore.chr.prefix = "chr"  # set "" if your CHROM is plain "22"
  )
  cat("✓ Rebuilt:", gds_file, "\n")
}

## Show summary (existing or rebuilt) and persist path for downstream steps
SNPRelate::snpgdsSummary(gds_file)
saveRDS(list(gds_file = gds_file), file.path(DIR_GLOBAL, "paths_gds.rds"))
cat("✓ Using GDS:", gds_file, "\n")
## ===== end [0b] =====
