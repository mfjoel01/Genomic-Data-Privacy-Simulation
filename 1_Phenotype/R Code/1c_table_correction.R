source(file.path("R", "config.R"))
suppressPackageStartupMessages(library(data.table))

ann <- readRDS(file.path(DIR_GLOBAL, "snp_annotations_chr22.rds"))
pgs_t2d_chr22 <- readRDS(file.path(DIR_GLOBAL, "pgs_t2d_chr22.rds"))
pgs_bmi_chr22 <- readRDS(file.path(DIR_GLOBAL, "pgs_bmi_chr22.rds"))

variant_counts <- data.frame(
  Step = c(
    "Raw GDS (chr22)",
    "MAF & missing %",
    "After LD pruning",
    "PGS003443 (T2D) chr22 all",
    "PGS003443 (T2D) Overlap",
    "PGS004994 (BMI) chr22 all",
    "PGS004994 (BMI) Overlap"
  ),
  Variants = c(
    ann$n_gds_raw,
    ann$n_maf_ok,
    ann$n_ld_keep,
    NA_integer_,                             # fill from files if present
    nrow(pgs_t2d_chr22),
    NA_integer_,
    nrow(pgs_bmi_chr22)
  ),
  stringsAsFactors = FALSE
)

# Try to fill "all" counts from 1a text if available (optional); otherwise leave NA
d <- stage_dirs("1_Phenotype")
f <- file.path(d$out_t, "pgs_overlap_counts.txt")
if (file.exists(f)) {
  lines <- readLines(f)
  get_num <- function(pat) { as.integer(gsub(".*: *", "", grep(pat, lines, value = TRUE))) }
  variant_counts$Variants[variant_counts$Step == "PGS003443 (T2D) chr22 all"] <- get_num("PGS003443 chr22 all")
  variant_counts$Variants[variant_counts$Step == "PGS004994 (BMI) chr22 all"] <- get_num("PGS004994 chr22 all")
}

write.csv(variant_counts, file.path(d$out_t, "variant_counts.csv"), row.names = FALSE)
cat("✓ Wrote:", file.path(d$out_t, "variant_counts.csv"), "\n")
