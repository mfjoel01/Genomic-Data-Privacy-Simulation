source(file.path("R", "config.R"))

# Use a project-local library so jobs share the same installed pkgs
LIB_DIR <- file.path(PROJECT_ROOT, "R Libraries")
dir.create(LIB_DIR, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(LIB_DIR, .libPaths()))

pkgs_cran <- c(
  "data.table","pROC","qqman","future.apply","parallel","bit64","matrixStats",
  "parallelly","logistf","broom","zoo","ragg"
)
pkgs_bioc <- c("SNPRelate","gdsfmt","VariantAnnotation","vcfR","PhenotypeSimulator","Rsamtools")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install_and_load <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) BiocManager::install(pkg, ask = FALSE, update = FALSE)
      else      install.packages(pkg, dependencies = TRUE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

install_and_load(pkgs_cran, bioc = FALSE)
install_and_load(pkgs_bioc, bioc = TRUE)

if (isTRUE(capabilities("cairo"))) {
  options(bitmapType = "cairo")
}

cat("✓ Packages ready in:", LIB_DIR, "\n")
