## -------- Project config (paths, options) --------
## Assumes working dir == "Data Privacy"
options(stringsAsFactors = FALSE)

# The project root is the "Data Privacy" directory itself
PROJECT_ROOT <- normalizePath(Sys.getenv("PROJECT_ROOT", "."), mustWork = FALSE)

# Shared data dirs
DIR_DATA   <- file.path(PROJECT_ROOT, "Data")
DIR_VCF    <- file.path(DIR_DATA,   "vcf")
DIR_GDS    <- file.path(DIR_DATA,   "gds")
DIR_GLOBAL <- file.path(DIR_DATA,   "global")
DIR_PGS    <- file.path(DIR_DATA,   "pgs")

# Stage output dirs (helpers)
stage_dirs <- function(stage_name) {
  base <- file.path(PROJECT_ROOT, stage_name)
  list(
    rcode = file.path(base, "R Code"),
    out_g = file.path(base, "output", "graphs"),
    out_t = file.path(base, "output", "text"),
    tmp   = file.path(base, "temp")
  )
}

# Ensure core dirs exist
dir.create(DIR_VCF,    showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_GDS,    showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_GLOBAL, showWarnings = FALSE, recursive = TRUE)

# CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- Safe PNG opener for headless clusters ----
png_open <- function(filename, width, height, res = 96) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = width, height = height, units = "px", res = res)
  } else {
    grDevices::png(filename, width = width, height = height,
                   type = if (capabilities("cairo")) "cairo" else "Xlib")
  }
}
