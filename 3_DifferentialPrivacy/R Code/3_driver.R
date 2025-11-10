#!/usr/bin/env Rscript

source(file.path("0_Prep", "R Code", "0a_env_setup.R"), local = FALSE)

cat("== Stage 3: Differential Privacy (shared R session) ==\n")
source(file.path("R", "config.R"))
d <- stage_dirs("3_DifferentialPrivacy")
dir.create(d$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d$out_t, showWarnings = FALSE, recursive = TRUE)
dir.create(d$tmp,   showWarnings = FALSE, recursive = TRUE)

source(file.path(d$rcode, "3a_dp_helpers.R"), local = FALSE)
source(file.path(d$rcode, "3b_dp_gwas.R"),    local = FALSE)
source(file.path(d$rcode, "3c_dp_plots.R"),   local = FALSE)
source(file.path(d$rcode, "3d_dp_pgs.R"),     local = FALSE)

if (tolower(Sys.getenv("RUN_EPS_SWEEP", "1")) %in% c("1","true","yes")) {
  source(file.path(d$rcode, "3e_epsilon_sweep.R"), local = TRUE)
}

cat("== Stage 3 done ==\n")
