#!/usr/bin/env Rscript
source(file.path("0_Prep", "R Code", "0a_env_setup.R"), local = FALSE)

cat("== Stage 4: Federated Learning (shared R session) ==\n")
source(file.path("R", "config.R"))
d <- stage_dirs("4_FederatedLearning")
dir.create(d$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d$out_t, showWarnings = FALSE, recursive = TRUE)
dir.create(d$tmp,   showWarnings = FALSE, recursive = TRUE)

source(file.path(d$rcode, "4a_fed_utils.R"), local = FALSE)
source(file.path(d$rcode, "4b_fed_gwas.R"),  local = FALSE)
source(file.path(d$rcode, "4c_fed_pgs.R"),   local = FALSE)
source(file.path(d$rcode, "4d_fed_plots.R"), local = FALSE)

# NEW: run K-sweep by default; set RUN_FED_SWEEP=0|false|no to skip
run_sweep <- !(tolower(Sys.getenv("RUN_FED_SWEEP", "1")) %in% c("0","false","no"))
if (run_sweep) {
  source(file.path(d$rcode, "4e_fed_k_sweep.R"), local = FALSE)
}

cat("== Stage 4 done ==\n")
