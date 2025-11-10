#!/usr/bin/env Rscript
source(file.path("0_Prep", "R Code", "0a_env_setup.R"), local = FALSE)

cat("== Stage 1: Phenotype & PGS (shared R session) ==\n")

source(file.path("R", "config.R"))
d <- stage_dirs("1_Phenotype")
dir.create(d$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d$out_t, showWarnings = FALSE, recursive = TRUE)
dir.create(d$tmp,   showWarnings = FALSE, recursive = TRUE)

source(file.path(d$rcode, "1a_covariates_pgs.R"), local = FALSE)
cat("== Part 1a done ==\n")

source(file.path(d$rcode, "1b_pheno_sim.R"),      local = FALSE)
cat("== Part 1b Phenotype Simulation done ==\n")

source(file.path(d$rcode, "1c_table_correction.R"), local = FALSE)
cat("== Part 1c Table done ==\n")

cat("== Stage 1 done ==\n")
