#!/usr/bin/env Rscript
cat("== Stage 0: Prep (shared R session) ==\n")

source(file.path("R", "config.R"))

d <- stage_dirs("0_Prep")
dir.create(d$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d$out_t, showWarnings = FALSE, recursive = TRUE)
dir.create(d$tmp,   showWarnings = FALSE, recursive = TRUE)

source(file.path(d$rcode, "0a_env_setup.R"), local = FALSE)
cat("== ENV setup done ==\n")

source(file.path(d$rcode, "0b_vcf_to_gds.R"), local = FALSE)
cat("== VCF/GDS done ==\n")

source(file.path(d$rcode, "0c_filter.R"),     local = FALSE)
cat("== Filter done ==\n")


cat("== Stage 0_Prep  done ==\n")
