#!/usr/bin/env Rscript
source(file.path("0_Prep", "R Code", "0a_env_setup.R"), local = FALSE)

cat("== Stage 2: Baseline (shared R session) ==\n")
source(file.path("R", "config.R"))
d <- stage_dirs("2_Baseline")
dir.create(d$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d$out_t, showWarnings = FALSE, recursive = TRUE)
dir.create(d$tmp,   showWarnings = FALSE, recursive = TRUE)

source(file.path(d$rcode, "2a_gwas.R"),         local = FALSE)
cat("== Part 2a GWAS done ==\n")

source(file.path(d$rcode, "2b_manhattan_qq.R"), local = FALSE)
cat("== Part 2b Manhattan and QQ Plot done ==\n")

source(file.path(d$rcode, "2c_pgs.R"),          local = FALSE)
cat("== Part 2c PGS done ==\n")

source(file.path(d$rcode, "2d_roc_scatter.R"),  local = FALSE)
cat("== Part 2d ROC and Scatter Plot done ==\n")

cat("== Stage 2 done ==\n")
