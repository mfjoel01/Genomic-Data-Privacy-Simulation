#!/usr/bin/env Rscript
# Stage 5 – Method comparison summary plots
# Produces:
#   - GWAS SNP effect heatmaps (Baseline vs DP for BMI; Baseline vs Federated for T2D)
#   - Predictive score distributions (violins) for Baseline / DP / Federated
#   - Predictive scores (pirate-style) with jitter + mean lines
#   - Per‑person PGS y=x scatter: Baseline vs DP
#   - Run‑to‑run variability (pirate-style, 100x) for Δβ
#
# This script reuses cached artifacts when available.
# If needed, it recomputes DP effect sizes for BMI (DP‑OLS) and T2D (DP logistic from noisy counts).

# Ensure packages & project config
source(file.path("0_Prep", "R Code", "0a_env_setup.R"), local = FALSE)
source(file.path("R", "config.R"))

d <- stage_dirs("5_MethodComparison")
dir.create(d$out_g, showWarnings = FALSE, recursive = TRUE)
dir.create(d$out_t, showWarnings = FALSE, recursive = TRUE)
dir.create(d$tmp,   showWarnings = FALSE, recursive = TRUE)

src <- function(f) source(file.path("5_MethodComparison", "R Code", f), local = FALSE)

cat("== Stage 5: Method Comparison Plots ==\n")

src("5a_compare_gwas_effects.R")
cat("== Part 5a: GWAS effect heatmaps done ==\n")

src("5b_predictive_scores_violin.R")
cat("== Part 5b: Predictive score violins done ==\n")

src("5c_predictive_scores_pirate.R")
cat("== Part 5c: Predictive score pirate plots done ==\n")

src("5d_pgs_per_person_scatter.R")
cat("== Part 5d: Per‑person y=x PGS scatter done ==\n")

src("5e_run_to_run_variability.R")
cat("== Part 5e: Run‑to‑run variability done ==\n")

cat("== Stage 5 complete ==\n")
