# Genomic Data Privacy Simulation

This repository contains scripts and workflows for simulating genomic data and evaluating privacy-preserving techniques across multiple analysis stages. Pipelines are organized by numbered folders that roughly follow the order of execution, from data preparation through method comparison.

## Repository structure
- `0_Prep/` – utilities for preparing simulation inputs and shared resources.
- `00_HapGen2Simulation/` – configuration and helpers for running HapGen2-based genotype simulations.
- `1_Phenotype/` – scripts for generating phenotypes from the simulated genotypes.
- `2_Baseline/` – baseline association analyses used as reference points.
- `3_DifferentialPrivacy/` – experiments applying differential privacy mechanisms to the analyses.
- `4_FederatedLearning/` – materials related to federated learning workflows.
- `5_MethodComparison/` – notebooks and scripts comparing results across approaches.
- `Scripts/` and `R/` – shared helper scripts and R utilities used throughout the project.

## Data availability
No datasets are hosted in this repository. To reproduce the simulations you will need to provide your own input data and configure paths accordingly. More information and download links for needed data can be found at the end.

## Quick start
1. Clone the repository and create a working directory for intermediate outputs.
2. Configure environment variables (e.g., reference data paths) expected by the shell drivers in `Scripts/`. These scripts orchestrate the numbered stages and can be adapted to your local toolchain (R, HapGen2, etc.).
3. Review the configuration files within each stage directory and update them to point to your local input data.

## How to run
Review the numbered folders in order to understand the expected workflow—each stage contains notes or scripts that describe required dependencies. After tailoring configuration variables to your environment:

- To execute the full pipeline sequentially, use the consolidated driver:

  ```bash
  bash Scripts/run_stages.sh
  ```
- You can also specify to run only certain stages

  ```bash
  bash Scripts/run_stages.sh 2 3     # runs only stages 2 then 3
  ```
  
- Individual stages can also be run with the dedicated helpers, for example:

  ```bash
  bash Scripts/0_prep_driver.sh      # prepare shared inputs
  bash Scripts/1_pheno_driver.sh     # derive phenotypes
  bash Scripts/2_baseline_driver.sh  # perform baseline analyses
  ```

  Adjust each script's configuration variables before running to point to your data locations.

  ## Outputs (high‑level)

* **Plots**: Manhattan/QQ, ROC curves, binned BMI vs PGS, K‑sweep panels, effect/score comparisons, per‑person y=x, and run‑to‑run variability. See each stage’s `output/graphs/`.
* **Text**: AUC/R² summaries, Bonferroni thresholds, DP selection summaries, sweep metrics (CSV). Global caches live under `Data/global/`. 
 
