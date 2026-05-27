# PolygramLCR Replication Code

This directory contains the R package and scripts used to reproduce the
simulation and ATUS analyses.

## Contents

- `polygramLCR/`: local R package and numbered scripts for data generation,
  model fitting, result summaries, and figures.
- `ATUS/`: ATUS public-use files used by the real-data example.
- `outputs/`: generated data, model fits, summaries, and figures.

## Setup

Install the required R packages and then install the local package from this
repository's root directory.

```bash
Rscript -e 'install.packages(c("Rcpp", "RTriangle", "mgcv", "e1071"), repos = "https://cloud.r-project.org")'
R CMD INSTALL replication/polygramLCR
```

## Run the Analyses

Run the scripts from the package root, not from inside `scripts/`, because the
scripts use paths relative to `polygramLCR/`.

```bash
cd replication/polygramLCR

Rscript scripts/01_generate_simulation_data.R

Rscript scripts/02_1_fit_simulation_logistic.R
Rscript scripts/02_2_fit_simulation_gam.R
Rscript scripts/02_3_fit_simulation_svm.R
Rscript scripts/02_4_fit_simulation_polygram.R
Rscript scripts/03_make_simulation_results.R

Rscript scripts/04_1_fit_real_data_logistic.R
Rscript scripts/04_2_fit_real_data_gam.R
Rscript scripts/04_3_fit_real_data_svm.R
Rscript scripts/04_4_fit_real_data_polygram.R
Rscript scripts/05_make_real_data_results.R
```

Generated outputs are written under `replication/outputs/`.

The simulation uses 30 replicates. SIM00--SIM02 use 5000 training
observations and 5000 test observations; SIM03 uses 10000 training observations
and 5000 test observations.

## Package Check

After installation, this command fits a small polygram model:

```r
set.seed(1)
X <- matrix(runif(90), ncol = 3)
X <- X / rowSums(X)
Z <- polygramLCR::comp_to_ilr(X)
dat <- data.frame(z1 = Z[, 1], z2 = Z[, 2], y = rbinom(nrow(Z), 1, 0.5))
fit <- polygramLCR::polygram(
  dat,
  centers = 10,
  number_lambdas = 2,
  lambda_max = 0.2,
  max_admm_iter = 5,
  max_newton_iter = 10,
  require_admm_convergence = FALSE,
  require_newton_convergence = FALSE
)
head(polygramLCR::predict_polygram(fit, dat))
```
